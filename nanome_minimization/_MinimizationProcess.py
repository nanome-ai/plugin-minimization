import nanome
from nanome.util import Logs, Octree
from nanome.api.structure import Complex
from nanome.util.stream import StreamCreationError
from nanome.util.enums import StreamType

import tempfile
import subprocess
from timeit import default_timer as timer
import traceback
from collections import deque
import os

SDFOPTIONS = Complex.io.SDFSaveOptions()
SDFOPTIONS.write_bonds = True

# hack to fix convert_to_frames killing atom indices:
_atom_shallow_copy = nanome._internal._structure._Atom._shallow_copy
def _atom_shallow_copy_fix(self, *args):
    atom = _atom_shallow_copy(self, *args)
    atom._index = self._index
    return atom
nanome._internal._structure._Atom._shallow_copy = _atom_shallow_copy_fix


class MinimizationProcess():
    def __init__(self, plugin, nanobabel_dir):
        self.__plugin = plugin
        self._is_running = False
        self.__stream = None
        self.__nanobabel_dir = nanobabel_dir

    def start_process(self, workspace, ff, steps, steepest):
        def on_stream_creation(stream, error):
            if error != StreamCreationError.NoError:
                Logs.error("Error while creating stream")
                return

            self.__stream = stream
            self.__data_queue = deque()
            cwd_path = self.__nanobabel_dir
            exe_path = os.path.join(cwd_path, 'nanobabel.exe')
            args = [exe_path, 'MINIMIZE', '-h', '-l', '1', '-n', str(steps), '-ff', ff, '-i', input_file.name, '-cx', constraints_file.name, '-o', output_file.name, '-dd', 'data']
            Logs.debug(args)
            if steepest:
                args.append('-sd')
            try:
                self.__process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=0, cwd=cwd_path, universal_newlines=True, encoding="utf-8")
                self._is_running = True
                Logs.debug("Nanobabel started")
            except:
                Logs.error("Couldn't execute nanobabel, please check if executable is in the plugin folder and has permissions:\n", traceback.format_exc())

        input_file = tempfile.NamedTemporaryFile(delete=False, suffix='.mol')
        constraints_file = tempfile.NamedTemporaryFile(delete=False, suffix='.txt')
        output_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
        self.__output_lines = []
        self.__timer = timer()
        self.__updating_stream = False
        self.__processed_output = 0
        self.__processed_error = 0

        (saved_atoms, indices) = self.__save__atoms(input_file.name, workspace)
        Logs.debug("Wrote input file:", input_file.name)
        self.__save__constraints(constraints_file.name, saved_atoms)
        Logs.debug("Wrote constraints file:", constraints_file.name)
        self.__stream = self.__plugin.create_writing_stream(indices, StreamType.position, on_stream_creation)

    def stop_process(self):
        self._is_running = False
        if self.__stream is not None:
            self.__stream.destroy()
        self.__stream = None
        self.__plugin.minimization_done()

    def update(self):
        if self._is_running == False:
            return

        output, error = self.__process.communicate()

        error = error[self.__processed_error:]
        self.__processed_error += len(error)
        error = error.strip()
        if error != '':
            Logs.error(error)  # Maybe abort in case of error?

        output = output[self.__processed_output:]
        self.__processed_output += len(output)
        output = output.strip()
        if output != '':
            split_output = output.split('\n')
            self.__processing_output(split_output)

        if len(self.__data_queue) > 0 and self.__updating_stream == False:
            data_chunk = self.__data_queue.popleft()
            complex = nanome.api.structure.Complex.io.from_pdb(lines=data_chunk)
            # if timer() - self.__timer >= 0.1:
            self.__match_and_move(complex)
            # self.__timer = timer()

        if self.__process.poll() is not None and len(self.__data_queue) == 0:
            self.stop_process()
            Logs.debug("Nanobabel done")
            return

    def __match_and_move(self, complex):
        positions = [0] * (len(self.__atom_position_index_by_serial) * 3)
        for atom in complex.atoms:
            if atom.serial in self.__atom_position_index_by_serial:
                (idx, complex) = self.__atom_position_index_by_serial[atom.serial]
                complex_absolute_to_relative = complex.get_workspace_to_complex_matrix()
                atom_relative_pos = complex_absolute_to_relative * atom.position
                positions[idx * 3] = atom_relative_pos.x
                positions[idx * 3 + 1] = atom_relative_pos.y
                positions[idx * 3 + 2] = atom_relative_pos.z
        if self.__stream == None:
            return
        self.__updating_stream = True
        self.__stream.update(positions, self.__update_done)

    def __update_done(self):
        self.__updating_stream = False

    def __processing_output(self, split_output):
        for line in split_output:
            if "Step update start" in line:
                self.__output_lines.clear()
            elif "Step update end" in line:
                self.__data_queue.append(self.__output_lines.copy())
            else:
                self.__output_lines.append(line)

    def __save__atoms(self, path, workspace):
        selected_atoms = Octree()

        for complex in workspace.complexes:
            complex_local_to_workspace_matrix = complex.get_complex_to_workspace_matrix()
            for atom in complex.atoms:
                if atom.selected is True:
                    atom_absolute_pos = complex_local_to_workspace_matrix * atom.position
                    selected_atoms.add(atom, atom_absolute_pos)

        self.__atom_position_index_by_serial = dict()
        atom_by_index = dict()
        saved_atoms = []
        saved_atoms_indices = []
        atom_position_index = 0
        serial = 1
        found_atoms = []
        bonds = []
        result_complex = nanome.structure.Complex()
        result_molecule = nanome.structure.Molecule()
        result_chain = nanome.structure.Chain()
        result_residue = nanome.structure.Residue()
        result_complex.add_molecule(result_molecule)
        result_molecule.add_chain(result_chain)
        result_chain.add_residue(result_residue)

        for complex in workspace.complexes:
            if not complex.visible:
                continue

            complex = complex.convert_to_frames()
            complex_local_to_workspace_matrix = complex.get_complex_to_workspace_matrix()

            molecule = complex._molecules[complex.current_frame]
            for atom in molecule.atoms:
                atom_absolute_pos = complex_local_to_workspace_matrix * atom.position
                selected_atoms.get_near_append(atom_absolute_pos, 7, found_atoms, 1)

                if len(found_atoms) > 0 and atom not in saved_atoms:
                    atom.serial = serial
                    serial += 1
                    self.__atom_position_index_by_serial[atom.serial] = (atom_position_index, complex)
                    atom_position_index += 1
                    atom.position = atom_absolute_pos
                    saved_atoms.append(atom)
                    saved_atoms_indices.append(atom.index)
                    result_residue.add_atom(atom)
                    atom_by_index[atom.index] = atom

                    for bond in atom.bonds:
                        if bond not in bonds \
                            and bond.atom1.index in atom_by_index \
                                and bond.atom2.index in atom_by_index:
                            bonds.append(bond)

                found_atoms.clear()

        for bond in bonds:
            result_residue.add_bond(bond)

        result_complex.io.to_sdf(path, SDFOPTIONS)
        return (saved_atoms, saved_atoms_indices)

    def __save__constraints(self, path, saved_atoms):
        file = open(path, 'w')
        for atom in saved_atoms:
            if atom.selected is False:
                file.write("ATOM:FIXED:" + str(atom.serial) + "\n")
        file.close()
