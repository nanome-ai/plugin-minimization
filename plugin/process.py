import time
import nanome
from nanome.util import Logs, Octree, Process
from nanome.api.structure import Complex
from nanome.util.stream import StreamCreationError
from nanome.util.enums import StreamType

import tempfile
from collections import deque
from functools import partial
import os
import sys

PACKET_QUEUE_LEN = 20

IS_WIN = sys.platform.startswith('win')
SDFOPTIONS = Complex.io.SDFSaveOptions()
SDFOPTIONS.write_bonds = True

# hack to fix convert_to_frames killing atom indices:
_atom_shallow_copy = nanome._internal.structure._Atom._shallow_copy


def _atom_shallow_copy_fix(self, *args):
    atom = _atom_shallow_copy(self, *args)
    atom._index = self._index
    return atom


nanome._internal.structure._Atom._shallow_copy = _atom_shallow_copy_fix


class MinimizationProcess():
    def __init__(self, plugin, nanobabel_dir):
        self.__plugin = plugin
        self.is_running = False
        self.__process_running = False
        self.__stream = None
        self.__data_queue = []
        self.__nanobabel_dir = nanobabel_dir
        self.temp_dir = tempfile.TemporaryDirectory()

    async def start_process(self, workspace, ff, steps, steepest):
        if sum(1 for _ in workspace.complexes) == 0:
            Logs.message('No structures to minimize')
            return
        input_file = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)
        constraints_file = tempfile.NamedTemporaryFile(delete=False, suffix='.txt', dir=self.temp_dir.name)
        output_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb', dir=self.temp_dir.name)
        self.__output_lines = []
        self.__updates_done = {}
        self.__packet_id = 0

        (saved_atoms, indices) = self.__save__atoms(input_file.name, workspace)
        Logs.debug("Wrote input file:", input_file.name)
        self.__save__constraints(constraints_file.name, saved_atoms)
        Logs.debug("Wrote constraints file:", constraints_file.name)
        self.__stream, error = await self.__plugin.create_writing_stream(indices, StreamType.position)

        if error == StreamCreationError.AtomNotFound:
            # User deleted atom in time between start_process() and create_writing_stream().
            # so lets update the workspace and try again
            Logs.warning(f"User deleted atoms while setting up process, retrying")
            updated_workspace = await self.__plugin.request_workspace()
            await self.start_process(updated_workspace, ff, steps, steepest)
            return

        elif error != StreamCreationError.NoError:
            Logs.error(f"Error while creating stream: {error}")
            return

        self.__data_queue = deque()
        cwd_path = self.__nanobabel_dir
        exe = 'nanobabel.exe' if IS_WIN else 'nanobabel'
        exe_path = os.path.join(cwd_path, exe)
        args = ['minimize', '-h', '-l', '20', '-n', str(steps), '-ff', ff, '-i', input_file.name, '-cx', constraints_file.name, '-o', output_file.name]
        if IS_WIN:
            args += ['-dd', 'data']
        if steepest:
            args.append('-sd')
        Logs.debug(args)

        p = Process(exe_path, args, True)
        p.on_error = self.__on_process_error
        p.on_output = self.__on_process_output
        p.on_done = self.__on_process_done
        self.calculation_start_time = time.time()
        log_data = {
            'exe_path': exe_path,
            'force_field': ff,
            'steps': steps,
            'steepest': steepest
        }
        Logs.message("Starting Minimization Process", extra=log_data)
        p.start()

        self.__process = p
        self.__process_running = True
        self.is_running = True

    def stop_process(self):
        if self.__process_running:
            self.__process.stop()
        if self.__stream is not None:
            self.__stream.destroy()
            self.__stream = None
        self.is_running = False
        self.__plugin.minimization_done()

    def update(self):
        if not self.is_running or \
                self.__packet_id > PACKET_QUEUE_LEN and \
                not self.__updates_done[self.__packet_id - PACKET_QUEUE_LEN]:
            return

        if len(self.__data_queue) > 0:
            data_chunk = self.__data_queue.popleft()
            complex = nanome.api.structure.Complex.io.from_pdb(lines=data_chunk)
            self.__match_and_move(complex)
        elif not self.__process_running:
            Logs.debug('Minimization complete')
            self.stop_process()

    def __on_process_error(self, error):
        Logs.warning('Error in nanobabel process:')
        Logs.warning(error)

    def __on_process_output(self, output):
        output = output.strip()
        if output != '':
            split_output = output.split('\n')
            self.__processing_output(split_output)

    def __on_process_done(self, code):
        self.__process_running = False

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
        self.__updates_done[self.__packet_id] = False
        self.__stream.update(positions, partial(self.__update_done, self.__packet_id))
        self.__packet_id += 1

    def __update_done(self, packet_id):
        self.__updates_done[packet_id] = True

    def __processing_output(self, split_output):
        for line in split_output:
            if "Step update start" in line:
                self.__output_lines.clear()
            elif "Step update end" in line:
                self.__data_queue.append(self.__output_lines.copy())
            else:
                self.__output_lines.append(line)

    def __save__atoms(self, path, workspace):
        visible_complexes = []
        for complex in workspace.complexes:
            if not complex.visible:
                continue

            current_frame = complex.current_frame
            if len(list(complex.molecules)) == 1:
                current_frame = next(complex.molecules).current_conformer

            complex = complex.convert_to_frames()
            complex.current_frame = current_frame
            visible_complexes.append(complex)

        selected_atoms = Octree()

        for complex in visible_complexes:
            complex_local_to_workspace_matrix = complex.get_complex_to_workspace_matrix()
            molecule = complex._molecules[complex.current_frame]

            for atom in molecule.atoms:
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

        for complex in visible_complexes:
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
