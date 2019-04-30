import nanome
from nanome.util import Logs, Octree
from nanome.util.stream import StreamCreationError

import tempfile
import subprocess
from timeit import default_timer as timer
import traceback

class MinimizationProcess():
    def __init__(self, plugin):
        self.__plugin = plugin
        self.__is_running = False
        self.__stream = None
    
    def start_process(self, workspace, ff, steps, steepest):
        def on_stream_creation(stream, error):
            if error != StreamCreationError.NoError:
                Logs.error("Error while creating stream")
                return

            self.__stream = stream
            args = ['./nanobabel/nanobabel.exe', 'MINIMIZE', '-h', '-l', '1', '-n', str(steps), '-ff', ff, '-i', input_file.name, '-cx', constraints_file.name, '-o', output_file.name, '-dd', 'data']
            Logs.debug(args)
            if steepest:
                args.append('-sd')
            try:
                self.__process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='nanobabel', text=True, encoding="utf-8")
                self.__is_running = True
                Logs.debug("Nanobabel started")
            except:
                Logs.error("Couldn't execute nanobabel, please check if executable is in the plugin folder and has permissions:\n", traceback.format_exc())

        input_file = tempfile.NamedTemporaryFile(delete=False, suffix='.mol')
        constraints_file = tempfile.NamedTemporaryFile(delete=False, suffix='.txt')
        output_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
        self.__output_lines = []
        self.__timer = timer()

        (saved_atoms, indices) = self.__save__atoms(input_file.name, workspace)
        Logs.debug("Wrote input file:", input_file.name)
        self.__save__constraints(constraints_file.name, saved_atoms)
        Logs.debug("Wrote constraints file:", constraints_file.name)
        self.__stream = self.__plugin.create_stream(indices, on_stream_creation)

    def stop_process(self):
        self.__is_running = False
        if self.__stream != None:
            self.__stream.destroy()
        self.__plugin.minimization_done()

    def update(self):
        if self.__is_running == False:
            return

        output, error = self.__process.communicate()
        if error != None:
            Logs.error(error) # Maybe abort in case of error?

        split_output = output.split('\n')
        data_chunk = self.__processing_output(split_output)
        if data_chunk != None:
            # Using internal functions here, we should expose a better API here
            content = nanome._internal._structure._io._pdb.parse_string(data_chunk)
            complex = nanome._internal._structure._io._pdb.structure(content)
            if timer() - self.__timer >= 0.1:
                self.__match_and_move(complex)
                self.__timer = timer()

        if self.__process.poll() == None:
            self.stop_process()
            Logs.debug("Nanobabel done")
            return

    @staticmethod
    def __generate_key(molecule_number, atom_serial):
        molecule = molecule_number << 32
        return molecule | atom_serial

    def __match_and_move(self, complex):
        positions = [0] * len(self.__atom_position_index_by_serial) * 3
        for atom in complex.atoms:
            if atom.molecular.serial in self.__atom_position_index_by_serial:
                (idx, complex) = self.__atom_position_index_by_serial[atom.molecular.serial]
                complex_absolute_to_relative = complex.transform.get_absolute_to_relative_matrix()
                atom_absolute_pos = complex_absolute_to_relative * atom.molecular.position
                positions[idx] = atom_absolute_pos.x
                positions[idx + 1] = atom_absolute_pos.y
                positions[idx + 2] = atom_absolute_pos.z
        self.__stream.update(positions)

    def __processing_output(self, split_output):
        result = None
        for line in split_output:
            if "Step update start" in line:
                self.__output_lines.clear()
            elif "Step update end" in line:
                result = self.__output_lines.copy()
            else:
                self.__output_lines.append(line)
        return result

    def __save__atoms(self, path, workspace):
        selected_atoms = Octree()

        for complex in workspace.complexes:
            complex_local_to_workspace_matrix = complex.transform.get_relative_to_absolute_matrix()
            for atom in complex.atoms:
                if atom.rendering.selected == True:
                    atom_absolute_pos = complex_local_to_workspace_matrix * atom.molecular.position
                    selected_atoms.add(atom, atom_absolute_pos)

        self.__atom_position_index_by_serial = dict()
        atom_by_index = dict()
        saved_atoms = []
        saved_atoms_indices = []
        atom_position_index = 0
        found_atoms = []
        result_complex = nanome.structure.Complex()
        result_molecule = nanome.structure.Molecule()
        result_chain = nanome.structure.Chain()
        result_residue = nanome.structure.Residue()
        result_complex.add_molecule(result_molecule)
        result_molecule.add_chain(result_chain)
        result_chain.add_residue(result_residue)

        for complex in workspace.complexes:
            complex_local_to_workspace_matrix = complex.transform.get_relative_to_absolute_matrix()

            molecule = complex._molecules[complex.rendering.current_frame]
            for atom in molecule.atoms:
                atom_absolute_pos = complex_local_to_workspace_matrix * atom.molecular.position
                selected_atoms.get_near_append(atom_absolute_pos, 7, found_atoms, 1)

                if len(found_atoms) > 0:
                    self.__atom_position_index_by_serial[atom.molecular.serial] = (atom_position_index, complex)
                    atom_position_index += 1
                    saved_atoms.append(atom)
                    saved_atoms_indices.append(atom.index)
                    result_residue.add_atom(atom)
                    atom_by_index[atom.index] = atom

                    for bond in atom.bonds:
                        if bond.atom1.index in atom_by_index and bond.atom2.index in atom_by_index:
                            result_residue.add_bond(bond)

                found_atoms.clear()

        result_complex.io.to_sdf(path)
        return (saved_atoms, saved_atoms_indices)

    def __save__constraints(self, path, saved_atoms):
        file = open(path, 'w')
        for atom in saved_atoms:
            if atom.rendering.selected == False:
                file.write("ATOM:FIXED:" + str(atom.molecular.serial))
        file.close()