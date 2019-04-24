import nanome
from nanome.util import Logs, Octree

import tempfile
import subprocess
from timeit import default_timer as timer

class MinimizationProcess():
    def __init__(self, plugin):
        self.__plugin = plugin
        self.__is_running = False
        self.__stream = None
    
    def start_process(self, workspace, ff, steps, steepest):
        def on_stream_creation(stream):
            self.__stream = stream
            args = ['./nanobabel.exe', 'MINIMIZE', '-h', '-l', '1', '-n', steps, '-ff', ff, '-i', input_file.name, '-cx', constraints_file.name, '-o', output_file.name, '-dd', data_dir.name]
            if steepest:
                args.append('-sd')
            try:
                self.__process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.__is_running = True
            except:
                Logs.error("Couldn't execute nanobabel, please check if executable is in the plugin folder and has permissions")

        input_file = tempfile.NamedTemporaryFile(delete=False)
        constraints_file = tempfile.NamedTemporaryFile(delete=False)
        output_file = tempfile.NamedTemporaryFile(delete=False)
        data_dir = tempfile.TemporaryDirectory()
        self.__output_lines = []
        self.__timer = timer()

        self.__workspace_absolute_to_relative_matrix = workspace.transform.get_absolute_to_relative_matrix()
        self.__workspace_relative_to_absolute_matrix = workspace.transform.get_relative_to_absolute_matrix()
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

        if self.__process.poll() == None:
            self.stop_process()
            return

        output, errors = self.__process.communicate()
        for error in errors:
            Logs.error(error) # Maybe abort in case of error?

        data_chunk = self.__processing_output(output)
        if data_chunk != None:
            # Using internal functions here, we should expose a better API here
            content = nanome._internal._structure._io._pdb.parse_string(data_chunk)
            complex = nanome._internal._structure._io._pdb.structure(content)
            if timer() - self.__timer >= 0.1:
                self.__match_and_move(complex)
                self.__timer = timer()

    @staticmethod
    def __generate_key(molecule_number, atom_serial):
        molecule = molecule_number << 32
        return molecule | atom_serial

    def __match_and_move(self, complex):
        positions = [0] * len(self.__atom_position_index_by_serial) * 3
        for atom in complex.atoms:
            if atom.serial in self.__atom_position_index_by_serial:
                (idx, complex) = self.__atom_position_index_by_serial[atom.serial]
                complex_absolute_to_relative = complex.transform.get_absolute_to_relative_matrix()
                atom_absolute_pos = self.__workspace_relative_to_absolute_matrix * complex_absolute_to_relative * atom.transform.position
                positions[idx] = atom_absolute_pos.x
                positions[idx + 1] = atom_absolute_pos.y
                positions[idx + 2] = atom_absolute_pos.z
        self.__stream.update(positions)

    def __processing_output(self, output):
        result = None
        for line in output:
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
            complex_local_to_world_matrix = complex.transform.get_relative_to_absolute_matrix()
            for atom in complex.atoms:
                if atom.selected == True:
                    atom_absolute_pos = complex_local_to_world_matrix * self.__workspace_absolute_to_relative_matrix * atom.position
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
            complex_local_to_world_matrix = complex.get_relative_to_absolute_matrix()

            molecule = complex.molecules[complex.rendering.current_frame]
            for atom in molecule.atoms:
                atom_absolute_pos = atom.position * complex_local_to_world_matrix * self.__workspace_absolute_to_relative_matrix
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
            if atom.selected == False:
                file.write("ATOM:FIXED:" + atom.molecular.serial)
        file.close()