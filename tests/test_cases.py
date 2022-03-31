import asyncio
import os
import unittest
from unittest.mock import patch
from random import randint

from unittest.mock import MagicMock
from nanome.api.structure import Chain, Complex, Molecule, Workspace
from nanome.util.stream import StreamCreationError
from plugin.Minimization import Minimization


fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')


def run_awaitable(awaitable, *args, **kwargs):
    loop = asyncio.get_event_loop()
    if loop.is_running:
        loop = asyncio.new_event_loop()
    loop.run_until_complete(awaitable(*args, **kwargs))
    loop.close()


class MinimizationTestCase(unittest.TestCase):
    """Test different combinations of args for calculate_interactions."""

    def setUp(self):
        tyl_pdb = f'{fixtures_dir}/1tyl.pdb'
        self.complex = Complex.io.from_pdb(path=tyl_pdb)
        self.workspace = Workspace()
        for atom in self.complex.atoms:
            atom.index = randint(1000000000, 9999999999)

        self.plugin_instance = Minimization()
        self.plugin_instance.start()
        self.plugin_instance._network = MagicMock()

        # Split out ligand into separate complex
        target_complex = self.complex
        chain_name = 'HC'
        residue_name = 'TYL'
        ligand_residue = next(res for res in self.complex.residues if res.name == residue_name)

        # Build new complex containing ligand residue
        ligand_complex = Complex()
        ligand_molecule = Molecule()
        ligand_chain = Chain()
        ligand_chain.name = chain_name

        ligand_chain.add_residue(ligand_residue)
        ligand_molecule.add_chain(ligand_chain)
        ligand_complex.add_molecule(ligand_molecule)

        target_complex.index = 98
        ligand_complex.index = 99
        self.workspace.add_complex(self.complex)
        self.workspace.add_complex(ligand_complex)

    @patch('nanome._internal._network._ProcessNetwork._instance')
    @patch('nanome.api.plugin_instance.PluginInstance.create_writing_stream')
    @patch('nanome.api.plugin_instance.PluginInstance.request_workspace')
    def test_start_minimization(self, request_workspace_mock, create_writing_stream_mock, mock_network):
        """Validate calculate_interactions call where ligand is on a separate complex."""
        ff = 'Uff'
        steps = 100
        steepest = True
        # Set up mocked result for create_writing_stream_mock
        stream_fut = asyncio.Future()
        stream_fut.set_result((MagicMock(), StreamCreationError.NoError))
        create_writing_stream_mock.return_value = stream_fut

        ws_fut = asyncio.Future()
        ws_fut.set_result(self.workspace)
        request_workspace_mock.return_value = ws_fut

        return run_awaitable(self.validate_start_minimization, ff, steps, steepest)

    async def validate_start_minimization(
            self, ff, steps, steepest):
        """Run plugin.calculate_interactions with provided args and make sure lines are added to LineManager."""
        await self.plugin_instance.start_minimization(ff, steps, steepest)
