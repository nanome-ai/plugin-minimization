import os
import sys
import nanome
from nanome.util import Logs, async_callback
from nanome.util.enums import Integrations

from .menu import MinimizationMenu
from .process import MinimizationProcess

NANOBABEL = os.environ.get('NANOBABEL', os.path.join(os.getcwd(), 'nanobabel'))
if not os.path.exists(NANOBABEL):
    NANOBABEL = None


class Minimization(nanome.AsyncPluginInstance):

    def start(self):
        self.__menu = MinimizationMenu(self)
        self._process = MinimizationProcess(self, NANOBABEL)
        self.__menu.build_menu()
        self.__integration_request = None
        self.integration.minimization_start = self.start_integration
        self.integration.minimization_stop = self.stop_integration

    @async_callback
    async def start_integration(self, request):
        (ff, steps, steepest, cutoff) = request.get_args()
        if self._process.is_running == True:
            if self.__integration_request != None:
                self.__integration_request.send_response(False)
            self._process.stop_process()
        ff = self.convert_forcefield_value(ff)

        self.__menu.change_running_status(True)
        await self.start_minimization(ff, steps, steepest)

    def stop_integration(self, request):
        self._process.stop_process()
        request.send_response(None)
        if self.__integration_request != None:
            self.__integration_request.send_response(True)
        self.__menu.change_running_status(False)

    def update(self):
        self._process.update()

    def on_run(self):
        self.__menu.toggle_minimization()

    def on_advanced_settings(self):
        self.__menu.open_menu()

    def on_stop(self):
        self.stop_minimization()

    async def start_minimization(self, ff, steps, steepest):
        ff = self.convert_forcefield_value(ff)
        workspace = await self.request_workspace()
        if sum(1 for _ in workspace.complexes) > 0:
            await self._process.start_process(workspace, ff, steps, steepest)
        else:
            Logs.message("No complexes found. nothing to minimize.")

    def stop_minimization(self):
        self._process.stop_process()

    def minimization_done(self):
        self.__menu.change_running_status(False)
        if self.__integration_request != None:
            self.__integration_request.send_response(True)

    def set_run_status(self, running):
        btn_type = nanome.util.enums.PluginListButtonType.run
        if running:
            self.set_plugin_list_button(btn_type, "Stop")
        else:
            self.set_plugin_list_button(btn_type, "Run")

    def convert_forcefield_value(self, value):
        if value == "General Amber" or value == 1:
            return 'Gaff'
        elif value == "Ghemical" or value == 2:
            return 'Ghemical'
        elif value == "MMFF94" or value == 3:
            return 'MMFF94'
        elif value == "MMFF94s" or value == 4:
            return 'MMFF94s'
        elif value == "Universal" or value == 0:
            return 'Uff'
        return 'Uff'


def main():
    if not NANOBABEL:
        Logs.error('Error: nanobabel not found, please set NANOBABEL env var')
        sys.exit(1)

    plugin = nanome.Plugin("Minimization", "Run minimization on selected structures. See Advanced Parameters for forcefield, number of steps, and steepest descent", "Minimization", True, integrations=[Integrations.minimization])
    plugin.set_plugin_class(Minimization)
    plugin.run()


if __name__ == "__main__":
    main()
