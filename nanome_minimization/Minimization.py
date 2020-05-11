import os
import sys
import nanome
from nanome.util import Logs

from ._MinimizationMenu import MinimizationMenu
from ._MinimizationProcess import MinimizationProcess

NANOBABEL = os.environ.get('NANOBABEL', os.path.join(os.getcwd(), 'nanobabel'))
if not os.path.exists(NANOBABEL):
    NANOBABEL = None

class Minimization(nanome.PluginInstance):
    def start(self):
        self.__menu = MinimizationMenu(self)
        self._process = MinimizationProcess(self, NANOBABEL)
        self.__menu.build_menu()
        self.__integration_request = None
        self.integration.minimization_start = self.start_integration
        self.integration.minimization_stop = self.stop_integration

    def start_integration(self, request):
        (ff, steps, steepest, cutoff) = request.get_args()
        if self._process._is_running == True:
            if self.__integration_request != None:
                self.__integration_request.send_response(False)
            self._process.stop_process()
        ff = self.convert_forcefield_value(ff)

        def workspace_received(workspace):
            self._process.start_process(workspace, ff, steps, steepest)
        self.request_workspace(workspace_received)
        self.__menu.change_running_status(True)

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

    def start_minimization(self, ff, steps, steepest):
        ff = self.convert_forcefield_value(ff)
        def workspace_received(workspace):
            self._process.start_process(workspace, ff, steps, steepest)

        self.request_workspace(workspace_received)

    def stop_minimization(self):
        self._process.stop_process()

    def minimization_done(self):
        self.__menu.change_running_status(False)
        if self.__integration_request != None:
            self.__integration_request.send_response(True)

    def set_run_status(self, running):
        if running:
            self.set_plugin_list_button(nanome.PluginInstance.PluginListButtonType.run, "Stop")
        else:
            self.set_plugin_list_button(nanome.PluginInstance.PluginListButtonType.run, "Run")

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

    plugin = nanome.Plugin("Minimization", "Run minimization on selected structures. See Advanced Parameters for forcefield, number of steps, and steepest descent", "Minimization", True)
    plugin.set_plugin_class(Minimization)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()
