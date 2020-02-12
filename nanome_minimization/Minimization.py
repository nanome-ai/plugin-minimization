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

    def update(self):
        self._process.update()

    def on_run(self):
        self.__menu.toggle_minimization()

    def on_advanced_settings(self):
        self.__menu.open_menu()

    def start_minimization(self, ff, steps, steepest):
        def workspace_received(workspace):
            self._process.start_process(workspace, ff, steps, steepest)

        self.request_workspace(workspace_received)

    def stop_minimization(self):
        self._process.stop_process()

    def minimization_done(self):
        self.__menu.change_running_status(False)

    def set_run_status(self, running):
        if running:
            self.set_plugin_list_button(nanome.PluginInstance.PluginListButtonType.run, "Stop")
        else:
            self.set_plugin_list_button(nanome.PluginInstance.PluginListButtonType.run, "Run")

def main():
    if not NANOBABEL:
        Logs.error('Error: nanobabel not found, please set NANOBABEL env var')
        sys.exit(1)

    plugin = nanome.Plugin("Minimization", "Run minimization on selected structures. See Advanced Parameters for forcefield, number of steps, and steepest descent", "Minimization", True)
    plugin.set_plugin_class(Minimization)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()
