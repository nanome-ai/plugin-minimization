import nanome

from _MinimizationMenu import MinimizationMenu
from _MinimizationProcess import MinimizationProcess

class Minimization(nanome.PluginInstance):
    def start(self):
        self.__menu = MinimizationMenu(self)
        self.__process = MinimizationProcess(self)
        self.__menu.build_menu()

    def update(self):
        self.__process.update()

    def on_run(self):
        self.__menu.start_minimization()

    def on_advanced_settings(self):
        self.__menu.open_menu()

    def start_minimization(self, ff, steps, steepest):
        def workspace_received(workspace):
            self.__process.start_process(workspace, ff, steps, steepest)

        self.request_workspace(workspace_received)

    def stop_minimization(self):
        self.__process.stop_process()

    def minimization_done(self):
        self.__menu.change_running_status(False)

if __name__ == "__main__":
    plugin = nanome.Plugin("Minimization", "Run minimization on selected structures. See Advanced Parameters for forcefield, number of steps, and steepest descent", "Minimization", True)
    plugin.set_plugin_class(Minimization)
    plugin.run('127.0.0.1', 8888)