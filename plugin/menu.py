import nanome
from nanome.util import async_callback
import os


class MinimizationMenu():
    def __init__(self, plugin):
        self.__plugin = plugin
        self.__menu = None
        self.__selected_ff_btn = None
        self.__nb_steps = 2500
        self.__steepest_descent = True
        self.__steepest_descent_btn = None
        self.__conjugate_gradient_btn = None
        self.__steps_label = None
        self.__start_btn = None
        self.__running = False

    def __update_start_btn(self, running):
        self.__start_btn.selected = running
        self.__plugin.update_content(self.__start_btn)
        self.__plugin.set_run_status(running)

    def change_running_status(self, running):
        self.__running = running
        self.__update_start_btn(running)

    def toggle_minimization(self):
        if self.__plugin._process.is_running:
            self.stop_minimization()
        else:
            self.start_minimization()

    @async_callback
    async def start_minimization(self):
        self.change_running_status(True)
        await self.__plugin.start_minimization(self.__get_selected_forcefield(), self.__nb_steps, self.__steepest_descent)

    def stop_minimization(self):
        self.change_running_status(False)
        self.__plugin.stop_minimization()

    def open_menu(self):
        self.__menu.enabled = True
        self.__plugin.update_menu(self.__menu)

    def build_menu(self):
        def ff_selected(btn):
            if self.__selected_ff_btn != None:
                self.__selected_ff_btn.selected = False
                self.__plugin.update_content(self.__selected_ff_btn)
            self.__selected_ff_btn = btn
            btn.selected = True
            self.__plugin.update_content(btn)

        def change_steps(btn):
            if btn.text.value_idle == "+":
                self.__nb_steps += 500
                self.__nb_steps = min(self.__nb_steps, 5000)
            elif btn.text.value_idle == "-":
                self.__nb_steps -= 500
                self.__nb_steps = max(self.__nb_steps, 500)
            self.__steps_label.text_value = str(self.__nb_steps)
            self.__plugin.update_content(self.__steps_label)

        def switch_steepest(btn):
            self.__steepest_descent = not self.__steepest_descent
            btn.selected = self.__steepest_descent
            self.__plugin.update_content(btn)
            self.__conjugate_gradient_btn.selected = not self.__steepest_descent
            self.__plugin.update_content(self.__conjugate_gradient_btn)

        def switch_conjugate_gradient(btn):
            self.__steepest_descent = not self.__steepest_descent
            btn.selected = not self.__steepest_descent
            self.__plugin.update_content(btn)
            self.__steepest_descent_btn.selected = self.__steepest_descent
            self.__plugin.update_content(self.__steepest_descent_btn)

        def switch_minimization(btn):
            self.toggle_minimization()

        # loading menus
        menu = nanome.ui.Menu.io.from_json(os.path.join(os.path.dirname(__file__), "minimization_menu.json"))
        self.__menu = menu
        self.__plugin.menu = menu

        # getting elements, for future update
        self.__steps_label = menu.root.find_node("steps_label", True).get_content()
        self.__start_btn = menu.root.find_node("start", True).get_content()

        # setting callbacks
        menu.root.find_node("general_amber", True).get_content().register_pressed_callback(ff_selected)
        menu.root.find_node("ghemical", True).get_content().register_pressed_callback(ff_selected)
        menu.root.find_node("mmff94", True).get_content().register_pressed_callback(ff_selected)
        menu.root.find_node("mmff94s", True).get_content().register_pressed_callback(ff_selected)
        universal_btn = menu.root.find_node("universal", True).get_content()
        universal_btn.register_pressed_callback(ff_selected)
        universal_btn.selected = True
        self.__selected_ff_btn = universal_btn
        menu.root.find_node("add_steps", True).get_content().register_pressed_callback(change_steps)
        menu.root.find_node("remove_steps", True).get_content().register_pressed_callback(change_steps)
        self.__steepest_descent_btn = menu.root.find_node("steepest", True).get_content()
        self.__steepest_descent_btn.register_pressed_callback(switch_steepest)
        self.__steepest_descent_btn.selected = True
        self.__conjugate_gradient_btn = menu.root.find_node("CG", True).get_content()
        self.__conjugate_gradient_btn.register_pressed_callback(switch_conjugate_gradient)
        self.__conjugate_gradient_btn.selected = False

        self.__start_btn.register_pressed_callback(switch_minimization)

    def __get_selected_forcefield(self):
        if self.__selected_ff_btn == None:
            return 'Uff'
        return self.__selected_ff_btn.text.value_idle
