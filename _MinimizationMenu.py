import nanome

class MinimizationMenu():
    def __init__(self, plugin):
        self.__plugin = plugin
        self.__menu = None
        self.__selected_ff_btn = None
        self.__nb_steps = 2500
        self.__steepest_descent = True
        self.__steps_label = None
        self.__start_btn = None
        self.__running = False

    def __update_start_btn(self, running):
        self.__start_btn.selected = running
        self.__plugin.update_content(self.__start_btn)

    def change_running_status(self, running):
        self.__running = running
        self.__update_start_btn(running)

    def start_minimization(self):
        self.change_running_status(True)
        self.__plugin.start_minimization(self.__get_selected_forcefield(), self.__nb_steps, self.__steepest_descent)

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
                self.__nb_steps = max(self.__nb_steps, 0)
            self.__steps_label.set_all_text(str(self.__nb_steps))
            self.__plugin.update_content(self.__steps_label)

        def switch_steepest(btn):
            self.__steepest_descent = not self.__steepest_descent
            btn.selected = self.__steepest_descent
            self.__plugin.update_content(btn)

        def switch_minimization(btn):
            if not self.__running:
                self.start_minimization()
            else:
                self.stop_minimization()

        # loading menus
        menu = nanome.ui.Menu.io.from_json("_MinimizationMenu.json")
        nanome.ui.Menu.set_plugin_menu(menu)
        self.__menu = menu

        # getting elements, for future update
        self.__steps_label = menu.root.find_node("steps_label", True).get_content()[0]
        self.__start_btn = menu.root.find_node("start_btn", True).get_content()[0]

        # setting callbacks
        menu.root.find_node("general_amber", True).get_content()[0].register_pressed_callback(ff_selected)
        menu.root.find_node("ghemical", True).get_content()[0].register_pressed_callback(ff_selected)
        menu.root.find_node("mmff94", True).get_content()[0].register_pressed_callback(ff_selected)
        menu.root.find_node("mmff94s", True).get_content()[0].register_pressed_callback(ff_selected)
        menu.root.find_node("universal", True).get_content()[0].register_pressed_callback(ff_selected)
        menu.root.find_node("add_steps", True).get_content()[0].register_pressed_callback(change_steps)
        menu.root.find_node("remove_steps", True).get_content()[0].register_pressed_callback(change_steps)
        menu.root.find_node("steepest", True).get_content()[0].register_pressed_callback(switch_steepest)
        menu.root.find_node("start", True).get_content()[0].register_pressed_callback(switch_minimization)

    def __get_selected_forcefield(self):
        if self.__selected_ff_btn:
            return 0
        ff_text = self.__selected_ff_btn.text.value_idle
        if ff_text == "General Amber":
            return 1
        elif ff_text == "Ghemical":
            return 2
        elif ff_text == "MMFF94":
            return 3
        elif ff_text == "MMFF94s":
            return 4
        elif ff_text == "Universal":
            return 0
        return 0