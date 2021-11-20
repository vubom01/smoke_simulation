#ifndef callback_h
#define callback_h

extern Grid grid;
extern bool mouse_down;
extern bool is_pause;
extern bool shift_pressed;
extern bool is_modify_vf;
extern bool reset;
extern Vector2D enter_cell;
extern Vector2D exit_cell;
extern int size_smoke;
extern double amount_smoke;
extern nanogui::Screen *screen;

void cursor_position_callback(GLFWwindow *window, double xpos, double ypos) {
    screen->cursorPosCallbackEvent(xpos, ypos);
    glfwGetCursorPos(window, &grid.cursor_pos.x, &grid.cursor_pos.y);
}

void error_callback(int error, const char *description) {
    puts(description);
}

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods) {
    screen->mouseButtonCallbackEvent(button, action, mods);
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            mouse_down = true;
        } else if (action == GLFW_RELEASE) {
            mouse_down = false;
        }
    }
}

void window_size_callback(GLFWwindow *main_window, int width, int height) {
    WINDOW_WIDTH = width;
    WINDOW_HEIGHT = height;
}

#endif