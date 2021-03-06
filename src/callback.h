#ifndef callback_h
#define callback_h

#include "grid.h"
#include "common.h"
#include "smoke_screen.h"

extern Grid grid;
extern SmokeScreen *screen;
extern vector<Vector2D> external_forces;

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods) {
    screen->mouseButtonCallbackEvent(button, action, mods);
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            Con::mouse_down = true;
        } else if (action == GLFW_RELEASE) {
            Con::mouse_down = false;
        }
    }
}

void window_size_callback(GLFWwindow *main_window, int width, int height) {
    Con::WINDOW_WIDTH = width;
    Con::WINDOW_HEIGHT = height;
}

void keyboard_callback(GLFWwindow *window, int key, int scancode, int action, int mods) {
    screen->keyCallbackEvent(key, scancode, action, mods);
    if (action == GLFW_PRESS) {
        switch (key) {
            case GLFW_KEY_LEFT_SHIFT :
                Con::shift_pressed = true;
                break;
            case GLFW_KEY_RIGHT_SHIFT:
                Con::shift_pressed = true;
                break;
            case GLFW_KEY_P:
                Con::is_pause = !Con::is_pause;
                break;
            case GLFW_KEY_M:
                Con::is_modify_vf = !Con::is_modify_vf;
                break;
            case GLFW_KEY_S:
                Con::reset = true;
                break;
            case GLFW_KEY_R:
                screen->reset_parameters();
            default:
                break;
        }
    } else if (action == GLFW_RELEASE) {
        switch (key) {
            case GLFW_KEY_LEFT_SHIFT :
                Con::shift_pressed = false;
                break;
            case GLFW_KEY_RIGHT_SHIFT:
                Con::shift_pressed = false;
                break;
            default:
                break;
        }
    }
}

void update_mouse() {
    double xpos = grid.cursor_pos.x;
    double ypos = grid.cursor_pos.y;
    int row = int(Con::NUMROW - Con::NUMROW * ypos / double(Con::WINDOW_HEIGHT));
    int col = int(Con::NUMCOL * xpos / double(Con::WINDOW_WIDTH));
    if (Con::is_modify_vf) {
        Con::exit_cell = Vector2D(col, row);
        if (Con::exit_cell.x != Con::enter_cell.x || Con::exit_cell.y != Con::enter_cell.y) {
            Vector2D direction_mouse_drag = Con::exit_cell - Con::enter_cell;
            for (int y = row - Con::size_mouse; y < row + Con::size_mouse; y++) {
                for (int x = col - Con::size_mouse; x < col + Con::size_mouse; x++) {
                    if (y < 1 || x < 1 || y >= grid.height - 1 || x >= grid.width - 1) {
                        continue;
                    }
                    external_forces[y * grid.width + x] = direction_mouse_drag.unit();
                }
            }
            Con::enter_cell = Con::exit_cell;
        }
    } else {
        Con::exit_cell = Vector2D(col, row);
        for (int y = row - Con::size_smoke; y <= row + Con::size_smoke; ++y) {
            for (int x = col - Con::size_smoke; x <= col + Con::size_smoke; ++x) {
                double dis2 = pow(y - row, 2.0) + pow(x - col, 2.0);
                if (y < 1 || y >= grid.height - 1 || x < 1 || x >= grid.width - 1 ||
                    (dis2 > Con::size_smoke * Con::size_smoke)) {
                    continue;
                }
                dis2 /= pow((Con::NUMCOL / 100.0), 2.0);
                double fall_off = 1.0 / max(dis2, 1.0);

                double den = grid.getDensity(x, y);
                double temp = grid.getTemperature(x, y);
                grid.setDensity(x, y, min(den + Con::amount_smoke * fall_off, 100.0));
                grid.setTemperature(x, y, min(temp + Con::amount_temperature * fall_off, 100.0));

                if (Con::exit_cell.x != Con::enter_cell.x || Con::exit_cell.y != Con::enter_cell.y) {
                    Vector2D mouse_velocity = (Con::exit_cell - Con::enter_cell) * 1.0;
                    grid.setVelocity(x, y,
                                     mouse_velocity * Con::ALPHA + grid.getVelocity(x, y) * (1 - Con::ALPHA));
                }
            }
        }
        Con::enter_cell = Con::exit_cell;
    }

}

void cursor_position_callback(GLFWwindow *window, double xpos, double ypos) {
    screen->cursorPosCallbackEvent(xpos, ypos);
    glfwGetCursorPos(window, &grid.cursor_pos.x, &grid.cursor_pos.y);
    if (Con::is_modify_vf && !Con::mouse_down) {
        double xpos = grid.cursor_pos.x;
        double ypos = grid.cursor_pos.y;

        int row = int(Con::NUMROW - Con::NUMROW * ypos / double(Con::WINDOW_HEIGHT));
        int col = int(Con::NUMCOL * xpos / double(Con::WINDOW_WIDTH));

        Con::enter_cell = Vector2D(col, row);
    }
    if (Con::mouse_down) {
        update_mouse();
    }
}

#endif