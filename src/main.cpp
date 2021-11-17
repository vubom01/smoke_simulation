#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include <math.h>
//#include <GLFW/glfw3.h>
#include <nanogui/nanogui.h>

//#include <OpenGL/gl.h>
//#include <OpenGL/OpenGL.h>
#include <CGL/CGL.h>

#include "grid.h"
#include "common.h"
#include "callback.h"
#include "color.h"

using namespace nanogui;
using namespace std;

static std::random_device rd;

static mt19937 rng(rd());

Grid grid;
bool mouse_down = false;
bool is_pause = false;
bool shift_pressed = false;
bool is_modify_vf = false;
bool reset = false;
Vector2D enter_cell = Vector2D(0,0);
Vector2D exit_cell = Vector2D(0,0);

int size_smoke = 5;
double amount_smoke = 90;
double amount_temperature = 50;
double ambient_temperature = 0;

int size_mouse = 1;

bool test = true;

GLuint VBOs[NUMCOL * NUMROW], VAOs[NUMCOL * NUMROW], EBOs[NUMCOL * NUMROW];

GLFWwindow *window = nullptr;
Screen *screen = nullptr;

void randomize_grid(Grid &grid, int num_speckle = 3, int size = 3) {
    uni_dis dis_x(0, NUMCOL - size);
    uni_dis dis_y(0, NUMROW - size);
    uni_dis dis_density(25, 75);
    uni_dis dis_size(1, size);
    while (num_speckle--) {
        int chosen_x = dis_x(rng);
        int chosen_y = dis_y(rng);
        int chosen_size = dis_size(rng);
        double chosen_density = dis_density(rng);
        for (int i = 0; i < chosen_size; ++i) {
            for (int j = 0; j < chosen_size; ++j) {
                grid.setDensity(chosen_x + i, chosen_y + j,
                                grid.getDensity(chosen_x + i, chosen_y + j) + chosen_density);
            }
        }
    }
}

void generate_vertices_array() {
    double width = 1 / (double) NUMCOL * 2;
    double height = 1 / (double) NUMROW * 2;
    GLuint elements[] = {
            0, 1, 2,
            2, 3, 0
    };
    for (int y = 0; y < NUMROW; ++y) {
        for (int x = 0; x < NUMCOL; ++x) {
            float bottom_left_x = -1 + width * x;
            float bottom_left_y = -1 + height * y;

            float rectangle_vertices[] = {
                    bottom_left_x, bottom_left_y + float(height), 0, // top left
                    bottom_left_x + float(width), bottom_left_y + float(height), 0, // top right
                    bottom_left_x + float(width), bottom_left_y, 0, // bottom right
                    bottom_left_x + float(width), bottom_left_y, 0, // bottom right
                    bottom_left_x, bottom_left_y, 0, // bottom left
                    bottom_left_x, bottom_left_y + float(height), 0, // top left
            };
            int index = y * NUMCOL + x;
            glBindVertexArray(VAOs[index]);
            glBindBuffer(GL_ARRAY_BUFFER, VBOs[index]);
            glBufferData(GL_ARRAY_BUFFER, sizeof(rectangle_vertices), rectangle_vertices,
                         GL_STATIC_DRAW); // TODO not sure if GL_DYNAMIC_DRAW is better
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *) 0);
            glEnableVertexAttribArray(0);
        }
    }
}

void set_callback() {
    glfwSetKeyCallback(window, [](GLFWwindow *, int key, int scancode, int action, int mods) {
        screen->keyCallbackEvent(key, scancode, action, mods);
    });

    glfwSetCharCallback(window, [](GLFWwindow *, unsigned int codepoint) {
        screen->charCallbackEvent(codepoint);
    });

    glfwSetDropCallback(window, [](GLFWwindow *, int count, const char **filenames) {
        screen->dropCallbackEvent(count, filenames);
    });

    glfwSetScrollCallback(window, [](GLFWwindow *, double x, double y) {
        screen->scrollCallbackEvent(x, y);
    });

    glfwSetFramebufferSizeCallback(window, [](GLFWwindow *, int width, int height) {
        screen->resizeCallbackEvent(width, height);
    });
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetKeyCallback(window, keyboard_callback);
}

GLuint build_shader_program() {
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertex_shader);
    int success;
    char info_log[512];
    glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertex_shader, 512, NULL, info_log);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << info_log << std::endl;
    }
    GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragment_shader);
    glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragment_shader, 512, NULL, info_log);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << info_log << std::endl;
    }
    GLuint shader_program = glCreateProgram();
    glAttachShader(shader_program, vertex_shader);
    glAttachShader(shader_program, fragment_shader);
    glLinkProgram(shader_program);
    glGetProgramiv(shader_program, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shader_program, 512, NULL, info_log);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << info_log << std::endl;
    }
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);

    return shader_program;
}

int main() {
    grid = Grid(NUMCOL + 2, NUMROW + 2);

    vector<Vector2D> external_forces(grid.width*grid.height, Vector2D(0, 0));

    glfwSetErrorCallback(error_callback);
    if (!glfwInit()) {
        return -1;
    }
//    glfwSetTime(0);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#endif

    window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Smoke Simulation", nullptr, nullptr);
    if (!window) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    GLuint shader_program = build_shader_program();
    glGenVertexArrays(NUMCOL * NUMROW, VAOs);
    glGenBuffers(NUMCOL * NUMROW, VBOs);
    glGenBuffers(NUMCOL * NUMROW, EBOs);
    generate_vertices_array();
    screen = new Screen();
    screen->initialize(window, true);
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);

    bool enabled = true;
    FormHelper *gui = new FormHelper(screen);
    nanogui::ref<Window> nanoguiWindow = gui->addWindow(Eigen::Vector2i(10, 10), "Adjustable parameters");
    gui->addGroup("Smoke");
    gui->addVariable("Size of smoke (1 to 20)", size_smoke);
    gui->addVariable("Density of smoke (0 to 100)", amount_smoke);
    gui->addVariable("Temperature of smoke (0 to 100)", amount_temperature);
    gui->addVariable("Ambient temperature (0 to 100)", ambient_temperature);

    screen->setVisible(true);
    screen->performLayout();
    nanoguiWindow->center();

    set_callback();


#if defined(_OPENMP)
#pragma omp parallel
    {
        int rank, rankn;
        rank = omp_get_thread_num();
        rankn = omp_get_num_threads();
    }
#endif

    glfwSwapInterval(1);
    glfwSwapBuffers(window);
    auto last_time = steady_clock::now();

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        if (reset) {
            reset = false;
            fill(external_forces.begin(), external_forces.end(), Vector2D(0, 0));
        }

        if (is_modify_vf && !mouse_down) {
            double xpos = grid.cursor_pos.x;
            double ypos = grid.cursor_pos.y;

            int row = int(NUMROW - NUMROW * ypos / double(WINDOW_HEIGHT));
            int col = int(NUMCOL * xpos / double(WINDOW_WIDTH));

            enter_cell = Vector2D(col, row);
        }

        if (mouse_down) {
            if (is_modify_vf) {
                double xpos = grid.cursor_pos.x;
                double ypos = grid.cursor_pos.y;

                int row = int(NUMROW - NUMROW * ypos / double(WINDOW_HEIGHT));
                int col = int(NUMCOL * xpos / double(WINDOW_WIDTH));

                exit_cell = Vector2D(col, row);
                if (exit_cell.x != enter_cell.x || exit_cell.y != enter_cell.y) {
                    Vector2D direction_mouse_drag = exit_cell - enter_cell;
                    for (int y = row - size_mouse; y < row + size_mouse; y++) {
                        for (int x = col - size_mouse; x < col + size_mouse; x++) {
                            if (y < 1 || x < 1 || y >= grid.height - 1 || x >= grid.width - 1) {
                                continue;
                            }
                            external_forces[y*grid.width + x] = direction_mouse_drag.unit();
                        }
                    }
                    enter_cell = exit_cell;
                }
            } else {
                double xpos = grid.cursor_pos.x;
                double ypos = grid.cursor_pos.y;

                int row = int(NUMROW - NUMROW * ypos / double(WINDOW_HEIGHT));
                int col = int(NUMCOL * xpos / double(WINDOW_WIDTH));

                for (int y = row - size_smoke; y <= row + size_smoke; ++y) {
                    for (int x = col - size_smoke; x <= col + size_smoke; ++x) {
                        double dis2 = pow(y - row, 2.0) + pow(x - col, 2.0);

                        if (y < 1 || y >= grid.height - 1 || x < 1 || x >= grid.width - 1 ||
                            (dis2 > size_smoke * size_smoke)) {
                            continue;
                        }

                        double fall_off = 2 * 1.0 / max(dis2, 1.0);

                        double den = grid.getDensity(x, y);
                        double temp = grid.getTemperature(x, y);
                        grid.setDensity(x, y, min(den + amount_smoke * fall_off, 100.0));
                        grid.setTemperature(x, y, min(temp + amount_temperature * fall_off, 100.0));

                    }
                }
            }
        }
        auto cur_time = steady_clock::now();
        auto elapsed = duration_cast<milliseconds>(cur_time - last_time);

        if (!is_pause) {
            if (FREQ * elapsed.count() >= 1000) {
                last_time = cur_time;
                auto start_time = steady_clock::now();
                grid.simulate(1, external_forces, ambient_temperature);
                auto end_time = steady_clock::now();
                auto simulate_time = duration_cast<milliseconds>(end_time - start_time);
            } else {
            }
        }

        auto start_time = steady_clock::now();
        glUseProgram(shader_program);

        if (is_modify_vf) {
            for (int y = 0; y < NUMROW; ++y) {
                for (int x = 0; x < NUMCOL; ++x) {
                    Vector2D accumulated_direction = Vector2D(0.0,0.0);
                    for (int ys = y; ys < y + 1; ys++) {
                        for (int xs = x; xs < x + 1; xs++) {
                            if (ys >= 0 && xs >= 0 && ys < NUMROW && xs < NUMCOL) {
                                accumulated_direction += external_forces[ys*grid.width + xs];
                            }
                        }
                    }
                    double hue = 0;
                    double saturate = 100;
                    double value = 100;
                    if (accumulated_direction.x == 0 && accumulated_direction.y == 0) {
                        value = 0;
                    } else {
                        accumulated_direction = accumulated_direction.unit();
                        double angle = accumulated_direction.y >= 0 ? acos(accumulated_direction.x) : acos(accumulated_direction.x) + PI;
                        angle = angle * 180 / PI;
                        hue = angle;
                    }

                    Vector3D rgb = hsv2rgb({hue, saturate, value});

                    int index = y * NUMCOL + x;

                    int vertexColorLocation = glGetUniformLocation(shader_program, "ourColor");
                    glUniform4f(vertexColorLocation, rgb.x, rgb.y, rgb.z, 1.0f);
                    glBindVertexArray(VAOs[index]);
                    glDrawArrays(GL_TRIANGLES, 0, 6);
                }
            }
        } else {
            for (int y = 0; y < NUMROW; ++y) {
                for (int x = 0; x < NUMCOL; ++x) {
                    double density = grid.getDensity(x, y);
                    double temperature = grid.getTemperature(x, y);

                    double hue = ((int)(450.0 - temperature * 1)) % 360;
                    double saturate = 100.0;
                    double value = density;

                    Vector3D rgb = hsv2rgb({hue, saturate, value});

                    int index = y * NUMCOL + x;

                    int vertexColorLocation = glGetUniformLocation(shader_program, "ourColor");
                    glUniform4f(vertexColorLocation, rgb.x, rgb.y, rgb.z, 1.0f);
                    glBindVertexArray(VAOs[index]);
                    glDrawArrays(GL_TRIANGLES, 0, 6);
                }
            }
        }
        auto end_time = steady_clock::now();
        auto display_time = duration_cast<milliseconds>(end_time - start_time);
        screen->drawContents();
        screen->drawWidgets();

        glfwSwapBuffers(window);

    }
    glDeleteVertexArrays(NUMCOL * NUMROW, VAOs);
    glDeleteBuffers(NUMCOL * NUMROW, VBOs);
    glfwTerminate();
    return 0;
}