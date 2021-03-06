#include <iostream>
#include <chrono>
#include <algorithm>
#include <thread>
#include <condition_variable>
#include <math.h>
#include <nanogui/nanogui.h>
#include <CGL/CGL.h>
#include "stb_image.h"

#include "smoke_screen.h"
#include "common.h"
#include "grid.h"
#include "color.h"
#include "shader.h"

using namespace nanogui;
using namespace std;

Grid grid;
GLFWwindow *window = nullptr;
Screen *screen = nullptr;
vector<Vector2D> external_forces;

static Vector3D rgb;
static bool test = true;

extern void set_callback(GLFWwindow *);

extern void error_callback(int error, const char *);

int main() {

    grid = Grid(Con::NUMCOL + 2, Con::NUMROW + 2);
    external_forces.resize(grid.width * grid.height, Vector2D(0, 0));

    glfwSetErrorCallback(error_callback);
    if (!glfwInit()) {
        return -1;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    window = glfwCreateWindow(Con::WINDOW_WIDTH, Con::WINDOW_HEIGHT, "Smoke Simulation", nullptr, nullptr);
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

    build_shader_program();
    generate_vertices_array();

    screen = new SmokeScreen(window);
    set_callback(window);

    glfwSwapInterval(1);
    glfwSwapBuffers(window);
    auto last_time = steady_clock::now();
    steady_clock::time_point rendering_start_time, rendering_end_time, simulation_start_time, simulation_end_time, cur_time;
    long long rendering_time, simulation_time;
    Vector2D prev_cursor = grid.cursor_pos;
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        glClear(GL_COLOR_BUFFER_BIT);
        if (Con::reset) {
            Con::reset = false;
            fill(external_forces.begin(), external_forces.end(), Vector2D(0, 0));
        }

        if (Con::mouse_down && prev_cursor.x == grid.cursor_pos.x && prev_cursor.y == grid.cursor_pos.y) {
            double xpos = grid.cursor_pos.x;
            double ypos = grid.cursor_pos.y;

            int row = int(Con::NUMROW - Con::NUMROW * ypos / double(Con::WINDOW_HEIGHT));
            int col = int(Con::NUMCOL * xpos / double(Con::WINDOW_WIDTH));

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

                }
            }
        } else {
            prev_cursor = grid.cursor_pos;
        }

        cur_time = steady_clock::now();
        auto elapsed = duration_cast<milliseconds>(cur_time - last_time);

        if (Con::debug) simulation_start_time = steady_clock::now();
        if (!Con::is_pause && (Con::FREQ * elapsed.count() >= 1000)) {
            last_time = cur_time;
            grid.simulate(1, external_forces, Con::ambient_temperature, Con::temperature_parameter,
                          Con::smoke_density_parameter,
                          Con::external_force_parameter, Con::num_iter);
        }

        if (Con::debug) {
            simulation_end_time = steady_clock::now();
            rendering_start_time = steady_clock::now();
        }

        unsigned char data[Con::NUMROW * Con::NUMCOL * 3] = {0};
        if (Con::is_modify_vf) {
            for (int y = 0; y < Con::NUMROW; ++y) {
                for (int x = 0; x < Con::NUMCOL; ++x) {
                    Vector2D accumulated_direction = Vector2D(0.0, 0.0);
                    for (int ys = -1 + y; ys <= y + 1; ++ys) {
                        for (int xs = -1 + x; xs <= x + 1; ++xs) {
                            if (ys >= 0 && xs >= 0 && ys < Con::NUMROW && xs < Con::NUMCOL) {
                                accumulated_direction += external_forces[ys * grid.width + xs];
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
                        double angle =
                                accumulated_direction.y >= 0 ? acos(accumulated_direction.x) :
                                acos(accumulated_direction.x) + PI;
                        angle = angle * 180 / PI;
                        hue = angle;
                    }

                    rgb = hsv2rgb({hue, saturate, value});

                    int index = y * Con::NUMCOL + x;
                    data[index * 3] = (unsigned char) (rgb.x * 255);
                    data[index * 3 + 1] = (unsigned char) (rgb.y * 255);
                    data[index * 3 + 2] = (unsigned char) (rgb.z * 255);
                }
            }
        } else {
            for (int y = 0; y < Con::NUMROW; ++y) {
                for (int x = 0; x < Con::NUMCOL; ++x) {
                    double density = grid.getDensity(x, y);
                    double temperature = grid.getTemperature(x, y);

                    double hue_center = 400;
                    double hue_halfspan = 50;
                    if (Con::picked_rgb.norm() > Con::EPS) {
                        Vector3D picked_hsv = rgb2hsv(Con::picked_rgb);
                        hue_center = picked_hsv.x;
                    }
                    double hue = (int) (hue_center - (temperature - hue_halfspan)) % 360;
                    double saturate = 100.0;
                    double value = density;

                    rgb = hsv2rgb({hue, saturate, value});

                    int index = y * Con::NUMCOL + x;
                    data[index * 3] = (unsigned char) (rgb.x * 255);
                    data[index * 3 + 1] = (unsigned char) (rgb.y * 255);
                    data[index * 3 + 2] = (unsigned char) (rgb.z * 255);
                }
            }
        }

        glBindTexture(GL_TEXTURE_2D, Con::texture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, Con::NUMCOL, Con::NUMROW, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, Con::texture);
        glUseProgram(Con::shader_program);
        glBindVertexArray(Con::VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        if (Con::debug) {
            rendering_end_time = steady_clock::now();
            rendering_time = duration_cast<milliseconds>(rendering_end_time - rendering_start_time).count();
            simulation_time = duration_cast<milliseconds>(simulation_end_time - simulation_start_time).count();
        }

        screen->drawContents();
        screen->drawWidgets();

        glfwSwapBuffers(window);

    }
    glDeleteVertexArrays(1, &Con::VAO);
    glDeleteBuffers(1, &Con::VBO);
    glDeleteBuffers(1, &Con::EBO);
    glDeleteTextures(1, &Con::texture);

    glfwTerminate();
    return 0;
}