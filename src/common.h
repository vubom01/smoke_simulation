#ifndef common_h
#define common_h

#include "chrono"
#include "grid.h"

#include <random>
#include <nanogui/nanogui.h>
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include "CGL/vector4D.h"

using std::mt19937;
using std::uniform_int_distribution;
using uni_dis = uniform_int_distribution<int>;
using namespace std::chrono;

class Con {

public:

    static constexpr int NUMCOL = 420;
    static constexpr int NUMROW = 250;

    static constexpr int FREQ = 30;
    static constexpr double EPS = 1e-3;

    static int WINDOW_WIDTH;
    static int WINDOW_HEIGHT;

    static bool mouse_down;
    static bool is_pause;
    static bool shift_pressed;
    static bool is_modify_vf;
    static bool reset;
    static bool debug;
    static int size_mouse;
    static double ALPHA;

    static int size_smoke;
    static double amount_smoke;
    static double amount_temperature;
    static double ambient_temperature;
    static double temperature_parameter;
    static double smoke_density_parameter;
    static double external_force_parameter;
    static double num_iter;

    static Vector3D picked_rgb;

    static const GLchar *vertexShaderSource;
    static const GLchar *fragmentShaderSource;

    static Vector2D enter_cell;
    static Vector2D exit_cell;

    static std::random_device rd;

    static mt19937 rng;

    static GLuint VAO, VBO, EBO;
    static GLuint texture;
    static GLuint shader_program;
};


#endif