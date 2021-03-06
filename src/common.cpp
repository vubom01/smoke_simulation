#include "common.h"

int Con::WINDOW_WIDTH = 1200;
int Con::WINDOW_HEIGHT = 720;

bool Con::mouse_down = false;
bool Con::is_pause = false;
bool Con::shift_pressed = false;
bool Con::is_modify_vf = false;
bool Con::reset = false;
bool Con::debug = false;
int Con::size_mouse = 3 * (Con::NUMROW / 100);

int Con::size_smoke = 6;
double Con::amount_smoke = 90;
double Con::amount_temperature = 50;
double Con::ambient_temperature = 0;
double Con::temperature_parameter = 0.010;
double Con::smoke_density_parameter = 0.005;
double Con::external_force_parameter = 0.5;
double Con::num_iter = 16;
double Con::ALPHA = 0.7;

Vector2D Con::enter_cell = Vector2D(0, 0);
Vector2D Con::exit_cell = Vector2D(0, 0);

Vector3D Con::picked_rgb = Vector3D(0, 0, 0);

std::random_device Con::rd;
mt19937 Con::rng = mt19937(Con::rd());

GLuint Con::VAO = 0;
GLuint Con::VBO = 0;
GLuint Con::EBO = 0;
GLuint Con::texture = 0;
GLuint Con::shader_program = 0;

const GLchar* Con::vertexShaderSource = R"glsl(
    #version 330 core
    layout (location = 0) in vec3 aPos;
    layout (location = 1) in vec3 aColor;
    layout (location = 2) in vec2 aTexCoord;
    out vec3 ourColor;
    out vec2 TexCoord;
    void main()
    {
        gl_Position = vec4(aPos, 1.0);
        ourColor = aColor;
        TexCoord = vec2(aTexCoord.x, aTexCoord.y);
    }
  )glsl";

const GLchar* Con::fragmentShaderSource = R"glsl(
  #version 330 core
  out vec4 FragColor;
  in vec3 ourColor;
  in vec2 TexCoord;
  // texture sampler
  uniform sampler2D texture1;
  void main()
  {
      FragColor = texture(texture1, TexCoord);
  }
)glsl";