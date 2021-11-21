#include <CGL/CGL.h>
#include "common.h"

class SmokeScreen : public nanogui::Screen {

public:
    void reset_parameters();

    SmokeScreen(GLFWwindow *glfw_window) {

        using namespace nanogui;
        using namespace std;

        this->initialize(glfw_window, true);
        int width, height;
        glfwGetFramebufferSize(glfw_window, &width, &height);
        glViewport(0, 0, width, height);

        bool enabled = true;
        FormHelper *gui = new FormHelper(this);
        nanogui::ref<Window> window = gui->addWindow(Eigen::Vector2i(Con::WINDOW_WIDTH - margin, margin),
                                                     "Parameters");

        s_size_smoke = new nanogui::Slider(window);
        s_ambient_temperature = new nanogui::Slider(window);
        s_temperature_parameter = new nanogui::Slider(window);
        s_smoke_density_parameter = new nanogui::Slider(window);
        s_num_iter = new nanogui::Slider(window);
        s_alpha = new nanogui::Slider(window);

        s_size_smoke->setCallback([&](double ret) { Con::size_smoke = int(ret); });
        s_ambient_temperature->setCallback([&](double ret) { Con::ambient_temperature = ret; });
        s_temperature_parameter->setCallback([&](double ret) { Con::temperature_parameter = ret; });
        s_smoke_density_parameter->setCallback([&](double ret) { Con::smoke_density_parameter = ret; });
        s_num_iter->setCallback([&](double ret) { Con::num_iter = ret; });
        s_alpha->setCallback([&](double ret) { Con::ALPHA = ret; });

        s_size_smoke->setValue(Con::size_smoke);
        s_ambient_temperature->setValue(Con::ambient_temperature);
        s_temperature_parameter->setValue(Con::temperature_parameter);
        s_smoke_density_parameter->setValue(Con::smoke_density_parameter);
        s_num_iter->setValue(Con::num_iter);
        s_alpha->setValue(Con::ALPHA);

        s_size_smoke->setRange({3.0, 30.0});
        s_ambient_temperature->setRange({0, 50.0});
        s_temperature_parameter->setRange({0, 0.02});
        s_smoke_density_parameter->setRange({0, 0.01});
        s_num_iter->setRange({4, 32});

        gui->addWidget("Size", s_size_smoke);
        gui->addWidget("Ambient", s_ambient_temperature);
        gui->addWidget("Heat", s_temperature_parameter);
        gui->addWidget("Density", s_smoke_density_parameter);
        gui->addWidget("Diffusion", s_num_iter);
        gui->addWidget("Velocity", s_alpha);

        nanogui::ref<Widget> layer = this->add<Widget>();
        layer->setLayout(new GroupLayout(margin));

        nanogui::ref<ColorWheel> color_wheel = layer->add<ColorWheel>();
        color_wheel->setCallback([&](const nanogui::Color &color) {
            Con::picked_rgb = Vector3D(color.r(), color.g(), color.b());
        });

        this->setVisible(true);
        this->performLayout();

        window->setPosition(window->position() - Vector2i(window->width() - 10, 10));
    }

private:

    double margin = 10;
    nanogui::ref<nanogui::Slider> s_size_smoke;
    nanogui::ref<nanogui::Slider> s_ambient_temperature;
    nanogui::ref<nanogui::Slider> s_temperature_parameter;
    nanogui::ref<nanogui::Slider> s_smoke_density_parameter;
    nanogui::ref<nanogui::Slider> s_num_iter;
    nanogui::ref<nanogui::Slider> s_alpha;

};