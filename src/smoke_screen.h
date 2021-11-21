#include <CGL/CGL.h>
#include "common.h"

class SmokeScreen : public nanogui::Screen {
public:

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
                                                     "Adjustable parameters");

        gui->addGroup("Smoke");

        nanogui::ref<Slider> s_size_smoke = new nanogui::Slider(window);
        nanogui::ref<Slider> s_amount_smoke = new nanogui::Slider(window);
        nanogui::ref<Slider> s_amount_temperature = new nanogui::Slider(window);
        nanogui::ref<Slider> s_ambient_temperature = new nanogui::Slider(window);
        nanogui::ref<Slider> s_temperature_parameter = new nanogui::Slider(window);
        nanogui::ref<Slider> s_smoke_density_parameter = new nanogui::Slider(window);
        nanogui::ref<Slider> s_external_force_parameter = new nanogui::Slider(window);
        nanogui::ref<Slider> s_num_iter = new nanogui::Slider(window);
        nanogui::ref<Slider> s_alpha = new nanogui::Slider(window);

        s_size_smoke->setCallback([&](double ret) { Con::size_smoke = int(ret); });
        s_amount_smoke->setCallback([&](double ret) { Con::amount_smoke = ret; });
        s_amount_temperature->setCallback([&](double ret) { Con::amount_temperature = ret; });
        s_ambient_temperature->setCallback([&](double ret) { Con::ambient_temperature = ret; });
        s_temperature_parameter->setCallback([&](double ret) { Con::temperature_parameter = ret; });
        s_smoke_density_parameter->setCallback([&](double ret) { Con::smoke_density_parameter = ret; });
        s_external_force_parameter->setCallback([&](double ret) { Con::external_force_parameter = ret; });
        s_num_iter->setCallback([&](double ret) { Con::num_iter = ret; });
        s_alpha->setCallback([&](double ret) { Con::ALPHA = ret; });

        s_size_smoke->setValue(Con::size_smoke);
        s_amount_smoke->setValue(Con::amount_smoke);
        s_amount_temperature->setValue(Con::amount_temperature);
        s_ambient_temperature->setValue(Con::ambient_temperature);
        s_temperature_parameter->setValue(Con::temperature_parameter);
        s_smoke_density_parameter->setValue(Con::smoke_density_parameter);
        s_external_force_parameter->setValue(Con::external_force_parameter);
        s_alpha->setValue(Con::ALPHA);

        s_size_smoke->setRange({3.0, 30.0});
        s_amount_smoke->setRange({50.0, 100.0});
        s_amount_temperature->setRange({0, 100.0});
        s_ambient_temperature->setRange({0, 50.0});
        s_temperature_parameter->setRange({0, 0.02});
        s_smoke_density_parameter->setRange({0, 0.01});
        s_external_force_parameter->setRange({0, 1.0});
        s_num_iter->setRange({4, 32});

        gui->addWidget("Size(0 to 30)", s_size_smoke);
        gui->addWidget("Density(50 to 100)", s_amount_smoke);
        gui->addWidget("Heat(0 to 100)", s_amount_temperature);
        gui->addWidget("Ambient", s_ambient_temperature);
        gui->addWidget("Heat parameter", s_temperature_parameter);
        gui->addWidget("Density parameter", s_smoke_density_parameter);
        gui->addWidget("Force parameter", s_external_force_parameter);
        gui->addWidget("Diffusion iterations", s_num_iter);
        gui->addWidget("Velocity ratio", s_alpha);

        nanogui::ref<TabWidget> tabWidget = this->add<TabWidget>();
        nanogui::ref<Widget> layer = tabWidget->createTab("Color Wheel");
        layer->setLayout(new GroupLayout(margin));

        layer->add<Label>("Color", "sans-bold");
        nanogui::ref<ColorWheel> color_wheel = layer->add<ColorWheel>();
        color_wheel->setCallback([&](const nanogui::Color &color) {
            Con::picked_rgb = Vector3D(color.r(), color.g(), color.b());
        });

        this->setVisible(true);
        this->performLayout();

        window->setPosition(window->position() - Vector2i(window->width(), 0));
    }

private:
    double margin = 10;
};