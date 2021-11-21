#include "common.h"
#include "smoke_screen.h"

void SmokeScreen::reset_parameters() {
    Con::size_smoke = 6;
    Con::ambient_temperature = 0;
    Con::temperature_parameter = 0.010;
    Con::smoke_density_parameter = 0.005;
    Con::num_iter = 16;
    Con::ALPHA = 0.7;

    s_size_smoke->setValue(Con::size_smoke);
    s_ambient_temperature->setValue(Con::ambient_temperature);
    s_temperature_parameter->setValue(Con::temperature_parameter);
    s_smoke_density_parameter->setValue(Con::smoke_density_parameter);
    s_alpha->setValue(Con::ALPHA);
}