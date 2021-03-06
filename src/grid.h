#ifndef grid_h
#define grid_h

#include <vector>

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector2D.h"

using namespace CGL;
using namespace std;

class Grid {

public:

    Grid() {}

    Grid(int width, int height);

    Grid(const Grid& grid);
    Grid& operator=(const Grid& grid);
    Grid(Grid&& grid);
    Grid& operator=(Grid&& grid);

    ~Grid() {}

    int height;
    int width;
    Vector2D cursor_pos;

    void simulate(double timestep, const vector<Vector2D>& external_forces, const double ambient_temperature, const double temperature_parameter, const double smoke_density_parameter, const double external_force_parameter, const double num_iter);

private:
    vector<Vector2D> simulate_velocity(double timestep, const vector<Vector2D>& external_forces, const double ambient_temperature, const double temperature_parameter, const double smoke_density_parameter, const double external_force_parameter, const double num_iter);
    vector<double> simulate_density(double timestep);
    vector<double> simulate_temperature(double timestep);

    void set_boundary_conditions(vector<Vector2D> &vec, int b);
    void set_boundary_conditions(vector<double> &vec, int b);

    int cell(int x, int y);

    vector<double> density;
    vector<double> temperature;
    vector<Vector2D> velocity;


public:

    double getDensity(int x, int y) const { return density[y * width + x]; }

    double getDensity(Vector2D vec) const {
        return density[vec.y * width + vec.x];
    };

    Vector2D getVelocity(int x, int y) const { return velocity[y * width + x]; }

    Vector2D getVelocity(Vector2D vec) const {
        return velocity[vec.y * width + vec.x];
    };

    double getTemperature(int x, int y) const { return temperature[y * width + x]; }

    double getTemperature(Vector2D vec) const { return temperature[vec.y * width + vec.x]; }

    void setDensity(int x, int y, double den);

    void setVelocity(int x, int y, Vector2D velocity);

    void setTemperature(int x, int y, double temp);

    void printGrid();
};


#endif