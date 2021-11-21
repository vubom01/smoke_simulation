#include <memory>
#include "common.h"
#include "grid.h"

double interpolate(double d1, double d2, double s);

Vector2D interpolate(Vector2D d1, Vector2D d2, double s);

Grid::Grid(int width, int height) {
    this->width = width;
    this->height = height;
    this->density.resize(width * height, 0.0);
    this->temperature.resize(width * height, 0.0);

    normal_distribution<double> dis_v_y(2, 5);
    normal_distribution<double> dis_v_x(0, 3);
    this->velocity.resize(width * height, Vector2D(0, 0));
    for (int i = 0; i < width * height; ++i) {
        this->velocity[i] = Vector2D(dis_v_x(Con::rng), dis_v_y(Con::rng));
    }
}

Grid::Grid(const Grid &grid) {
    height = grid.height;
    width = grid.width;
    density = grid.density;
    velocity = grid.velocity;
    temperature = grid.temperature;
}

Grid &Grid::operator=(const Grid &grid) {
    height = grid.height;
    width = grid.width;
    density = grid.density;
    velocity = grid.velocity;
    temperature = grid.temperature;
    return *this;
}

Grid::Grid(Grid &&grid) {
    height = grid.height;
    width = grid.width;
    density = move(grid.density);
    velocity = move(grid.velocity);
    temperature = move(grid.temperature);
}

Grid &Grid::operator=(Grid &&grid) {
    height = grid.height;
    width = grid.width;
    density = move(grid.density);
    velocity = move(grid.velocity);
    temperature = move(grid.temperature);
    return *this;
}

void Grid::simulate(double timestep, const vector<Vector2D>& external_forces, const double ambient_temperature, const double temperature_parameter, const double smoke_density_parameter, const double external_force_parameter, const double num_iter) {

    vector<double> new_density = simulate_density(timestep);

    vector<double> new_temperature = simulate_temperature(timestep);

    vector<Vector2D> new_velocity = simulate_velocity(timestep, external_forces, ambient_temperature, temperature_parameter, smoke_density_parameter, external_force_parameter, num_iter);

    this->density = new_density;
    this->velocity = new_velocity;
    this->temperature = new_temperature;
}

vector<double> Grid::simulate_density(const double timestep) {
    vector<double> combined_density(width * height, 0.0);
    vector<double> advection_grid(width * height, 0.0);
#pragma omp parallel for
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            Vector2D reverse_velocity = -getVelocity(x, y) * timestep;
            if (x + reverse_velocity.x < 0 || x + reverse_velocity.x > width - 2 || y + reverse_velocity.y < 0 ||
                y + reverse_velocity.y > height - 2) {
                advection_grid[y * width + x] = 0.0;
            } else {
                double newx = x + reverse_velocity.x;
                double newy = y + reverse_velocity.y;
                Vector2D newloc = Vector2D(newx, newy);

                Vector2D bl = Vector2D(int(newloc.x), int(newloc.y));
                Vector2D br = bl + Vector2D(1, 0);
                Vector2D tl = bl + Vector2D(0, 1);
                Vector2D tr = bl + Vector2D(1, 1);

                double s = newx - tl.x;
                double t = newy - bl.y;
                double tlerp = interpolate(getDensity(tl), getDensity(tr), s);
                double blerp = interpolate(getDensity(bl), getDensity(br), s);
                double vlerp = interpolate(blerp, tlerp, t);

                advection_grid[y * width + x] += vlerp;
            }
        }
    }

    combined_density = advection_grid;
    return combined_density;
}

vector<double> Grid::simulate_temperature(const double timestep) {
    vector<double> advection_grid(width * height, 0.0);
#pragma omp parallel for
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            Vector2D reverse_velocity = -getVelocity(x, y) * timestep;
            if (x + reverse_velocity.x < 0 || x + reverse_velocity.x > width - 2 || y + reverse_velocity.y < 0 ||
                y + reverse_velocity.y > height - 2) {
                advection_grid[y * width + x] = 0.0;
            } else {
                double newx = x + reverse_velocity.x;
                double newy = y + reverse_velocity.y;
                Vector2D newloc = Vector2D(newx, newy);

                Vector2D bl = Vector2D(int(newloc.x), int(newloc.y));
                Vector2D br = bl + Vector2D(1, 0);
                Vector2D tl = bl + Vector2D(0, 1);
                Vector2D tr = bl + Vector2D(1, 1);

                double s = newx - tl.x;
                double t = newy - bl.y;
                double tlerp = interpolate(getTemperature(tl), getTemperature(tr), s);
                double blerp = interpolate(getTemperature(bl), getTemperature(br), s);
                double vlerp = interpolate(blerp, tlerp, t);

                advection_grid[y * width + x] += vlerp;
            }
        }
    }

    return advection_grid;
}

vector<Vector2D> Grid::simulate_velocity(double timestep, const vector<Vector2D>& external_forces, const double ambient_temperature, const double temperature_parameter, const double smoke_density_parameter, const double external_force_parameter, const double num_iter) {
    vector<Vector2D> combined_velocity(width * height, Vector2D(0, 0));
    vector<Vector2D> self_advection_grid(width * height, Vector2D(0, 0));
#pragma omp parallel for
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            Vector2D reverse_velocity = -getVelocity(x, y) * timestep;
            if (x + reverse_velocity.x < 0 || x + reverse_velocity.x > width - 2 || y + reverse_velocity.y < 0 ||
                y + reverse_velocity.y > height - 2) {
                continue;
            } else {
                double newx = x + reverse_velocity.x;
                double newy = y + reverse_velocity.y;
                Vector2D newloc = Vector2D(newx, newy);

                Vector2D bl = Vector2D(int(newloc.x), int(newloc.y));
                Vector2D br = bl + Vector2D(1, 0);
                Vector2D tl = bl + Vector2D(0, 1);
                Vector2D tr = bl + Vector2D(1, 1);

                double s = newx - tl.x;
                double t = newy - bl.y;

                Vector2D tlerp = interpolate(getVelocity(tl), getVelocity(tr), s);
                Vector2D blerp = interpolate(getVelocity(bl), getVelocity(br), s);
                Vector2D vlerp = interpolate(blerp, tlerp, t);

                self_advection_grid[y * width + x] += vlerp;

            }
        }
    }

    set_boundary_conditions(self_advection_grid, -1);

    vector<Vector2D> viscous_velocity_grid(width * height, Vector2D(0,0));
    vector<Vector2D> tem = self_advection_grid;

    double alpha = 1 / (timestep * num_iter);
    double beta = 4 + alpha;
    for (int iter = 0; iter < num_iter; ++iter) {
#pragma omp parallel for
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                // Ignore boundaries for now
                if (x == 0 || x == width - 1 || y == 0 || y == height - 1) {
                    continue;
                }
                Vector2D &l = tem[y * width + x - 1];
                Vector2D &r = tem[y * width + x + 1];
                Vector2D &u = tem[(y + 1) * width + x];
                Vector2D &b = tem[(y - 1) * width + x];
                Vector2D &center = tem[y * width + x];

                Vector2D new_velocity = (l + r + u + b + alpha * center) / beta;
                viscous_velocity_grid[y * width + x] = new_velocity;
            }
        }
        tem.swap(viscous_velocity_grid);
    }

    set_boundary_conditions(viscous_velocity_grid, -1);

    Vector2D buoyant_direction = Vector2D(0, 1);
#pragma omp parallel for
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            Vector2D buoyant_force = (-smoke_density_parameter * getDensity(x, y) + (getTemperature(x, y) - ambient_temperature)*timestep*temperature_parameter)*buoyant_direction;
            viscous_velocity_grid[y*width + x] += buoyant_force;
            viscous_velocity_grid[y*width + x] += external_force_parameter * external_forces[y*width + x];
        }
    }

    set_boundary_conditions(viscous_velocity_grid, -1);

    vector<double> divergence(width * height, 0.0);
    double halfrdx = 0.5;
#pragma omp parallel for
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            Vector2D wL = viscous_velocity_grid[x - 1 + y*width];
            Vector2D wR = viscous_velocity_grid[x + 1 + y*width];
            Vector2D wB = viscous_velocity_grid[x + (y-1)*width];
            Vector2D wT = viscous_velocity_grid[x + (y+1)*width];

            divergence[y * width + x] = halfrdx * ((wR.x - wL.x) + (wT.y - wB.y));
        }
    }

    vector<double> pressure(width * height, 0.0);
    vector<double> tem_2(width * height, 0.0);
    alpha = -1;
    beta = 4;
    for (int iter = 0; iter < num_iter; ++iter) {
#pragma omp parallel for
        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                double &l = tem_2[y * width + x - 1];
                double &r = tem_2[y * width + x + 1];
                double &u = tem_2[(y + 1) * width + x];
                double &b = tem_2[(y - 1) * width + x];
                double &center = divergence[y * width + x];

                double new_pressure = (l + r + u + b + alpha * center) / beta;
                pressure[y * width + x] = new_pressure;
            }
        }
        tem_2.swap(pressure);
    }

    set_boundary_conditions(pressure, 1);

#pragma omp parallel for
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            double &pL = pressure[y * width + x - 1];
            double &pR = pressure[y * width + x + 1];
            double &pB = pressure[(y - 1) * width + x];
            double &pT = pressure[(y + 1) * width + x];

            viscous_velocity_grid[y * width + x] -= halfrdx * Vector2D(pR - pL, pT - pB);
        }
    }

    set_boundary_conditions(viscous_velocity_grid, -1);

    return viscous_velocity_grid;
}

void Grid::set_boundary_conditions(vector<Vector2D> &vec, int b) {
    for (int i = 1; i < width - 1; i++) {
        vec[cell(i, 0)] = Vector2D(0, 1);
        vec[cell(i, height - 1)] = Vector2D(0, 1);
    }
    for (int i = 1; i < height - 1; i++) {
        vec[cell(0, i)] = Vector2D(0, 1);
        vec[cell(width - 1, i)] = Vector2D(0, 1);
    }
}

void Grid::set_boundary_conditions(vector<double> &vec, int b) {
    for (int i = 1; i < width - 1; i++) {
        vec[cell(i, 0)] = b * vec[cell(i, 1)];
        vec[cell(i, height - 1)] = b * vec[cell(i, height - 2)];
    }
    for (int i = 1; i < height - 1; i++) {
        vec[cell(0, i)] = b * vec[cell(1, i)];
        vec[cell(width - 1, i)] = b * vec[cell(width - 2, i)];
    }
}

double interpolate(double d1, double d2, double s) {
    return (1 - s) * d1 + s * d2;
}

Vector2D interpolate(Vector2D d1, Vector2D d2, double s) {
    return (1 - s) * d1 + s * d2;
}

int Grid::cell(int x, int y) {
    return y*width + x;
}

void Grid::setDensity(int x, int y, double den) {
    this->density[y * width + x] = den;
}

void Grid::setVelocity(int x, int y, Vector2D vel) {
    this->velocity[y * width + x] = vel;
}

void Grid::setTemperature(int x, int y, double temp) {
    this->temperature[y * width + x] = temp;
}

void Grid::printGrid() {
    for (int y = height - 1; y >= 0; y--) {
        string s = "";
        for (int x = 0; x < width; x++) {
            s.append(to_string(getDensity(x, y)).substr(0, 5));
            s.push_back(' ');
        }
        cout << s << endl;
    }
    cout << "____________" << endl;
}

void randomize_grid(Grid &grid, int num_speckle = 3, int size = 3) {
    uni_dis dis_x(0, Con::NUMCOL - size); // uniform distribution in C++11
    uni_dis dis_y(0, Con::NUMROW - size); // uniform distribution in C++11
    uni_dis dis_density(25, 75); // uniform distribution in C++11
    uni_dis dis_size(1, size);
    while (num_speckle--) {
        int chosen_x = dis_x(Con::rng);
        int chosen_y = dis_y(Con::rng);
        int chosen_size = dis_size(Con::rng);
        double chosen_density = dis_density(Con::rng);
        for (int i = 0; i < chosen_size; ++i) {
            for (int j = 0; j < chosen_size; ++j) {
                grid.setDensity(chosen_x + i, chosen_y + j,
                                grid.getDensity(chosen_x + i, chosen_y + j) + chosen_density);
            }
        }
    }
}