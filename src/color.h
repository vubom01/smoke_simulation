#ifndef color_h
#define color_h

#include "CGL/vector3D.h"
#include "CGL/vector4D.h"

using namespace CGL;
using std::min;
using std::max;

static inline Vector3D clamp(Vector3D col, double mi, double mx) {
    return Vector3D(min(max(mi, col.x), mx), min(max(mi, col.y), mx), min(max(mi, col.z), mx));
}

static inline Vector3D mix(Vector3D x, Vector3D y, double a) {
    return x * (1 - a) + y * a;
}

static inline Vector3D fract(Vector3D col) {
    return Vector3D(fmod(col.x, 1.0), fmod(col.y, 1.0), fmod(col.z, 1.0));
}

Vector3D hsv2rgb(Vector3D HSV) {

    Vector3D RGB;
    double H = HSV.x, S = HSV.y, V = HSV.z, P, Q, T, fract;

    (H == 360.) ? (H = 0.) : (H /= 60.);
    fract = H - floor(H);

    P = V * (1. - S);
    Q = V * (1. - S * fract);
    T = V * (1. - S * (1. - fract));

    if (0. <= H && H < 1.)
        RGB = {V, T, P};
    else if (1. <= H && H < 2.)
        RGB = {Q, V, P};
    else if (2. <= H && H < 3.)
        RGB = {P, V, T};
    else if (3. <= H && H < 4.)
        RGB = {P, Q, V};
    else if (4. <= H && H < 5.)
        RGB = {T, P, V};
    else if (5. <= H && H < 6.)
        RGB = {V, P, Q};
    else
        RGB = {0., 0., 0};

    return RGB / 256.0;
}

#endif