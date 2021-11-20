#ifndef color_h
#define color_h

#include "CGL/vector3D.h"
#include "CGL/vector4D.h"

using namespace CGL;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

static inline Vector3D clamp(Vector3D col, double mi, double mx) {
    return Vector3D(std::min(std::max(mi, col.x), mx),
                    std::min(std::max(mi, col.y), mx),
                    std::min(std::max(mi, col.z), mx));
}

static inline Vector3D mix(Vector3D x, Vector3D y, double a) {
    return x * (1 - a) + y * a;
}

static inline Vector3D fract(Vector3D col) {
    return Vector3D(fmod(col.x, 1.0), fmod(col.y, 1.0), fmod(col.z, 1.0));
}

Vector3D hsv2rgb(Vector3D in)
{
    double      hh, p, q, t, ff;
    long        i;
    Vector3D         out;

    in.y /= 100.0, in.z /= 100.0;

    if(in.y <= 0.0) {
        out.x = in.z;
        out.y = in.z;
        out.z = in.z;
        return out;
    }
    hh = in.x;
    if(hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = in.z * (1.0 - in.y);
    q = in.z * (1.0 - (in.y * ff));
    t = in.z * (1.0 - (in.y * (1.0 - ff)));

    switch(i) {
        case 0:
            out.x = in.z;
            out.y = t;
            out.z = p;
            break;
        case 1:
            out.x = q;
            out.y = in.z;
            out.z = p;
            break;
        case 2:
            out.x = p;
            out.y = in.z;
            out.z = t;
            break;

        case 3:
            out.x = p;
            out.y = q;
            out.z = in.z;
            break;
        case 4:
            out.x = t;
            out.y = p;
            out.z = in.z;
            break;
        case 5:
        default:
            out.x = in.z;
            out.y = p;
            out.z = q;
            break;
    }
    return out;
}

Vector3D hsv2rgbb(Vector3D HSV) {

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