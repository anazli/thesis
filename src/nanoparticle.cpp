#include "../include/nanoparticle.h"

double zeeman(const Nanoparticle& p)
{
    double prod = dot(p.easy_axis, p.dipole_moment);
    double m1 = p.easy_axis.magnitude();
    double m2 = p.dipole_moment.magnitude();
    double theta = acos(prod/(m1*m2));
    double prod2 = dot(p.easy_axis, field);
    double m3 = field.magnitude();
    double phi = acos(prod2/(m1*m3));
    return -p.dipole_length() * cos(theta - phi);
}
double anisotropy(const Nanoparticle& p)
{
    double prod = dot(p.easy_axis, p.dipole_moment);
    double m1 = p.dipole_moment.magnitude();
    double m2 = p.easy_axis.magnitude();
    double theta = acos(prod/(m1*m2));
    return KANIS * p.particleVolume() * sin(theta) * sin(theta);
}