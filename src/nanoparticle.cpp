#include "../include/nanoparticle.h"

double zeeman(const Nanoparticle& p)
{
    return -dot(p.dipole_moment * p.dipole_length(), field);
}
double anisotropy(const Nanoparticle& p)
{
    double e = dot(p.dipole_moment, p.easy_axis);
    e /= dot(p.dipole_moment * p.dipole_length(), p.dipole_moment * p.dipole_length());
    return KANIS * p.particleVolume() * (1. - e*e);
}