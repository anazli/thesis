#ifndef DIPOLAR_H
#define DIPOLAR_H

#include "potential.h"

class Dipolar : public Potential {

public:

    Dipolar(const double& param, const double& cutoff = RCM):
        Potential(param, cutoff){}

    virtual double interactionBetween(const Nanoparticle& selected, const Nanoparticle& other) const
    {
        double energy = 0.;

        if(surfaceDistanceBetween(selected, other) > cutOffDistance())
            return energy;

        Vec3 v = distanceVectorBetween(selected, other);
        double r = distanceBetween(selected, other);
        double r3 = r * r * r;
        double r5 = r3 * r * r;

        double prod1 = dot(selected.dipole_moment * selected.dipole_length(),
                           other.dipole_moment * other.dipole_length());
        double prod2 = dot(selected.dipole_moment * selected.dipole_length(), v);
        double prod3 = dot(other.dipole_moment * other.dipole_length(), v);

        double e1 = prod1/r3;
        double e2 = (prod2 * prod3)/r5;

        energy = e1 - 3. * e2;

        return parameter() * energy;

    }

};

#endif // DIPOLAR_H
