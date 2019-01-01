#ifndef STERIC_H
#define STERIC_H

#include<cmath>
#include "potential.h"

class Steric : public Potential {

public:

    Steric(const double& param, const double& cutoff = 1.):
        Potential(param, cutoff){}

    virtual double interactionBetween(const Nanoparticle& selected, const Nanoparticle& other) const
    {
        double energy = 0.;

        double s = surfaceDistanceBetween(selected, other);
        double ts = s/(selected.slw() + other.slw());

        if ( ts <= cutOffDistance() )
        {
            double Dij = (selected.diameter() + other.diameter())/2.;
            double l = 2. * s/Dij;
            double t = (selected.slw() + other.slw())/Dij;

            double e1 = 2. - l/t;
            double e2 = (l + 2.)/t;
            double e3 = (1. + t)/(1. + l/2.);

            energy = Dij * Dij * (e1 - e2 * log(e3));
        }

        return parameter() * energy;
    }

};

#endif
