#ifndef VANDERWAALS_H
#define VANDERWAALS_H

#include<cmath>
#include "potential.h"

class Vanderwaals : public Potential {

public:

    Vanderwaals(const double& param, const double& cutoff = RCW):
        Potential(param, cutoff){}

    virtual double interactionBetween(const Nanoparticle& selected, const Nanoparticle& other) const
    {
            double energy = 0.;

            if ( surfaceDistanceBetween(selected, other) > cutOffDistance() )
                return energy;

            double r = distanceBetween(selected, other);
            double r2 = r * r;
            double Dij = (selected.diameter() + other.diameter())/2.;
            double Dij2 = Dij * Dij;

            double e1 = Dij2/r2;
            double e2 = Dij2/(r2 - Dij2);
            double e3 = fabs((r2 - Dij2)/r2);

            energy = e1 + e2 + 2. * log(e3);

            return parameter() * energy;
    }


};

#endif 
