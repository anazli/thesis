#ifndef PARTICLE_H
#define PARTICLE_H

#include<cmath>

#include "parameters.h"
#include "vec3.h"

class Particle {

public:

    /****************************************************************
     *  INTERFACE
     * *************************************************************/

    Particle(const Vec3& origin = Vec3(), const double& d = 0.);
    Particle(const Particle& p); 
    ~Particle(){}

    void setDiameter(double d);
    double diameter()const;
    double particleVolume()const;


    /****************************************************************
     *  IMPLEMENTATION
     * *************************************************************/

    Vec3 center;/*!< \brief Particle's center, public Vec3.*/

private:
    double m_diameter;/*!< \brief Particle's diameter.*/

};

double distanceBetween(const Particle&, const Particle&);
Vec3 distanceVectorBetween(const Particle&, const Particle&);
double surfaceDistanceBetween(const Particle&, const Particle&);
bool theyOverlap(const Particle& p1, const Particle& p2);

#endif