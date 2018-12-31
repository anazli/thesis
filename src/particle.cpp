#include "../include/particle.h"

using std::ofstream; using std::endl;

//! The default constructor.
/*! When no arguments are provided, a particle at (0,0,0) with 0 diameter
is created.*/
Particle::Particle(const Vec3& origin, const double& d):
    center(origin),
    m_diameter(d)
{

}

//! The copy constructor
Particle::Particle(const Particle& p):
    center(p.center),
    m_diameter(p.diameter())
{

}

/*!< \brief Sets particle's diameter.*/
void Particle::setDiameter(double d){m_diameter = d;}
/*!< \brief Returns particle's diameter.*/
double Particle::diameter() const {return m_diameter;}
/*!< \brief Returns particle's (Sphere) volume.*/
double Particle::particleVolume() const 
{return (4./3.) * PI * pow(diameter()/2., 3);}

/*! Returns the center to center distance between to particles.*/
double distanceBetween(const Particle& p1, const Particle& p2)
{
    return (p1.center - p2.center).magnitude();
}

/*! Returns the center to center vector between to particles.*/
Vec3 distanceVectorBetween(const Particle& p1, const Particle& p2)
{
    return p1.center - p2.center;
}

/*! Returns the surface to surface distance between to particles.*/
double surfaceDistanceBetween(const Particle& p1, const Particle& p2)
{
    return distanceBetween(p1, p2) - (p1.diameter() + p2.diameter())/2.;
}

/*! Checks if two particles overlap. True if they do, false otherwise.*/
bool theyOverlap(const Particle& p1, const Particle& p2)
{
    return surfaceDistanceBetween(p1, p2) < 0.;
}