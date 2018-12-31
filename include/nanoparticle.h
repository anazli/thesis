#ifndef NANOPARTICLE_H
#define NANOPARTICLE_H

#include "particle.h"

class Nanoparticle : public Particle {


public:

    /****************************************************************
     *  INTERFACE
     * *************************************************************/

    
    Nanoparticle(){};//! The default constructor.
    Nanoparticle(Vec3 origin, double d, Vec3 dip, Vec3 easy, double slw):
                Particle(origin, d), dipole_moment(dip), easy_axis(easy),
                m_slw(slw), m_FoF(false)
                {dipole_magnitude = particleVolume() * MSAT;}

    /*!< \brief Sets the magnitude of dipole moment.*/
    void setDipoleMagnitude(double p){dipole_magnitude = p;}
    /*!< \brief Sets the surfactant layer width.*/
    void setSLW(double n){m_slw = n;}
    void setFoF(bool b){m_FoF = b;}/*!< \brief Sets FoF value.*/

    /*!< \brief Returns the magnitude of dipole moment.*/
    double dipole_length() const {return dipole_magnitude;}
    /*!< \brief Returns the surfactant layer width.*/
    double slw() const {return m_slw;}
    bool FoF() const {return m_FoF;}/*!< \brief Returns the FoF value.*/

    /****************************************************************
     *  IMPLEMENTATION
     * *************************************************************/

    Vec3 dipole_moment;/*!< \brief Nanoparticles' dipole moment, public Vec3.*/
    Vec3 easy_axis;/*!< \brief Nanoparticles' easy axis, public Vec3.*/

private:
    double dipole_magnitude;/*!< \brief Magnitude of the dipole_moment vector.*/
    double m_slw;/*!< \brief Surfactant layer width.*/
    bool m_FoF;/*!< \brief A value that is used in agglomerate function 
                of the metropolis class.*/
};

#endif 