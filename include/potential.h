#ifndef POTENTIAL_H
#define POTENTIAL_H

#include<vector>
#include<omp.h>
#include "nanoparticle.h"

class Potential {

public:

    Potential(double param, double cutoff):
        m_total_energy{0.},
        m_cut_off{cutoff},
        m_parameter{param}
        {   }


    /*!< \brief Sets the cutoff distance for the interaction.*/
    double cutOffDistance() const {return m_cut_off;}
    /*!< \brief Sets the interaction's parameter.*/
    double parameter() const {return m_parameter;}

    /*!< \brief Sets the total energy of the system for a given interaction.*/
    void setEnergy(const double& e){m_total_energy = e;}
    /*!< \brief Returns the total energy of the system for a given interaction.*/
    double totalEnergy() const {return m_total_energy;}

    /*!< \brief Computes the interaction between 2 particles. A pure virtual
            function which every interaction class implements.*/
    virtual double interactionBetween(const Nanoparticle& selected,
                                      const Nanoparticle& other) const = 0;
    /*!< \brief Computes the interaction between an indexed particle with
     all the other particles of the system that are in a vector(array).*/
    double appliedTo(size_t index, const std::vector<Nanoparticle>& vec,
                                              size_t startIndex = 0) const;
    /*!< \brief Computes the total energy of the system for one interaction.*/
    double totalInteraction(const std::vector<Nanoparticle>&) const;

private:

    double m_total_energy; /*!< brief The total energy of the system for a given interaction*/
    double m_cut_off; /*!< brief The cut off distance for a given interaction*/
    double m_parameter; /*!< brief The interaction's parameter (scale)*/
};

#endif 
