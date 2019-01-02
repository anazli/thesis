#ifndef METROPOLIS_H
#define METROPOLIS_H

#include<algorithm>
#include<cmath>
#include<vector>
#include<fstream>
#include<string>

#include "interactions.h"
#include "random.h"
#include "nanoparticle.h"

class Metropolis : public Random {

public:

    Metropolis(const double& d1, const double& d2,
               const unsigned int& nst, int seed, int Nsize):
        delta1{d1}, delta2{d2}, Nsteps{nst}, Random(seed, Nsize)
        {
            Naccpt = 0;
        }

        unsigned int steps() const {return Nsteps;}
        unsigned int acceptedSteps() const {return Naccpt;}
        double logNormal(const double& x) const
        {return exp(-pow(log(x/MODAL),2)/(2. * pow(SIGMA,2)))/(sqrt(2. * PI) * SIGMA * x);}
        double boltzmann(const double& x) const {return exp(-x * BETA);}
        bool theyDontOverlap(const std::vector<Nanoparticle>&);
        bool theyDontOverlap(size_t, const std::vector<Nanoparticle>&);
        void setPositions(std::vector<Nanoparticle>&);
        void agglomerates( std::ofstream&, std::ofstream&, std::vector<Nanoparticle>, double );
        double particleVolumeFraction( std::vector<Nanoparticle>& );
        double mrt2(double);
        std::vector<Nanoparticle> mrt2(size_t, const std::vector<Nanoparticle>&,
                                      const std::vector<std::string>&);

private:

    double delta1;
    double delta2;
    unsigned int Nsteps;
    unsigned int Naccpt;
};


double totalZeemann(const std::vector<Nanoparticle>&);
double totalAnisotropy(const std::vector<Nanoparticle>&);

double total_potential(size_t, const std::vector<Nanoparticle>&,
                       const std::vector<std::string>&);

double total_potential(const std::vector<Nanoparticle>&,
                       const std::vector<std::string>&);

#endif 
