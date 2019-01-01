#ifndef GUARD_RANDOM_H
#define GUARD_RANDOM_H

#include<random>

/*! \brief The class Random is used for the generation of random numbers used in the simulation.
 *
 *  This class uses the Mersenne Twister generator for the generation of random numbers, it
 *  consists of 3 data members, the mt19937 generator, one variable (rnd) for uniform real
 *  random numbers and one (randint) for uniform integer random numbers.
 *  It is inherited by the \ref Metropolis class which uses extensively random numbers.
 */
class Random {

public:

    Random(int seed, int Nsize):generator(seed), rnd(0., 1.), randint(0, Nsize) { }

    ~Random(){}

    double random()
    {return rnd(generator);}

    int randomInt()
    {return randint(generator);}


private:

    std::mt19937 generator; 		    /*!< \brief The mt19937 generator.*/
    std::uniform_real_distribution<double> rnd; /*!< \brief real random numbers.*/
    std::uniform_int_distribution<int> randint; /*!< \brief real integer numbers.*/

};


#endif	// Random

