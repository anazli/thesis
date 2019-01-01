#include "../include/potential.h"

// Interaction between the indexed particle and all the others in the system
double Potential::appliedTo(size_t index, const std::vector<Nanoparticle>& v, size_t startIndex) const
{
	double energy = 0.;
	size_t n = v.size();
	size_t i;

#pragma omp parallel shared(v) private(i) \
	firstprivate(n, index) reduction(+:energy)

{
#pragma omp for schedule(static, chunkSZ)

	//for all particles before the indexed one
	for ( i = startIndex ; i < index ; ++i )
		energy += interactionBetween(v[index], v[i]);

#pragma omp for schedule(static, chunkSZ)

	//for all particles after the indexed one
	for ( i = (index + 1) ; i < n ; ++i )
		energy += interactionBetween(v[index], v[i]);

} //end of parallel region

	return energy;
}




// Interaction applied to the system
double Potential::totalInteraction(const std::vector<Nanoparticle>& v) const
{
	double energy = 0.;
	size_t n = v.size();

	//they have already been parallelized in function calInteraction
	for (size_t i = 0 ; i < n ; ++i) {

		energy += appliedTo(i, v, i + 1);
	}

	return energy;
}
