#include "../include/metropolis.h"
using std::vector;  using std::ofstream;
using std::endl;    using std::max;
using std::string;

bool Metropolis::theyDontOverlap(const vector<Nanoparticle>& v)
{
    size_t n = v.size();

    for ( size_t i = 0 ; i < (n-1) ; ++i )
    {
        for ( size_t j = (i+1) ; j < n ; ++j )
        {
            if ( theyOverlap(v[i], v[j]) )
            {
                return false;
            }
        }
    }

    return true;
}





bool Metropolis::theyDontOverlap(size_t index, const vector<Nanoparticle>& v)
{
    size_t n = v.size();

    for ( size_t i = 0 ; i < index ; ++i )
    {
        if ( theyOverlap(v[index], v[i]) )
            return false;
    }

    for ( size_t i = (index+1) ; i < n ; ++i )
    {
        if ( theyOverlap(v[index], v[i]) )
            return false;
    }

    return true;

}

void Metropolis::setPositions(vector<Nanoparticle>& v)
{
    size_t n = v.size();

	while( !(theyDontOverlap(v)) )
	{ // Assure that they don't overlap (All of them)
		for(size_t i = 0 ; i < n ; ++i)
        {
			while(!theyDontOverlap(i, v))
			{ //while they are still overlaping try to
			//place them in random positions until there are no overlaps.
			//but one particle at a time because it's faster.
                Vec3 r = randomVector();
                r = -L/2. + L * r;
				v[i].center = r;
			}

		}
	}
}

void Metropolis::agglomerates(ofstream& out, ofstream& aggrsize, vector<Nanoparticle> vec, double Tvol)
{

	size_t pnumber = 0, k = 0, N = vec.size();


	//initialize all indexes with 0 value
	for (size_t i = 0 ; i != N ; ++i)
		vec[i].setFoF(0);

	vector< vector<Nanoparticle> > group(N);

	while (pnumber != N) {


		//for each particle, n, not yet associated with any group
		//assign n to group k, initialize a new member list group[k]
		//for group k with particle n as a first entry
		for (size_t i = 0 ; i != N ; ++i) {

			if(vec[i].FoF() == 0) {

				group[k].push_back(vec[i]);
				vec[i].setFoF(1);
				break; //add the first with 0 index and stop
			}
		}


		//Recursively, for each new particle p in group[k]
		//find neighbores of p lying within distance less than or equal
		//to linking length, add to group[k] those not already assigned to group[k]
		for(size_t i = 0 ; i != group[k].size() ; ++i) {

			for(size_t j = 0 ; j != N ; ++j) {

				if( (vec[j].FoF() == 0) &&
				   surfaceDistanceBetween(group[k][i], vec[j]) <=
				   (group[k][i].slw() + vec[j].slw()) ) {

				group[k].push_back(vec[j]);
				vec[j].setFoF(1);

	 			}

			}

		}


		//record group[k] for group k, set k = k + 1
		pnumber += group[k].size();
		k++;
	}

	int mon = 0, dimer = 0, aggl = 0;
	double Mvol = 0., Dvol = 0., Avol = 0.;
	vector<int> asize;

	for(vector< vector<Nanoparticle> >::size_type i = 0 ; i != group.size() ; ++i) {

		if(group[i].size() == 1) {

			mon++;
			Mvol += group[i][0].particleVolume();
		}
		else if(group[i].size() == 2) {

			dimer++;
			for(size_t j = 0 ; j != group[i].size() ; ++j)
				Dvol += group[i][j].particleVolume();
		}
		else if(group[i].size() > 2) {

			aggl++;
			asize.push_back(group[i].size());
			for(size_t j = 0 ; j != group[i].size() ; ++j)
				Avol += group[i][j].particleVolume();

		}
	}

	out << Mvol/Tvol << " " << Dvol/Tvol << " " << Avol/Tvol << endl;

	int maxsize = 0;
	if(asize.size() != 0) {
		maxsize = asize[0];
		for(vector<int>::size_type i = 0 ; i != asize.size() ; ++i)
			maxsize = max(maxsize, asize[i]);
	}

	aggrsize << maxsize << endl;

}

double Metropolis::particleVolumeFraction(vector<Nanoparticle>& v)
{

	double sum = 0.;
	double n = v.size();
	for(size_t i = 0 ; i != n ; ++i)
		sum += v[i].particleVolume();

	return sum/pow(L,3);
}

double Metropolis::mrt2(double x)
{
    double a = -delta1;
    double b =  delta1;

    double dx = a + (b - a) * random();
    double xtrial = x + dx;

    if ( (xtrial < LEFT) || (xtrial > RIGHT) )
        return x;

    double r = logNormal(xtrial)/logNormal(x);

    if ( r > 1.  || r > random() )
    {
        x = xtrial;
        Naccpt++;
        return x;
    }
    else
    {
        return x;
    }

}


vector<Nanoparticle> Metropolis::mrt2(size_t index,
		const std::vector<Nanoparticle>& v, const vector<string>& interactions)
{

	vector<Nanoparticle> temp = v;
	Vec3 displacement = (2. * randomVector() - 1.) * delta1;
    Vec3 trial = temp[index].center + displacement;

	Vec3 mtrial = temp[index].dipole_moment + delta2 * randomVectorOnUnitSphere();
    mtrial.normalize();

	temp[index].center = trial;
	temp[index].dipole_moment = mtrial;


	if(!theyDontOverlap(index, temp))
    {
		return v;
	}


	Dipolar dip(Mo/(4. * PI));
	Vanderwaals van(-ALPHA/12.);
	Steric st(PI * Kb * grafting * TEMP/2.);
	Vec3 b_field(0., 1., 0.);

	int num_inter = interactions.size();
	double Eold = total_potential(index, v, interactions);
	double Enew = total_potential(index, temp, interactions);

	double dE = Enew - Eold;

	if (dE <= 0.) {

		Naccpt++;
		return temp;

	}
	else if (boltzmann(dE) >= random()) {

		Naccpt++;
		return temp;

	}
	else {

		return v;

	}

}

double totalZeemann(const vector<Nanoparticle>& v, const Vec3& b)
{
    double energy = 0.;
    size_t vsize = v.size();
    for( size_t i = 0 ; i < vsize ; ++i)
        energy += -dot(v[i].dipole_moment * v[i].dipole_length(), b);

    return energy;
}

double total_potential(size_t index, const std::vector<Nanoparticle>& v,
                       const std::vector<std::string>& interactions)
{
	Dipolar dip(Mo/(4. * PI));
	Vanderwaals van(-ALPHA/12.);
	Steric st(PI * Kb * grafting * TEMP/2.);
	Vec3 b_field(0., 1., 0.);

	int num_inter = interactions.size();
	double energy(0.);
	for(int i = 0 ; i < num_inter ; ++i)
	{
		if(interactions[i] == "dipolar")
			energy += dip.appliedTo(index, v);
		else if(interactions[i] == "vanderwaals")
			energy += van.appliedTo(index, v);
		else if(interactions[i] == "steric")
			energy += st.appliedTo(index, v);
		else if(interactions[i] == "zeeman")
			energy += -dot(v[index].dipole_moment * v[index].dipole_length(), b_field);
	}

	return energy;
}


double total_potential(const std::vector<Nanoparticle>& v,
                       const std::vector<std::string>& interactions)
{
	Dipolar dip(Mo/(4. * PI));
	Vanderwaals van(-ALPHA/12.);
	Steric st(PI * Kb * grafting * TEMP/2.);
	Vec3 b_field(0., 1., 0.);

	int num_inter = interactions.size();
	double energy(0.);
	for(int i = 0 ; i < num_inter ; ++i)
	{
		if(interactions[i] == "dipolar")
			energy += dip.totalInteraction(v);
		else if(interactions[i] == "vanderwaals")
			energy += van.totalInteraction(v);
		else if(interactions[i] == "steric")
			energy += st.totalInteraction(v);
		else if(interactions[i] == "zeeman")
			energy += totalZeemann(v, b_field);
	}

	return energy;
}