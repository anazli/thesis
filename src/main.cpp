#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>

#include "../include/particle.h"
#include "../include/vec3.h"
#include "../include/nanoparticle.h"
#include "../include/parameters.h"
#include "../include/interactions.h"
#include "../include/metropolis.h"
#include "../include/random.h"

using namespace std;

const size_t Nparticles = 100;

double totalZeemann(const vector<Nanoparticle>&, const Vec3&);
Vec3 systemMagnetization(const vector<Nanoparticle>&);


int main()
{
    vector<Nanoparticle> v;

    Metropolis init_sim(1.5, 0.2, 1000, 1984, Nparticles-1);

    // Setting diameters
    double diam = LEFT;
    for ( size_t i = 0 ; i < Nparticles ; ++i )
    {
        diam = init_sim.mrt2(diam);
        Vec3 cen = randomVector() * NM;
        Vec3 dip =  randomVectorOnUnitSphere();
        Vec3 easy = Vec3(0.,0.,1.);//randomVectorOnUnitSphere();
        Nanoparticle np(cen, diam * NM, dip, easy, 1.2 * NM);
        v.push_back(np);
    }

    ofstream energy_stream, magnet_stream,
             coords_stream, dipolCoords_stream,
             aglomerates_stream, aglSize_stream;

    energy_stream.open("data/energy.dat");
    magnet_stream.open("data/magnetization.dat");
    coords_stream.open("data/coords.dat");
    dipolCoords_stream.open("data/m_coords.dat");
    aglomerates_stream.open("data/aglomerates.dat");
    aglSize_stream.open("data/aglomerates_size.dat");


    init_sim.setPositions(v);

    double norm = init_sim.particleVolumeFraction(v) * MSAT;
    Vec3 magnetization = systemMagnetization(v);
    cout << "Initial magnetization of the system:" << magnetization/norm << endl;
    vector<string> interactions{"zeeman"};//{"dipolar", "anisotropy", "zeeman"};
    Metropolis main_sim(1.* NM, 0.2, 500000, 999, Nparticles-1);
    double energy = total_potential(v, interactions); 
    energy_stream << energy << endl;
    magnet_stream << magnetization.z()/norm << endl;

    for ( size_t t = 0 ; t < main_sim.steps() ; ++t )
    {
        size_t index = main_sim.randomInt();
        v = main_sim.mrt2(index, v, interactions);

        if ( t % v.size() == 0 )
        {
            energy = total_potential(v, interactions);
            magnetization = systemMagnetization(v);

            energy_stream << energy << endl;
            magnet_stream << magnetization.z()/norm << endl;
        }
    }

    cout << "Final magnetization of the system:" << magnetization/norm << endl;

    for ( size_t i = 0 ; i < Nparticles ; ++i )
    {
        coords_stream << v[i].center << endl;
        dipolCoords_stream << v[i].dipole_moment << endl;
    }

    energy_stream.close();
    magnet_stream.close();
    coords_stream.close();
    dipolCoords_stream.close();
    aglomerates_stream.close();
    aglSize_stream.close();

    return 0;
}


Vec3 systemMagnetization(const vector<Nanoparticle>& v)
{

	Vec3 mag(0.,0.,0.);
	size_t n = v.size();
	size_t i;
	for( i = 0 ; i < n ; ++i ) 
        mag = mag + v[i].dipole_moment * v[i].dipole_length();

	mag = mag/(L*L*L);
	return mag;
}
