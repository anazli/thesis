#include<iostream>
#include "../include/particle.h"
#include "../include/vec3.h"

using namespace std;

int main()
{
    Particle p(Vec3(0.,1.,0.), 15*1E-9);
    Particle p1(p); 

    return 0;
}