#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "vector3D.h"
////////////////////////////
//nbody-particle:basic class for simple nbody implementation
//J.Makino 1998/11/29
//Modified 2005/01/16
//
//Copied&Changed by Inata Yusei
//Copy date:2024/06/17
//Change Point:real type -> double type
/////////////////////////////

class Particle{
	public:
		Vector3D pos;	//position vector
		Vector3D vel;	//velocity vector
		Vector3D acc; //accelerate vector
		double phi;		//phisical
		double mass;	//mass
		Particle* next;	//next particle pointer
		Particle* subnext;	//

		Particle(){
			pos=0.0;
			vel=0.0;
			acc=0.0;
			phi=mass=0.0;
		}

		void predict(double dt){
			double dt2=dt*dt*0.5;
			pos+=(dt*vel+dt2*acc);
			vel+=(dt*0.5)*acc;
		}

		void correct(double dt){
			vel+=(dt*0.5)*acc;
		}
	//end public
};

#endif //ifndef _PARTICLE_H
