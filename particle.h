#ifndef  PARTICLE_H
#  define   PARTICLE_H
/*-----------------------------------------------------------------------------
 *  nbody-particle : basic class for simple nbody implementation
 *  J. Makino 1998/11/29
 *  Modified  2005/01/16
 *-----------------------------------------------------------------------------
 */
#include "vector.h"

class particle
    {
    public:
        myvector pos;
        myvector vel;
	myvector acc;
	real phi;
	real mass;
        particle *  next;
	particle(){
	    pos = 0.0;
	    vel = 0.0;
	    acc = 0.0;
	    phi =  mass =  0.0;
	}
	void predict(real dt){
	    real dt2 = dt*dt*0.5;
	    pos = pos + dt*vel + dt2*acc;
	    vel +=  (dt*0.5)*acc;
	}
	void correct(real dt){
	    vel +=  (dt*0.5)*acc;
	}
};

void create_uniform_sphere(particle * pb, int nbody, real power_index, real r0);
void clear_acc_and_phi(particle * parray, int n);


#endif
