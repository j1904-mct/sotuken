#ifndef  BHTREE_H
#  define   BHTREE_H
/*-----------------------------------------------------------------------------
 *  BHtree : basic class for C++ implementation of BH treecode
 *  J. Makino 2001/1/23
 *  Modified  2005/1/16
 *-----------------------------------------------------------------------------
 */


#include "vector.h"
#include "particle.h"

class bhnode
{
private:
    myvector cpos;
    real l;
    bhnode * child[8];
    particle * pfirst;
    int nparticle;
    myvector pos;
    real mass;
    
public:
    bhnode(){
	cpos = 0.0;
	l = 0.0;
	for(int i = 0; i<8;i++)child[i] = NULL;
	pfirst = NULL;
	nparticle = 0;
	pos = 0.0;
	mass = 0.0;
    }
    void assign_child(int subindex, bhnode * & heap_top, int & heap_remainder);
    void create_tree_recursive(bhnode * & heap_top, int & heap_remainder);
    void assign_root(myvector root_pos, real length, particle * p, int nparticle);
    void dump(int indent = 0);
    void dumptree(int indent = 0);
    int sanity_check();

    void set_cm_quantities();
    void accumulate_force_from_tree(myvector & ipos, real eps2, real theta2,
				   myvector & acc,
				   real & phi);
};

void calculate_gravity_using_tree(particle * p, bhnode * bn, real eps2, real theta2);

#endif
