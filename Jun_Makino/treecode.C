// 
// treecode.C

#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

//#include  "../cpgplot.h"
#include  <stdlib.h>
#include  <math.h>
#include  <iostream>

#define real double
#include "BHtree.h"

/*void initgraph()
{
    
    if(cpgopen("") !=1 ) exit(-1);
    cpgask(0);
}*/

/*void plot_particles(particle * p,
		    int n,
		    double xmax,
		    int first)
{
    static float *  xp;
    static float *  yp;
    static int nbuf = -1;
    cpgbbuf();
    if (first == 1){
         cpgpage();
         cpgenv(-xmax, xmax, -xmax, xmax,1,-2);
    }
    if ( nbuf < n){
	if (nbuf > 0){
	    delete [] xp;
	    delete [] yp;
	}
	xp = new float[n];	
	yp = new float[n];
	nbuf = n;
    }
    for(int i = 0;i<n;i++){
	xp[i] = p[i].pos[0];
	yp[i] = p[i].pos[1];
    }
    cpgeras();
    cpgbox("bcnst", 0.0, 0, "bcnst", 0.0, 0);
    cpglab("X", "Y", " ");
    cpgpt(n,xp,yp,-1);
    cpgebuf();
}*/


real calculate_size(particle * p,
		     int n)
{
    real rsize = 1;
    for(int i=0;i<n;i++){
	myvector ppos = p[i].pos;
	for(int k=0;k<3;k++) {
	    if (fabs(ppos[k])>rsize) {
		rsize *= 2;
	    }
	}
    }
    return rsize;
}

void calculate_gravity(bhnode* bn,
		       int nnodes,
		       particle * pp,
		       int n,
		       real eps2,
		       real theta)
{
    real rsize = calculate_size(pp,n);
    bn->assign_root(myvector(0.0), rsize*2, pp, n);
    int heap_remainder = nnodes-1;
    bhnode * btmp = bn + 1;
    bn->create_tree_recursive(btmp,  heap_remainder);
    //    bn->dump();
    //    PRL(bn->sanity_check());
    bn->set_cm_quantities();
    clear_acc_and_phi(pp, n);
    
    particle * p = pp;
    for(int i = 0; i<n; i++){
	calculate_gravity_using_tree(p, bn, eps2, theta);
	//	PR(i); PR(p->pos);PR(p->phi); PRL(p->acc);
	p++;
    }
}


void integrate(bhnode * bn,
	       int nnodes,
	       particle * pp,
	       int n,
	       real eps2,
	       real theta,
	       real dt)
{
    for(int i = 0;i<n;i++) pp[i].predict(dt);
    calculate_gravity(bn,nnodes,pp,n,eps2,theta);
    for(int i = 0;i<n;i++) pp[i].correct(dt);
}    

real kinetic_energy(particle * pp,
		    int n)
{
    real ke = 0;
    for(int i = 0;i<n;i++) {
	ke += (pp[i].vel*pp[i].vel)*pp[i].mass*0.5;
    }
    return ke;
}

void print_cm(particle * pp,
	      int n)
{
    myvector cmpos = 0.0;
    myvector cmvel = 0.0;
    real cmmass = 0.0;
    for(int i = 0;i<n;i++) {
	cmpos += pp[i].mass*pp[i].pos;
	cmvel += pp[i].mass*pp[i].vel;
	cmmass += pp[i].mass;
    }
    cmpos /= cmmass;
    cmvel /= cmmass;
    cerr << "CM pos = " << cmpos << " vel = " << cmvel <<endl;
}

real potential_energy(particle * pp,
		      int n)
{
    real pe = 0;
    for(int i = 0;i<n;i++) {
	pe += pp[i].phi*pp[i].mass;
    }
    return pe*0.5;
}

void print_energies(particle * pp,
		    int n,
		    int check)
{
    static real etot0;
    real pe = potential_energy(pp,n); 
    real ke = kinetic_energy(pp,n);
    real etot = pe+ke;
    if (check==0 ) etot0 = etot;
    cerr << "etot = " << pe+ke << " PE=" << pe << " KE="<<ke << " V.R.= " << -ke/pe;
    if(check)cerr <<  " de=" << (etot-etot0)/etot0;
    cerr <<endl;
}
    


int main()
{
    particle * pp;
    int n;
    cerr << "Enter n:";
    cin >> n ;
    pp = new particle[n];
    double rsize = 1.0;
    cerr << "Enter power index:";
    real power_index;
    cin >> power_index ;
    cerr << "power index = " << power_index <<endl;
    create_uniform_sphere(pp, n, power_index , rsize);
    //    pp[0].pos=myvector(1,0,0);
    //
    //    pp[1].pos=myvector(-1,0,0);
    //    pp[0].vel=myvector(0,0.5,0);
    //    pp[1].vel=myvector(0,-0.5,0);
    //    for(int i =0;i <n;i++){
    //PRC(i); PRC(pp[i].pos);
    //    }
    bhnode * bn = NULL;

    int nnodes = n*2+100;
    bn = new bhnode[nnodes];
    real eps2;
    real theta;
    real dt;
    real tend;
    int iout;
    cerr << "Enter eps2, theta, dt, tend, iout:";
    cin >> eps2 >>theta >>dt >>tend >> iout;
    cerr << "eps2=" << eps2 << " theta=" <<theta <<" dt=" << dt 
	 << " tend=" << tend << " iout=" << iout <<endl;
    real xmax_plot;
    cerr << "Enter xmax for plot:";
    cin >> xmax_plot;
    cerr << "xmax for plot = " << xmax_plot <<endl;
    //initgraph();
    //plot_particles(pp,n,xmax_plot,1);
    calculate_gravity(bn,nnodes,pp,n,eps2, theta);
    cerr << "Initial data"<<endl;
    print_energies(pp,n,0);
    int istep = 0;
    for(real t=0;t<tend; t+=dt){
	integrate(bn, nnodes, pp, n, eps2, theta, dt);
	istep++;
	//plot_particles(pp,n,xmax_plot,0);
	if (istep % iout == 0){
	    cerr << "Time=" << t+dt <<endl;
	    print_cm(pp, n);
	    print_energies(pp,n,1);
	}
    }
    return 0;
}
