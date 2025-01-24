// 
// BHtree.C

#define PR(x)  cerr<<#x<<"="<<x<<" "
#define PRC(x) cerr<<#x<<"="<<x<<","
#define PRL(x) cerr<<#x<<"="<<x<<"\n"

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <chrono>

#define real double
#include "BHtree.h"

void bhnode::assign_root(myvector root_pos,real length,particle* p,int np){
	pos=root_pos;
	l=length;
	pfirst=p;
	nparticle=np;
	for(int i=0;i<np-1;i++){
		p->next=p+1;
		p++;
	}
	p->next=NULL;
}

int childindex(myvector pos,myvector cpos){
	int subindex=0;
	for(int k=0;k<3;k++){
		subindex<<=1;
		if(pos[k] > cpos[k]){
		subindex+=1;
		}
	}
	return subindex;
}

void bhnode::assign_child(int subindex,bhnode*& heap_top,int& heap_remainder){
	if(heap_remainder<=0){
		cerr<<"create_tree: no more free node... exit\n";
		exit(1);
	}
	if(child[subindex]==NULL){
		child[subindex]=heap_top;
		heap_top++;
		heap_remainder--;
		child[subindex]->cpos=cpos+myvector(((subindex&4)*0.5-1)*l/4,((subindex&2)-1)*l/4,((subindex&1)*2-1)*l/4);
		child[subindex]->l=l*0.5;
		child[subindex]->nparticle=0;
	}
}


void bhnode::create_tree_recursive(bhnode*& heap_top,int& heap_remainder){
	for(int i=0;i<8;i++)child[i]=NULL;
	particle* p=pfirst;
	for(int i=0;i<nparticle;i++){
		particle* pnext=p->next;
		int subindex=childindex(p->pos,cpos);
		assign_child(subindex,heap_top,heap_remainder);
		child[subindex]->nparticle++;
		p->next=child[subindex]->pfirst;
		child[subindex]->pfirst=p;
		p=pnext;
	}
	for(int i=0;i<8;i++) if(child[i]!=NULL){
		if(child[i]->nparticle>1){
			child[i]->create_tree_recursive(heap_top,heap_remainder);
		}else{
			child[i]->pos=child[i]->pfirst->pos;
			child[i]->mass=child[i]->pfirst->mass;
		}
	}
}

void spc(int indent){
	for(int i=0;i<indent;i++) cerr<< " ";
}

void bhnode::dump(int indent){
	int i;
	spc(indent);
	cerr<<"node center pos "<<cpos;
	cerr<<endl;
	spc(indent);cerr<<"node cm  "<<pos << " m " <<mass;
	if(nparticle == 1){
		cerr<< " IS LEAF" ;
		PRL(nparticle);
		particle* p=pfirst;
		for(i=0;i<nparticle;i++){
			for(int j=0;j<indent+2;j++) cerr<< " ";
			PRL(p->pos);
			p=p->next;
		}
	}else{
		cerr<< " IS _not_ LEAF ";
		PRL(nparticle);
		for(i=0;i<8;i++){
			if(child[i]!=NULL){
			child[i]->dump(indent+2);
			}
		}
	}
}

void bhnode::dumptree(int indent){
	int i;
	spc(indent); cerr<< "node center pos " << cpos ;
	if(nparticle == 1){
		cerr<< " IS LEAF" ;PRL(nparticle);
		particle* p=pfirst;
		for(i=0; i < nparticle; i++){
			for(int j=0;j<indent+2;j++)cerr<< " ";
			PRL(p->pos);
			p=p->next;
		}
	}else{
		cerr<< " IS _not_ LEAF ";PRL(nparticle);
		for(i=0;i<8;i++){
			if(child[i]!=NULL){
				child[i]->dumptree(indent+2);
			}
		}
	}
}

//cpos:center of the box
// position of the particle
// length of one side of the box
int inbox(myvector& cpos,myvector& pos,real l){
	for(int i=0;i<ndim;i++){
		if(fabs(pos[i]-cpos[i])>l*0.5) return 1;
	}
	return 0;
}

int bhnode::sanity_check(){
	int i;
	int iret=0;
	if(nparticle == 1 ){
		// this is the lowest level node. Things to check:
		// particle is in the cell
		particle* p=pfirst;
		if(inbox(pos,p->pos,l)){
			cerr<< "Error,particle out of box ... \n";
			dump();
			return 1;
		}
	}else{
		// This is the non-leaf node. Check the position and side
		// length of the child cells and then check recursively..
		for(i=0;i<8;i++){
			if(child[i] !=NULL){
				int err=0;
				err=child[i]->sanity_check();
				if(l*0.5 !=child[i]->l) err +=2;
				myvector relpos=cpos-child[i]->cpos;
				for (int k=0 ; k<ndim;k++){
					if(fabs(relpos[k]) !=l*0.25)err +=4;
				}
				if(err){
					cerr<< "Child " << i << " Error type=" << err << endl;
					dump();
				}
				iret +=err;
			}
		}
	}
	return iret;
}

void  bhnode::set_cm_quantities(){
	if(nparticle>1){
		int i;
		pos=0.0;
		mass=0.0;
		for(i=0;i<8;i++){
			if(child[i]!=NULL){
				child[i]->set_cm_quantities();
				real mchild=child[i]->mass;
				pos+=mchild*child[i]->pos;
				mass+=mchild;
			}
		}
		pos/=mass;
	}
}

void accumulate_force_from_point(myvector dx,real r2,real eps2,myvector& acc,real& phi,real jmass){
	double r2inv=1/(r2+eps2);
	double rinv=sqrt(r2inv);
	double r3inv=r2inv*rinv;
	phi -=jmass*rinv;
	acc +=jmass*r3inv*dx;
}

void bhnode::accumulate_force_from_tree(myvector& ipos,real eps2,real theta2,myvector& acc,real& phi){
	myvector dx=pos-ipos;
	real r2=dx*dx;
	if((r2*theta2>l*l)||(nparticle==1)){
		// node and position is well separated;
		accumulate_force_from_point(dx,r2,eps2,acc,phi,mass);
	}else{
		for(int i=0;i<8;i++){
			if(child[i]!=NULL){
				child[i]->accumulate_force_from_tree(ipos,eps2,theta2,acc,phi);
			}
		}
	}
}

void calculate_gravity_using_tree(particle* p,bhnode* bn,real eps2,real theta2){
	p->acc=0;
	p->phi=p->mass/sqrt(eps2);
	bn->accumulate_force_from_tree(p->pos,eps2,theta2,p->acc,p->phi);
}

const myvector uniform_random_position_in_sphere(double power_index){
	myvector x;
	do{
		for(int i=0; i<3;i++) x[i]=drand48()*2-1;
	}while (x*x >= 1);
	x*=pow(x*x,3.0/(power_index+3)-1);
	return x;
}

void create_uniform_sphere(particle* pb,int nbody,real power_index,real r0){
	PRC(nbody); PRL(power_index);
	particle* p=pb;
	for(int i=0; i<nbody;i++){
		p->pos=uniform_random_position_in_sphere(power_index)*r0;
		p->vel=0;
		p->mass=1.0/nbody;
		p++;
	}
}

void accumulate_mutual_gravity(particle& p1,particle& p2,real eps2){
	myvector dx=p1.pos-p2.pos;
	double r2inv=1/(dx*dx+eps2);
	double rinv=sqrt(r2inv);
	double r3inv=r2inv*rinv;
	p1.phi -=p2.mass*rinv;
	p2.phi -=p1.mass*rinv;
	p1.acc -=p2.mass*r3inv*dx;
	p2.acc +=p1.mass*r3inv*dx;
}

void clear_acc_and_phi(particle* parray,int n){
	particle* p=parray;
	for (int i=0;i<n;i++){
		p->acc=0;
		p->phi=0;
		p++;
	}
}

void calculate_uncorrected_gravity_direct(particle* parray,int n,double eps2){
	clear_acc_and_phi(parray,n);
	int i,j;
	particle* pi;
	particle* pj;
	for(i=0,pi=parray; i<n-1; i++,pi++){
		for(j=i+1,pj=pi+1; j<n; j++,pj++){
			accumulate_mutual_gravity(*pi,*pj,eps2);
		}
	}
}


#ifdef TEST
int main(){
	particle* pp;
	int n;
	cerr<<"Enter n:";
	cin>>n;
	pp=new particle[n];
	double rsize=1.0;
	create_uniform_sphere(pp,n,0,rsize);
	for(int i=0;i<n;i++){
		PRC(i);
		PRL(pp[i].pos);
	}
	bhnode* bn=NULL;

	int nnodes=n*2;
	bn=new bhnode[nnodes];
	bn->assign_root(myvector(0.0),rsize*2,pp,n);
	int heap_remainder=nnodes-1;
	bhnode* btmp=bn+1;
	bn->create_tree_recursive(btmp,heap_remainder);
	PRL(bn->sanity_check());
	bn->set_cm_quantities();
	bn->dump();
	double  eps2=0.01;
	calculate_uncorrected_gravity_direct(pp,n,eps2);
	cerr<< "Direct force \n";
	particle* p=pp;
	for(int i=0; i<n; i++){
		PR(i);
		PR(p->pos);
		PR(p->phi);
		PRL(p->acc);
		p++;
	}
	cerr<< "Tree   force \n";
	clear_acc_and_phi(pp,n);
	p=pp;
	for(int i=0; i<n; i++){
		calculate_gravity_using_tree(p,bn,eps2,0.5);
		PR(i);
		PR(p->pos);
		PR(p->phi);
		PRL(p->acc);
		p++;
	}
	return 0;
}
#endif

#ifdef SOTUKEN
int main(){
	chrono::steady_clock::time_point begin,end;
	chrono::microseconds time;

	particle* pp;
	int n;
	cerr<<"Enter n:";
	cin>>n;
	pp=new particle[n];
	double rsize=1.0;
	create_uniform_sphere(pp,n,0,rsize);
	bhnode* bn=NULL;
	int nnodes=n*2;
	double  eps2=0.01;
	bn=new bhnode[nnodes];
	for(int j=0;j<=5;j++){
		begin=chrono::steady_clock::now();
		bn->assign_root(myvector(0.0),rsize*2,pp,n);	//kore
		int heap_remainder=nnodes-1;
		bhnode* btmp=bn+1;
		bn->create_tree_recursive(btmp,heap_remainder);
		//PRL(bn->sanity_check());
		bn->set_cm_quantities();
	//	bn->dump();
	//	calculate_uncorrected_gravity_direct(pp,n,eps2);
	//	clear_acc_and_phi(pp,n);
		particle* p=pp;
		for(int i=0; i<n; i++){
			calculate_gravity_using_tree(p,bn,eps2,0.5);
			p++;
		}
		end=chrono::steady_clock::now();
		time=chrono::duration_cast<chrono::microseconds>(end-begin);
		cout<<j<<":"<<time.count()<<"μs"<<endl;
	}
/*
	for(int j=0;j<5;j++){
		clear_acc_and_phi(pp,n);
		p=pp;
		begin=chrono::steady_clock::now();
		for(int i=0;i<n;i++){
			calculate_gravity_using_tree(p,bn,eps2,0.5);
			p++;
		}
		end=chrono::steady_clock::now();
		time=chrono::duration_cast<chrono::microseconds>(end-begin);
		cout<<j+1<<":"<<time.count()<<"μs"<<endl;
	}
*/
	return 0;
}
#endif

#ifdef TREETEST
int main(){
	particle* pp;
	int n;
	cerr<<"Enter n:";
	cin>>n ;
	pp=new particle[n];
	double rsize=1.0;
	create_uniform_sphere(pp,n,0 ,rsize);
	for(int i=0;i<n;i++){
		PRC(i);
		PRL(pp[i].pos);
	}
	bhnode* bn=NULL;

	int nnodes=n*2;
	bn=new bhnode[nnodes];
	bn->assign_root(myvector(0.0),rsize*2,pp,n);
	int heap_remainder=nnodes-1;
	bhnode* btmp=bn+1;
	bn->create_tree_recursive(btmp,heap_remainder);
	PRL(bn->sanity_check());
	bn->dumptree();
	return 0;
}
#endif
