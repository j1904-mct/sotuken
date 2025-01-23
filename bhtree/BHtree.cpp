#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "bhtree.h"

#define NO_MORE_FREE_NODE_ERROR -1

#define PR(x) std::cout<<#x<<"="<<x<<" "
#define PRC(x) std::cout<<#x<<"="<<x<<", "
#define PRL(x) std::cout<<#x<<"="<<x<<"\n"

void BHnode::assign_root(Vector3D root_pos,double length,Particle* p,int np){
	pos=root_pos;
	l=length;
	pfirst=p;
	nParticle=np;
	for(int i=0;i<np-1;i++){
		p->next=p+1;
		p++;
	}
	p->next=NULL;
}

int childindex(Vector3D pos,Vector3D cpos){
	int subindex=0;
	for(int i=0;i<3;i++){
		subindex<<=1;
		if(pos[i]>cpos[i]) subindex+=1;
	}
	return subindex;
}

void BHnode::assign_child(int subindex,BHnode*& heap_top,int& heap_remainder){
	if(heam_reminder<=0){
		std:cerr<<"create tree:no more free node... exit"<<std::endl;
		exit(NO_MORE_FREE_NODE_ERROR);
	}
	if (child[subindex]==NULL){
		child[subindex]=heap_top;
		heap_top++;
		heap_reminder--;
		child[subindex]->cpos=cpos+Vector3D(((subindex&4)*0.5-1)*l/4,((subindex&2)-1)*l/4,((subindex&1)*2-1)*l/4);
		child[subindex]->l=l*0.5;
		child[subindex]->nParticle=0;
	}
}

void BHnode::create_tree_recursive(BHnode*& heap_top,int& heap_reminder){
	for(int i=0;i<8;i++) child[i]=NULL;
	Particle* p=pfirst;
	for(int i=0;i<nParticle;i++){
		Particle* pnext=p->next;
		int subindex=childindex(p->pos,cpos);
		assign_child(subindex,heap_top,heap_reminder);
		child[subindex]->nParticle++;
		p->next=child[subindex]->pfirst;
		child[subindex]->pfirst=p;
		p->next;
	}
	for(int i=0;i<8;i++){
		if(child[i]!=NULL){
			if(child[i]->nParticle>1){
				child[i]->create_tree_recursive(heap_top,heap_remainder);
			}else{
				child[i]->pos=child[i]->pfirst->pos;
				child[i]->mass=child[i]->pfist->mass;
			}
		}
	}
}

void spc(int indent){
for(int i=0;i<indent;i++) cout<<" ";
}

void BHnode::dump(int indent){
	int i,j;
	i=j=0;
	spc(indent);
	std::cout<<"node center pos "<<cpos<<std::endl;
	spc(indent);
	std::cout<<"node cm "<<pos<<" m "<<mass;
	if(nParticle==1){
		std::cout<<"IS LEAF";
		PRL(nParticle);
		Particle* p=pfirst;
		for(i=0;i<nParticle;i++){
			for(j=0;j<indent+2;j++) std::cout<<" ";
			PRL(p->pos);
			p=p->next;
		}
	}else{
		std::cout<<"IS _not_ LEAF";
		PRL(nParticle);
		for(int i=0;i<8;i++){
			if(child[i]!=NULL) child[i]->dump(indent+2);
		}
	}
}

void BHnode::dumptree(int indent){
	int i,j;
	i=j=0;
	spc(indent);
	std::cout<<"node center pos "<<cpos;
	if(nParticle==1){
		std::cout<<" IS LEAF";
		PRL(nParticle);
		Particle* p=pfirst;
		for(i=0;i<nParticle;i++){
			for(j=0;j<indent+2;j++) std::cout<<" ";
			PRL(p->pos);
			p=p->next;
		}
	}else{
		sttd::cout<<" IS _not_ LEAF";
		PRL(nParticle);
		for(i=0;i<8;i++){
			if(child[i]!=NULL) child[i]->dumptree(indent+2);
		}
	}
}
//cpos:center of box
//pos:position of the Particle
//l:length of one side of the box
int inbox(Vector3D& cpos,Vector3D& pos,double l){
	for(int i=0;i<ndim;i++){
		if(fabs(pos[i]-cpos[i])>l*0.5)return 1;
	}
	return 0;
}

int BHnode::sanity_check(){
	int i;
	int iret=0;
	if(nParticle==1){
		//Leaf node (Lowest level node)
		//Check:Particle is in the cell
		Particle* p->pfirst;
		if(inbox(pos,p->pos,l)){
			std::cerr<<"Error,Particle out of box..."<<std::endl;
			dump();
			return 1;
		}
	}else{
		//Non-leaf node
		//Check:position,side length of the child cells,recursively
		for(int i=0;i<8;i++){
		if(child[i]!=NULL){
			int err=0;
			err=child[i]->sanity_check;
			if(l*0.5!=child[i]->l) err+=2;
			Vector3D relpos=cpos-child[i]->cpos;
			for(int k=0;k<ndim;k++){
				if(fabs(relpos[k])!=l*0.25) err+=4;
			}
			if(err){
				std::cerr<"Child "<<i<<" Error type="<<err<<endl;
				dump();
			}
			iret+=err;
		}
	}
	return iret;
}

void BHnode::set_cm_quantities(){
	if(nParticle>1){
		int i;
		pos=0.0;
		mass=0.0;
		for(int i=0;i<8;i++){
			if(child[i]!=NULL){
				child[i]->set_cm_quantities();
				double mchild=child[i]->mass;
				pos+=mchild*child[i]->pos;
				mass+=mchild;
			}
		}
		pos/=mass;
	}
}

void accumulate_force_from_point(Vector3D dx,double r2,real eps2,Vector3D& acc,double& phi,double jmass){
	double r2inv=1/(r2+eps2);
	double rinv=sqrt(r2inv);
	double r3inv=r2inv*rinv;
	phi-=jmass+rinv;
	acc+=jmass*r3inv*dx;
}

void accumulate_force_from_tree(Vector3D& ipos,double eps2,double theta2,Vector3D& acc,double& phi){
	Vector3D dx=pos-ipos;
	double r2=dx*dx;
	if((r2*theta2>l*l)||(nParticle==1)){ //node and position is well separated;
		accumulate_force_from_point(dx,r2,eps2,acc,phi,mass);
	}else{
		for(int i=0;i<8;i++){
			if(child[i]!=NULL) child[i]->accumulate_force_from_tree(ipos,eps2,theta2,acc,phi);
		}
	}
}

void calculate_gravity_using_tree(Particle* p,BHnode* bn,double eps2,double theta2){
	p->acc=0;
	p->phi=p->mass/sqrt(eps2);
	bn->accumulate_force_from_tree(p->pos,eps2,theta2,p->acc,p->phi);
}

const Vector3D uniform_random_position_in_sphere(double power_index){
	Vector3D x;
	int i;
	do{
		for(i=0;i<3;i++) x[i]=drand48()*2-1;
	}while(x*x>=1);
	x*=pow(x*x,3.0/(power_index+3)-1);
	return x;
}

void create_uniform_sphere(Particle* pb,int nbody,double power_index,double r0){
	PRC(nbody);
	PRL(power_index);
	Particle* p=pb;
	for(int i=0;i<nbody;i++){
		p->pos=uniform_random_position_in_sphere(power_index)*r0;
		p->vel=0;
		p->mass=1.0/nbody;
		p++;
	}
}

void accumulate_mutual_gravity(Particle& p1,Particle& p2,double eps2){
	Vector3D dx=p1.pos-p2.pos;
	double r2inv=1/(dx*dx+eps2);
	double rinv=sqrt(r2inv);
	double r3inv=r2inv*rinv;
	p1.phi-=p2.mass*rinv;
	p2.phi-=p1.mass*rinv;
	p1.acc-=p2.mass*r3inv*dx;
	p2.acc+=p1.mass*r3inv*dx;
}

void clear_acc_and_phi(parray,n);
	int i,j;
	Particle* pi;
	Particle* pj;
	for(i=0,pi=parray;i<n-1;i++,pi++){
		for(j=i+1,pj=pi+1;j<n;j++,pj++){
			accumulate_mutual_gravity(*pi,*pj,eps2);
		}
	}
}

#ifdef TEST
int main(void){
	Particle* pp;
	int n;
	std::cout<<"Enter n:";
	std::cin>>n;
	pp=new Particle[n];
	double rsize=1.0;
	create uniform_sphere(pp,n,0,rsize);
	for(int i=0;i<n;i++){
		PRC(i);
		PRL(pp[i].pos);
	}
	BHnode*bn=NULL;

	int nodes=n*2;
	bn=new BHnode[nnodes];
	bn->assign_root(Vector3D(0.0),rsize*2,pp,n);
	int heap_remainder=nnodes-1;
	BHnode* btmp=bn+1;
	bn->create_tree_recursive(btmp,heap_remainder);
	PRL(bn->sanity_check());
	bn->set_cm_quantities();
	bn->dump;
	double eps2=0.01;
	calculate_uncorrected_gravity_direct(pp,n,eps2);
	std::cout<<"Direct force"<<std::endl;
	Particle* p=pp;
	int i;
	for(i=0;i<n;i++){
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


