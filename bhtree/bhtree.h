#ifndef _BHTREE_H
#define _BHTREE_H

///////////////////////////////////
//BHtree:basic class for C++ implementation of BH treecode
//J.Makino 2001/1/23
//Modified 2005/1/16
//
//Copied&Changed by Inata Yusei
//Copied Date:2024/6/17
//Changed Date:2024/6/17
//Changed Point:real type->double type
///////////////////////////////////

#include "vector3D.h"
#include "particle.h"
#include <iostream>

class BHnode{
	private:
		Vector3D cpos;
		double l;
		BHnode* child;
		Particle* pfirst;
		int nparticle;
		Vector3D pos;
		double mass;
	//end private
	public:
		BHnode(){
			cpos=0.0;
			l=0.0;
			for(int i=0;i<8;i++){
				child[i]=NULL;
			}
			pfirst=NULL;
			nparticle=0;
			pos=0.0;
			mass=0.0;
		}

		void assignRoot(Vector3D root_pos,double length,Particle* p,int np);
		void create_tree_recursive(BHnode*& heap_top,int& heap_remainder);
		int childindex(Vector3D pos,Vector3D cpos);

	//end public
};

void BHnode::assign_root(Vector3D root_pos,double length,Particle* p,int np){
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

void BHnode::create_tree_recursive(BHnode*& heap_top,int& heap_remainder){
	for(int i=0;i<8;i++){
		child[i]=NULL;
	}
	Particle* p=this->pfirst;
	for(int i=0;i<this->nparticle;i++){
		Particle* pnext=p->next;
		int subindex=childindex(p->pos,cpos);
		assign_child(subindex,heap_top,heap_remainder);
		child[subindex]->nparticle++;
		p->next=child[subindex]->pfirst;
		child[subindex]->pfirst=p;
		p=pnext;
	}

	for (int i=0;i<8;i++){
		if(child[i]->nparticle>1){
			child[i]->create_tree_recursive(heap_top,heap_remainder);
		}else{
			child[i]->pos=child[i]->pfirst->pos;
			child[i]->mass=child[i]->pfirst->mass;
		}
	}
}

int BHnode::childindex(Vector3D pos,Vector3D cpos){
	int subindex=0;
	for(int k=0;k<3;k++){
		subindex<<=1;
		if(pos[k]>cpos[k]){
			subindex+=1;
		}
	}
	return subindex;
}

void BHnode::assign_child(int subindex,BHnode*& heap_top,int& heap_remainder){
	if(heap_remainder<=0){
		std::cerr<<"create_tree: no more free node... exit\n";
		exit(1);
	}
	if(child[subindex]==NULL){
		child[subindex]=heap_top;
		heap_top++;
		heap_reminder--;
		child[subindex]->cpos=cpos+myvector(((subindex&4)*0.5-1)*l/4,((subindex&2)-1)*l/4, ((subindex&1)*2-1)*l\4);
		child[subindex]->l=l*0.5;
		child[subindex]->nparticle=0;
	}
}

#endif //ifndef _BHTREE_H
