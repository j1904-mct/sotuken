#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "bhtree.h"

#define NO_MORE_FREE_NODE_ERROR -1

#define PR(x) std::cerr<<#x<<"="<<x<<" "
#define PRC(x) std::cerr<<#x<<"="<<x<<", "
#define PRL(x) std::cerr<<#x<<"="<<x<<"\n"

void BHnode:assign_root(Vector3D root_pos,double length,particle* p,int np){
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

int childindex(Vector3D pos,Vector3D cpos){
	int subindex=0;
	for(int i=0;i<3;i++){
		subindex<<=1;
		if(pos[i]>cpos[i]) subindex+=1;
	}
	return subindex;
}

void BHnode::assign_child(int subindex,bhnode*& heap_top,int& heap_remainder){
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
		child[subindex]->nparticle=0;
	}
}

void BHnode::create_tree_recursive(BHnode*& heap_top,int& heap_reminder){
	for(int i=0;i<8;i++) child[i]=NULL;
	particle* p=pfirst;
	for(int i=0;i<nparticle;i++){
		particle* pnext=p->next;
		int subindex=childindex(p->pos,cpos);
		assign_child(subindex,heap_top,heap_reminder);
		child[subindex]->nparticle++;
		p->next=child[subindex]->pfirst;
		child[subindex]->pfirst=p;
		p->next;
	}
	for(int i=0;i<8;i++){
		if(child[i]!=NULL){
			if(child[i]->nparticle>1)
