///////////////////////////
//calculate center of gravity	with Vector class
//Created:Inata Yusei
//Date:2024/6/17
///////////////////////////

#include <iostream>
#include <iosfwd>

#define VLEN 3
#include "vector.h"

int main(void){
	Vector v1(0.0,0.0,0.0);
	Vector v2(1.0,1.0,1.0);
	Vector v3(2.0,2.0,2.0);
	Vector result;
//	std::cout<<v1;

	result=grvPtTop(v1,v2);
	std::cout<<"center of gravity about vec1,vec2 \n";
	std::cout<<result;

	result=grvPtTop(v1,v2,v3);
	std::cout<<"center of gravity about vec1,vec2,vec3 \n";
	std::cout<<result;
	return 0;
}
