/*
*  vector.h: 3D vector operations include file
*.............................................................................
*    version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
*    version 2:
*.............................................................................
*     This file includes:
*  1) definition of class vector
*.............................................................................
*/

// This is slightly modified version of vector header file
// taken from STARLAB
// -- J. Makino

///////////////////////////////////////////////////////////
//Copied from https://jun-makino.sakuraweb.com/kougi/keisan_tenmongakuII/programs/simpletree/vector.h
//Copied&Changed:Inata Yusei
//Changed at 2024/06/19
//Changed point: changed this to make it easier to read
//               fixed to prevent compilation errors
///////////////////////////////////////////////////////////

#ifndef  STARLAB_VECTOR_H
#define  STARLAB_VECTOR_H
using namespace std;

#define vector myvector
//#include "stdinc.h"

/*-----------------------------------------------------------------------------
*  vector  --  a class for 3-dimensional vectors
*-----------------------------------------------------------------------------
*/

const int ndim = 3;

class vector{
	private:
		double element[3];
	//end private					//added comment by Inata Yusei
	public:
	//	Default: initialize to zero.
		vector(double c = 0){
			element[0] = element[1] = element[2] = c;
		}

		vector(double x, double y, double z){
			element[0] = x;
			element[1] = y;
			element[2] = z;
		}

		//  []: the return type is declared as a reference (&), so that it can be used
		//  on the left-hand side of an asignment, as well as on the right-hand side,
		//  i.e.  v[1] = 3.14  and  x = v[2]  are both allowed and work as expected.

		double & operator [] (int i) {
			return element[i];
		}

		inline void print() {
		cout << element[0] << " " << element[1] << " "<< element[2] << std::endl;
		}
#ifdef PERIODIC
// PERIODIC basic reajustment menber function
		vector readjust(){
			return vector(readjust_r(element[0]),readjust_r(element[1]),readjust_r(element[2]));
		}
#else
		vector  readjust(){
			return vector(*this);
		}
#endif

		//	Unary -

		vector operator - (){
			return vector(-element[0], -element[1], -element[2]);
		}

		//	Dot product.
		double operator * (const vector& b){
			return element[0]*b.element[0]+element[1]*b.element[1]+element[2]*b.element[2];
		}

		//	Outer product.
		vector operator ^ (const vector &b){
			return vector(element[1]*b.element[2] - element[2]*b.element[1],element[2]*b.element[0] - element[0]*b.element[2],element[0]*b.element[1] - element[1]*b.element[0]);
		}

		//	Vector +, -
		vector operator + (const vector &b){
			return vector(element[0]+b.element[0],element[1]+b.element[1],element[2]+b.element[2]);
		}
		vector operator - (const vector &b){
			return vector(element[0]-b.element[0],element[1]-b.element[1],element[2]-b.element[2]);
		}

		friend vector operator + (double, const vector & );
		friend vector operator + (const vector &, double);

		//	Scalar *, /
		friend vector operator * (double, const vector & );
		friend vector operator * (const vector &, double);
		friend vector operator / (const vector &, double);

		//	Vector +=, -=, *=, /=

		vector& operator += (const vector& b){
			element[0] += b.element[0];
			element[1] += b.element[1];
			element[2] += b.element[2];
			return *this;
		}

		vector& operator -= (const vector& b){
			element[0] -= b.element[0];
			element[1] -= b.element[1];
			element[2] -= b.element[2];
			return *this;
		}

		vector& operator *= (const double b){
			element[0] *= b;
			element[1] *= b;
			element[2] *= b;
			return *this;
		}

		vector& operator /= (const double b){
			double binv = 1.0/b;	//C++17:removed register keyword 	//Changed by Inata Yusei 
			element[0] *= binv;
			element[1] *= binv;
			element[2] *= binv;
			return *this;
		}

		//      Input / Output
		friend ostream & operator << (ostream & , const vector & );
		friend istream & operator >> (istream & , vector & );

	//end public			//added by Inata Yusei
};

inline  ostream & operator << (ostream & s, const vector & v){
	return s << v.element[0] << "  " << v.element[1]<< "  " << v.element[2];
}

inline  istream & operator >> (istream & s, vector & v){
	s >> v.element[0] >> v.element[1] >> v.element[2];
	return s;
}

inline  double square(vector v){
	return v*v;
}
inline  double abs(vector v){
	return sqrt(v*v);
}

inline  vector operator + (double b, const vector & v){
	return vector(b+v.element[0],b+v.element[1],b+v.element[2]);
}

inline  vector operator + (const vector & v, double b){
	return vector(b+v.element[0],b+v.element[1],b+v.element[2]);
}

inline  vector operator * (double b, const vector & v){
	return vector(b*v.element[0],b*v.element[1],b*v.element[2]);
}

inline  vector operator * (const vector & v, double b){
	return vector(b*v.element[0],b*v.element[1],b*v.element[2]);
}

inline  vector operator / (const vector & v, double b){
	return vector(v.element[0]/b,v.element[1]/b,v.element[2]/b);
}

#endif

//=======================================================================
//  |  the end of:  |         /|\         |  inc/vector.h
//========================= STARLAB =====================================

