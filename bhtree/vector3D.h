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

/////////////////////////////////////////////////////////
//Copied from https://jun-makino.sakuraweb.com/kougi/keisan_tenmongakuII/programs/simpletree/Vector3D.h
//Copied&Changed by Inata Yusei
//Copied at 2024/6/21
//Changed Point: add comment,change this to it easier to read,class name (vector -> Vector3D),add constant DIM(means number of dimensions)
//
/////////////////////////////////////////////////////////

#ifndef  Vector3D_H
#define Vector3D_H
//using namespace std;

typedef std::ostream ostream;
typedef std::istream istream;

//#include "stdinc.h"

/*-----------------------------------------------------------------------------
 *  Vector3D  --  a class for 3-dimensional Vector3Ds
 *-----------------------------------------------------------------------------
 */

#define DIM 3;
const int ndim = 3;

class Vector3D{
	private:
		double element[DIM];
	//private end
	public:
		//Constructor
		//Default: initialize to zero.
		Vector3D(double c = 0){
			element[0] = element[1] = element[2] = c;
		}
		Vector3D(double x, double y, double z){
			element[0] = x; element[1] = y; element[2] = z;
		}
		Vector3D(double data[DIM]){
			for(int i=0;i<DIM;i++){
				this->element[i]=data[i];
			}
		}
		template<typename... Data>
		Vector3D(Data... data){
			int i=0;
			for(double d: std:initializer_list<double>{data...}){
				this->element[i]=d;
				i++;
			}
		}
		//Constructor end

		//  []: the return type is declared as a reference (&), so that it can be used
		//  on the left-hand side of an asignment, as well as on the right-hand side,
		//  i.e.  v[1] = 3.14  and  x = v[2]  are both allowed and work as expected.
		double& operator [] (int i){
			return element[i];
		}
		inline void print(){
			//std::cout << element[0] << " " << element[1] << " "<< element[2] << "\n";
			std::cout<<"["<<element[0]<<","<<element[1]<<","<<element[2]<<"]"<<std::endl;
		}

		#ifdef PERIODIC
		// PERIODIC basic reajustment menber function
		Vector3D readjust(){
			return Vector3D(readjust_r(element[0]),readjust_r(element[1]),readjust_r(element[2]));
		}
		#else
		Vector3D  readjust(){
			return Vector3D(*this);
		}
		#endif //ifdef PERIODIC



		//	Unary -
		Vector3D operator - (){
			return Vector3D(-element[0], -element[1], -element[2]);
		}

		//	Dot product.

		double operator * (const Vector3D& b){
			return element[0]*b.element[0]+element[1]*b.element[1]+element[2]*b.element[2];
		}

		//	Outer product.

		Vector3D operator ^ (const Vector3D& b){
//			return Vector3D(element[1]*b.element[2]-element[2]*b.element[1],element[2]*b.element[0]-element[0]*b.element[2],element[0]*b.element[1]-element[1]*b.element[0]);
			Vector3D ret;
			ret.x=this->element[1]*b.element[2]-this->element[2]*b.element[1];
			ret.y=this->element[2]*b.element[0]-this->element[0]*b.element[2];
			ret.z=this->element[0]*b.element[1]-this->element[1]*b.element[0];
			return ret;
		}

		//	Vector +, -

		Vector3D operator + (const Vector3D &b){
			return Vector3D(element[0]+b.element[0],element[1]+b.element[1],element[2]+b.element[2]);
		}
		Vector3D operator - (const Vector3D &b){
			return Vector3D(element[0]-b.element[0],element[1]-b.element[1],element[2]-b.element[2]);
		}
		friend Vector3D operator + (double, const Vector3D & );
		friend Vector3D operator + (const Vector3D &, double);

		//	Scalar *, /
		friend Vector3D operator * (double, const Vector3D & );
		friend Vector3D operator * (const Vector3D &, double);
		friend Vector3D operator / (const Vector3D &, double);

		//	Vector3D +=, -=, *=, /=
		Vector3D& operator += (const Vector3D& b)	{
			this->element[0] += b.element[0];
			this->element[1] += b.element[1];
			this->element[2] += b.element[2];
			return *this;
		}

		Vector3D& operator -= (const Vector3D& b){
			this->element[0] -= b.element[0];
			this->element[1] -= b.element[1];
			this->element[2] -= b.element[2];
			return *this;
		}

		Vector3D& operator *= (const double b){
			this->element[0] *= b;
			this->element[1] *= b;
			this->element[2] *= b;
			return *this;
		}

		Vector3D& operator /= (const double b){
		register double binv = 1.0/b;
			this->element[0] *= binv;
			this->element[1] *= binv;
			this->element[2] *= binv;
			return *this;
		}

		//      Input / Output
		friend ostream& operator << (ostream& s,const Vector3D& v);
		friend istream& operator >> (istream& s,Vector3D& v);

	//public end
};

//Input/Output
inline ostream& operator << (ostream& s, const Vector3D & v){
//	return s << v.element[0] << "  " << v.element[1]<< "  " << v.element[2];
	return s<<"["<<v.element[0]<<","<<v.element[1]<<","<<v.element[2]<<"]"<<std::endl;
}

inline istream & operator >> (istream & s, Vector3D & v){
	s >> v.element[0] >> v.element[1] >> v.element[2];
	return s;
}

inline double square(Vector3D v) {return v*v;}
inline double abs(Vector3D v)    {return sqrt(v*v);}

//Vector +
inline Vector3D operator + (double b, const Vector3D & v){
	return Vector3D(b+v.element[0],b+v.element[1],b+v.element[2]);
	}

inline  Vector3D operator + (const Vector3D & v, double b){
	return Vector3D(b+v.element[0],b+v.element[1],b+v.element[2]);
	}

//Scalar *,/
inline  Vector3D operator * (double b, const Vector3D & v){
	return Vector3D(b*v.element[0],b*v.element[1],b*v.element[2]);
}

inline  Vector3D operator * (const Vector3D & v, double b){
	return Vector3D(b*v.element[0],b*v.element[1],b*v.element[2]);
}

inline  Vector3D operator / (const Vector3D & v, double b){
	return Vector3D(v.element[0]/b,v.element[1]/b,v.element[2]/b);
}

#endif

//=======================================================================
//  |  the end of:  |         /|\         |  inc/vector.h
//========================= STARLAB =====================================
