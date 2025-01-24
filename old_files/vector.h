#ifndef VECTOR_H
#define VECTOR_H

/////////////////////
//Vector: a class for VLEN-dimensional vectors
//VLEN must be defined as constant before this header
//created by MAKINO Junichiro
//copy & make change by INATA Yusei
//copy date:2024/6/13
//copy from https://jun-makino.sakuraweb.com/kougi/keisan_tenmongakuII/index.html
/////////////////////

#include <iosfwd>
#include <iostream>

typedef std::istream istream;
typedef std::ostream ostream;

class Vector{
	public:
		double element[VLEN];

		//constructor
		Vector(double data=0){
			for(int i=0;i<VLEN;i++){
				this->element[i]=data;
			}
		}
		Vector(double data[VLEN]){
			for(int i=0;i<VLEN;i++){
				this->element[i]=data[i];
			}
		}

		template<typename... Data>
		Vector(Data... data){
			int i=0;
			for(double d: std::initializer_list<double>{data...}){
				this->element[i]=d;
				i++;
			}
		}
		double& operator [](int i){
			return element[i];
		}

		friend const Vector operator + (const Vector& v1,const Vector& v2);
		friend const Vector operator - (const Vector& v1,const Vector& v2);
		friend const double operator * (const Vector& v1,const Vector& v2);
		friend const Vector operator * (const Vector& v,const double& d);
		friend const Vector operator * (const Vector& v,const int& d);
		friend const Vector operator / (const Vector& v,const double& d);
		friend const Vector operator / (const Vector& v,const int& d);

		Vector& operator += (const Vector& v){
			for(int i=0;i<VLEN;i++){
				this->element[i]+=v.element[i];
			}
			return *this;
		}

		Vector& operator -= (const Vector& v){
			for(int i=0;i<VLEN;i++){
				this->element[i]-=v.element[i];
			}
			return *this;
		}

		Vector& operator *= (const double d){
			for(int i=0;i<VLEN;i++){
				this->element[i]*=d;
			}
			return *this;
		}
		Vector& operator *= (const int d){
			for(int i=0;i<VLEN;i++){
				this->element[i]*=d;
			}
			return *this;
		}
		Vector& operator /= (const double d){
			for(int i=0;i<VLEN;i++){
				this->element[i]/=d;
			}
			return *this;
		}
		Vector& operator /= (const int d){
			for(int i=0;i<VLEN;i++){
				this->element[i]/=d;
			}
			return *this;
		}

		friend ostream& operator << (ostream& s,const Vector& v);
		friend istream& operator >> (istream& s,Vector& v);


//		friend std::ostream& operator << (std::ostream& s,const Vector& v);
//		friend std::istream& operator >> (std::istream& s,Vector& v);
	//public end
};

inline ostream& operator << (ostream& s,const Vector& v){
	s<<"[";
	for(int i=0;i<VLEN;i++){
		if(i!=0){
			s<<",";
		}
		s<<v.element[i];
	}
	s<<"]"<<std::endl;
	return s;
}

inline istream& operator >> (istream& s,Vector& v){
	for(int i=0;i<VLEN;i++){
		s>>v.element[i];
	}
	return s;
}


inline const Vector operator + (const Vector& v1,const Vector& v2){
	Vector ret;
	for(int i=0;i<VLEN;i++){
		ret.element[i]=v1.element[i]+v2.element[i];
	}
	return ret;
}

inline const Vector operator - (const Vector& v1,const Vector& v2){
	Vector ret;
	for(int i=0;i<VLEN;i++){
		ret.element[i]=v1.element[i]-v2.element[i];
	}
	return ret;
}

inline const double operator * (const Vector& v1,const Vector& v2){
	double ret=0;
	for(int i=0;i<VLEN;i++){
		ret+=v1.element[i]*v2.element[2];
	}
	return ret;
}

inline const Vector operator * (const Vector& v,const double d){
	Vector ret;
	for(int i=0;i<VLEN;i++){
		ret.element[i]=v.element[i]*d;
	}
	return ret;
}

inline const Vector operator * (const Vector& v,const int d){
	Vector ret;
	for(int i=0;i<VLEN;i++){
		ret.element[i]=v.element[i]*d;
	}
	return ret;
}

inline const Vector operator / (const Vector& v,const double d){
	Vector ret;
	for(int i=0;i<VLEN;i++){
		ret.element[i]=v.element[i]/d;
	}
	return ret;
}

inline const Vector operator / (const Vector& v,const int d){
	Vector ret;
	for(int i=0;i<VLEN;i++){
		ret.element[i]=v.element[i]/d;
	}
	return ret;
}

template<class head,class... tails>
Vector grvPtTop(head vec1,tails... vecs){
	int cnt=1;
	Vector ret=grvPt(cnt,vec1,vecs...);
	ret/=cnt;
	return ret;
}

template<class head,class... tails>
Vector grvPt(int& cnt,head vec1,tails... vecs){
	Vector ret=grvPt(cnt,vecs...);
	ret+=vec1;
	cnt++;
	return ret;
}

Vector grvPt(int& cnt,Vector v1,Vector v2){
	Vector ret=v1+v2;
	cnt++;
	return ret;
}

typedef const Vector (vfunc)(const Vector&);

#endif	//ifndef VECTOR_H

