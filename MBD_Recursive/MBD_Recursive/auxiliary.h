#pragma once

#include <Dense>
#include <Eigen>
#include <iostream>
#include <fstream>
#include <string>
#include "json.h"
#ifndef IN
#define IN
#endif // !IN

#ifndef OUT
#define OUT
#endif // !OUT


constexpr double PI=3.1415926535898;

#define INI_FAILURE(str) string(str)+string(" Initialization Failure")


using namespace Eigen;
using namespace std;
typedef Matrix<double, 3, -1> MatR3CX;

enum class RCORDS
{
	EULERANGLE=0x01,
	EULERQUATERNION=0x02,
	CARDANANGLE=0x03,
};
namespace AUX
{
	extern double Tol;
	/*some auxiliary functions*/
	bool A(IN double* q,OUT Matrix3d& M);    //EQtoA
	bool EQtoA(IN double* q, OUT Matrix3d& M);
	bool EAtoA(IN double* q, OUT Matrix3d& M);
	bool CAtoA(IN double* q, OUT Matrix3d& M);
	bool R(IN double* q, OUT Matrix<double, 3, 4>& M);
	bool tilde(IN Vector3d& v, OUT Matrix3d& M);
	bool tilde(IN double* q, OUT Matrix3d& M);
	bool isEqual(const double* v, const double* w, const unsigned int& len);
	/*transform A to Euler Quaternion*/
	bool AtoEQ(IN const Matrix3d& A, OUT double* q);
	/*transform A to Euler Angle*/
	bool AtoEA(IN const Matrix3d& A, OUT double* q);
	/*transform A to Cardan Angle*/
	bool AtoCA(IN const Matrix3d& A, OUT double* q);
	/*calculate Kr matrix: transformation matrix from Rotation cords to angular velocity*/
	bool Kr(IN double* q,IN RCORDS rtype, OUT MatR3CX& K);

}

class MBException
{
private:
	string msg;
public:
	MBException(const string& s) { msg = s; };
	~MBException() = default;
	virtual const string& what() const { return msg; };
};

