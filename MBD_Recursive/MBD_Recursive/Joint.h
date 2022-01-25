#pragma once
#include "Body.h"
using namespace Eigen;
class JointBase
{
	typedef Matrix<double, 3, -1> MatR3CX;
	typedef Matrix<double, 3, 4> MatR3C4;
protected:
	Body* Bi;                              //pointer to outer body.
	Body* Bj;                              //pointer to inner body.
	double* yi;                            //generalized coordinates.
	double* dyi;                           //generalized velocity.
	int NP;                                //P elements ID, -1 represents Bi is not flexible body.
	int NQ;                                //Q elements ID, -1 represents Bj is not flexible body.
	Vector3d rhoi;                         //rhop vector
	Vector3d rhoj;                         //rhoq vector
	Matrix3d BiP;                          //transform from P0 to P
	Matrix3d CiP;                          //transform from P to h.   Constant
	Matrix3d BjQ;                          //transform from Q0 to Q.
	Matrix3d CjQ;                          //transform from Q to h0.  Constant
	Matrix3d Dih;                          //transform from h0 to h.
	MatrixXd Ui;                           //need derived class to resize it.
	VectorXd betai;          
	MatrixXd Tij;
protected:
	Matrix3d& calCiP();
	Matrix3d& calCjQ();
	virtual Matrix3d& calBiP();
	virtual Matrix3d& calBjQ();
	virtual bool calAh0(OUT Matrix3d& A);
	virtual Vector3d Hi()const;
	virtual Vector3d Wri()const;
	virtual Vector3d Vri()const;
private:
	/*need derived class to complete it.*/
	virtual MatR3CX HhT()const = 0;
	virtual MatR3CX HoT()const = 0;
	virtual Vector3d Yita()const = 0;                                             
	virtual Matrix3d& calDih() = 0;
public:
	enum
	{
		REVOLUTIONAL = 0x0,
		UNIVERSAL = 0x1,
		SPHERICAL=0x2,
		PRISMATIC=0x3,
		CYLINDRICAL=0x4,
		VIRTUAL=0x5,
	};
	friend class TreeSystem;
	friend class MBFileParser;
	JointBase(Body* Bi_ptr, Body* Bj_ptr, const Vector3d& rho_i, const Vector3d& rho_j);
	virtual ~JointBase();
	/*need derived class to complete it.*/
	inline virtual unsigned int DOF()const = 0;
	inline virtual unsigned int type()const = 0;
	/*some common function for the Joint class */
	virtual MatrixXd& calUi(IN double t, IN double* y);
	virtual VectorXd& calbetai(IN double t, IN double* y, IN double* dy);
	virtual MatrixXd& calTij(IN double t,IN double* y);
	virtual bool calytoq(IN double t,IN double* y);                             //calculate from y to q
	virtual bool caldytodq(IN double t, IN double* dy);                         //calculate from dy to dq
	virtual bool setCiP(const Matrix3d& c);                 //CiP CjQ is transformation matrix from
	virtual bool setCjQ(const Matrix3d& c);                 //element cords to joint cords
	/*function to write some joint config*/
	virtual bool Write(Json::Value& joint)const;
};

/*the rotation axis is along x-axis by default.*/
class Revolute :public JointBase
{
	typedef Matrix<double, 3, -1> MatR3CX;
private:
	virtual MatR3CX HhT()const;
	virtual MatR3CX HoT()const;
	virtual Vector3d Yita()const;
	virtual Matrix3d& calDih();
public:
	Revolute(Body* Bi_ptr, Body* Bj_ptr,Vector3d& rho_i,Vector3d& rho_j);
	~Revolute();
	inline unsigned int DOF()const;
	inline unsigned int type()const;
};

class Universe :public JointBase
{
	typedef Matrix<double, 3, -1> MatR3CX;
private:
	virtual MatR3CX HhT()const;
	virtual MatR3CX HoT()const;
	virtual Vector3d Yita()const;
	virtual Matrix3d& calDih();
public:
	Universe(Body* Bi_ptr, Body* Bj_ptr, Vector3d& rho_i, Vector3d& rho_j);
	~Universe();
	inline unsigned int DOF()const;
	inline unsigned int type()const;
};

class Virtual :public JointBase
{
	typedef Matrix<double, 3, -1> MatR3CX;
private:
	virtual MatR3CX HhT()const;
	virtual MatR3CX HoT()const;
	virtual Vector3d Yita()const;
	virtual Matrix3d& calDih();
public:
	Virtual(Body* Bi_ptr, Body* Bj_ptr, Vector3d& rho_i, Vector3d& rho_j);
	~Virtual();
	inline unsigned int DOF()const;
	inline unsigned int type()const;
};

class Sphere :public JointBase
{
	typedef Matrix<double, 3, -1> MatR3CX;
private:
	virtual MatR3CX HhT()const;
	virtual MatR3CX HoT()const;
	virtual Vector3d Yita()const;
	virtual Matrix3d& calDih();
public:
	Sphere(Body* Bi_ptr, Body* Bj_ptr, Vector3d& rho_i, Vector3d& rho_j);
	~Sphere();
	inline unsigned int DOF()const;
	inline unsigned int type()const;
};

class Prism :public JointBase
{
	typedef Matrix<double, 3, -1> MatR3CX;
private:
	virtual MatR3CX HhT()const;
	virtual MatR3CX HoT()const;
	virtual Vector3d Yita()const;
	virtual Matrix3d& calDih();
public:
	Prism(Body* Bi_ptr, Body* Bj_ptr, Vector3d& rho_i, Vector3d& rho_j);
	~Prism();
	inline unsigned int DOF()const;
	inline unsigned int type()const;
};

class Cylinder :public JointBase
{
	typedef Matrix<double, 3, -1> MatR3CX;
private:
	virtual MatR3CX HhT()const;
	virtual MatR3CX HoT()const;
	virtual Vector3d Yita()const;
	virtual Matrix3d& calDih();
public:
	Cylinder(Body* Bi_ptr, Body* Bj_ptr, Vector3d& rho_i, Vector3d& rho_j);
	~Cylinder();
	inline unsigned int DOF()const;
	inline unsigned int type()const;
};