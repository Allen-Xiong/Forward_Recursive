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
		REVOLUTIONAL = 0x00,
		UNIVERSAL = 0x01,
		SPHERICAL=0x02,
		PRISMATIC=0x03,
		CYLINDRICAL=0x04,
		VIRTUAL=0x05,

		REVOLUTIONALDRIVE = 0x10,
		UNIVERSALDRIVE = 0x11,
		SPHERICALDRIVE = 0x12,
		PRISMATICDRIVE = 0x13,
		CYLINDRICALDRIVE = 0x14,
		VIRTUALDRIVE = 0x15,
	};
	friend class TreeSystem;
	friend class MBFileParser;
	JointBase(IN Body* Bi_ptr,IN Body* Bj_ptr,IN const Vector3d& rho_i,IN const Vector3d& rho_j);
	virtual ~JointBase();
	virtual VectorXd Acceleration(IN double t);
	virtual VectorXd Velocity(IN double t);
	virtual VectorXd Position(IN double t);
	bool operator<(IN const JointBase& other)const;
	/*served for drive joints*/
	virtual void Acceleration(IN double t,OUT double* ddy)const;
	virtual void Velocity(IN double t,OUT double* dy)const;
	virtual void Position(IN double t,OUT double* y)const;
	/*need derived class to complete it.*/
	inline virtual unsigned int DOF()const = 0;
	inline virtual unsigned int type()const = 0;
	/*some common function for the Joint class */
	virtual MatrixXd& calUi(IN double t, IN double* y);
	virtual VectorXd& calbetai(IN double t, IN double* y, IN double* dy);
	virtual MatrixXd& calTij(IN double t,IN double* y);
	virtual bool calytoq(IN double t,IN double* y);                             //calculate from y to q
	virtual bool caldytodq(IN double t, IN double* dy);                         //calculate from dy to dq
	virtual bool setCiP(IN const Matrix3d& c);                                  //CiP CjQ is transformation matrix from
	virtual bool setCjQ(IN const Matrix3d& c);                                  //element cords to joint cords
	/*function to write some joint config*/
	virtual bool Write(OUT Json::Value& joint)const;
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
	Revolute(IN Body* Bi_ptr,IN Body* Bj_ptr,IN Vector3d& rho_i,IN Vector3d& rho_j);
	~Revolute();
	inline unsigned int DOF()const;
	virtual inline unsigned int type()const;
};

/*the rotation axis is along x-axis and y-axis by default*/
class Universe :public JointBase
{
	typedef Matrix<double, 3, -1> MatR3CX;
private:
	virtual MatR3CX HhT()const;
	virtual MatR3CX HoT()const;
	virtual Vector3d Yita()const;
	virtual Matrix3d& calDih();
public:
	Universe(IN Body* Bi_ptr,IN Body* Bj_ptr,IN Vector3d& rho_i,IN Vector3d& rho_j);
	~Universe();
	inline unsigned int DOF()const;
	virtual inline unsigned int type()const;
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
	Virtual(IN Body* Bi_ptr,IN Body* Bj_ptr,IN Vector3d& rho_i,IN Vector3d& rho_j);
	~Virtual();
	inline unsigned int DOF()const;
	virtual inline unsigned int type()const;
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
	Sphere(IN Body* Bi_ptr,IN Body* Bj_ptr,IN Vector3d& rho_i,IN Vector3d& rho_j);
	~Sphere();
	inline unsigned int DOF()const;
	virtual inline unsigned int type()const;
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
	Prism(IN Body* Bi_ptr,IN Body* Bj_ptr,IN Vector3d& rho_i,IN Vector3d& rho_j);
	~Prism();
	inline unsigned int DOF()const;
	virtual inline unsigned int type()const;
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
	Cylinder(IN Body* Bi_ptr,IN Body* Bj_ptr,IN Vector3d& rho_i,IN Vector3d& rho_j);
	~Cylinder();
	inline unsigned int DOF()const;
	virtual inline unsigned int type()const;
};



class RevoluteDrive :public Revolute
{
private:
	double (*pf1)(double);
public:
	RevoluteDrive(IN Body* Bi_ptr,IN Body* Bj_ptr,IN Vector3d& rho_i,IN Vector3d& rho_j,IN double (*f)(double));
	~RevoluteDrive();
	virtual inline unsigned int type()const;
	virtual bool calytoq(IN double t, IN double* y);                            
	virtual bool caldytodq(IN double t, IN double* dy);
	virtual VectorXd Acceleration(IN double t);
	virtual VectorXd Velocity(IN double t);
	virtual VectorXd Position(IN double t);
	virtual void Acceleration(IN double t,OUT double* ddy)const;
	virtual void Velocity(IN double t,OUT double* dy)const;
	virtual void Position(IN double t,OUT double* y)const;
};

class UniverseDrive :public Universe
{
private:
	double (*pf1)(double);
	double (*pf2)(double);
public:
	UniverseDrive(IN Body* Bi_ptr,IN Body* Bj_ptr,IN Vector3d& rho_i,IN Vector3d& rho_j,IN double (*q1)(double),IN double (*q2)(double));
	~UniverseDrive();
	virtual bool calytoq(IN double t, IN double* y);
	virtual bool caldytodq(IN double t, IN double* dy);
	virtual VectorXd Acceleration(IN double t);
	virtual VectorXd Velocity(IN double t);
	virtual VectorXd Position(IN double t);
	virtual void Acceleration(IN double t,OUT double* ddy)const;
	virtual void Velocity(IN double t,OUT double* dy)const;
	virtual void Position(IN double t,OUT double* y)const;
};