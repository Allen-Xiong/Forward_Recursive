#pragma once

#include "auxiliary.h"
#include <vector>
using namespace std;
using namespace Eigen;
class Body
{
protected:
	int id = -1;                               // id 
	double tM = 0.;                            //total Mass
	double* pos = nullptr;
	double* vel = nullptr;
	double* acc = nullptr;
public:
	enum
	{
		RIGID = 0x00,
		FLEXIBLE = 0x01,
		BASE = 0x02
	};
	static RCORDS m_s_rtype;                   //type of rotation coordinates.
	static unsigned int NC;                    //number of the absolute coordinates.
	static bool setRcords(IN RCORDS rtype);
	Body(IN double m)noexcept;
	Body()noexcept;
	bool operator==(IN Body& other)const;
	bool operator<(IN Body& other)const;
	bool operator<=(IN Body& other)const;
	bool operator>(IN Body& other)const;
	bool A(OUT Matrix3d& M)const;
	virtual bool setID(IN int i);
	virtual unsigned int type()const = 0;
	virtual bool calMass(OUT MatrixXd& M,IN int k) = 0;                   // M is the system mass matrix, k represents the begin location.
	virtual unsigned int nMode()const = 0;
	virtual Vector3d angularVel()const;
	virtual Vector3d angularVel(IN int EID)const;                      // angular velocity of the element on the body.
	virtual Vector3d rhoP(IN Vector3d& _rp,IN int EID)const;             //the conjoined vector of element P on the body
	virtual Vector3d translationalVel(IN Vector3d& _rp)const;          //the translational velocity of P on the body.
	virtual Vector3d translationalVel(IN Vector3d& _rp,IN int EID)const; //the translational velocity of element P on the body.
	virtual VectorXd inertiaForce();                                //the inertia force of the body.
	virtual ~Body();
	virtual bool Write(OUT Json::Value& body)const;
protected:
	virtual Vector3d uiP(IN int EID)const;                             //the relative velocity of the element P about floating coordinate system.
public:
	friend class TreeSystem;
	friend class JointBase;
};

class BaseBody :public Body
{
	VectorXd(*pFun)(double) = nullptr;
	VectorXd(*vFun)(double) = nullptr;
	VectorXd(*aFun)(double) = nullptr;
public:
	BaseBody();
	BaseBody(IN VectorXd(*p)(double),IN VectorXd(*v)(double),IN VectorXd(*a)(double));
	~BaseBody();
	virtual bool Write(OUT Json::Value& body)const;
	unsigned int type()const;
	bool calMass(OUT MatrixXd& M,IN int k);
	unsigned int nMode()const;
	VectorXd acceleration(IN double t);                     /*remain completed.*/
	bool update(IN double t);
};

class RigidBody :public Body
{
protected:
	Matrix3d M11;
	Matrix3d M12;
	Matrix3d M22;
	Matrix3d Jc;                                        //the inertia of moment of the centroids
	VectorXd InerF;                                       //the inertia force.
protected:
	virtual Matrix3d& calM11();
	virtual Matrix3d& calM12();
	virtual Matrix3d& calM22();
public:
	RigidBody();
	RigidBody(IN double m);
	RigidBody(IN double m,IN Matrix3d& I);
	virtual unsigned int type()const;
	virtual bool calMass(OUT MatrixXd& M,IN int k);            
	virtual unsigned int nMode()const;
	virtual VectorXd inertiaForce();
	virtual ~RigidBody();
	virtual bool Write(OUT Json::Value& body)const;
};

class FlexibleBody :public RigidBody
{
	typedef Matrix<double, 3, -1> MatR3CX;
private:
	/*
	* g5,G4,G5 is related to the mode coordinates.
	* so every time they should be recalculated.
	* the rest g1,g2,g3,g4,G1,G2,G3 are independent
	* to the mode coordinates.
	*/
	Vector3d g1;
	MatR3CX g2;
	MatR3CX g3;
	vector<vector<Vector3d>> g4;
	MatR3CX g5;
	Matrix3d G1;
	vector<Matrix3d> G2;
	vector<vector<Matrix3d>> G3;
	vector<Matrix3d> G4;
	vector<Matrix3d> G5;
protected:
	unsigned int s;          //number of mode
	unsigned int ne;         //number of elements
	VectorXd me;             //ne*1
	MatR3CX rho;            //3*ne
	vector<MatR3CX> PHI;    //3*s*ne
	vector<MatR3CX> PSI;    //3*s*ne
	MatrixXd Ka;             //s*s
	MatrixXd Ca;             //s*s
	MatR3CX M13;             //3*s
	MatR3CX M23;             //3*s
	MatrixXd M33;            //s*s
protected:
	virtual Matrix3d& calM11();
	virtual Matrix3d& calM12();
	virtual Matrix3d& calM22();
	virtual MatR3CX& calM13();
	virtual MatR3CX& calM23();
	virtual MatrixXd& calM33();
	virtual Vector3d uiP(IN int EID)const;
private:
	bool calg5();
	bool calG4();
	bool calG5();
public:
	friend class JointBase;
	friend class ActiveForce;
	FlexibleBody(IN int NE,IN int nmode);
	bool setMe(IN VectorXd& Me);
	bool setRho(IN MatR3CX& Rho);
	bool setPHI(IN vector<MatrixXd>& phi);
	bool setPSI(IN vector<MatrixXd>& psi);
	bool setKa(IN MatrixXd& ka);
	bool setCa(IN MatrixXd& ca);
	virtual unsigned int type()const;
	virtual bool calMass(OUT MatrixXd& M,IN int k);
	virtual unsigned int nMode()const;
	virtual Vector3d angularVel(IN int EID)const;                      // angular velocity of the element on the body.
	virtual Vector3d rhoP(IN Vector3d& _rp,IN int EID)const;             //the conjoined vector of element P on the body
	virtual Vector3d translationalVel(IN Vector3d& _rp,IN int EID)const;
	virtual VectorXd inertiaForce();
	virtual ~FlexibleBody();
	virtual bool Write(OUT Json::Value& body)const;
};
