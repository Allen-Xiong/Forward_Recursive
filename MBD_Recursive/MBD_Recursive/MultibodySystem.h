#pragma once
#include <vector>
#include <list>
#include "tree.h"
#include "Body.h"
#include "Joint.h"
#include "auxiliary.h"
using namespace std;
using namespace Eigen;


class TreeSystem
{
public:
	struct BJpair
	{
		Body* bptr;
		JointBase* jptr;
		BJpair(IN Body* b = nullptr,IN JointBase* j = nullptr) { bptr = b; jptr = j; }
		~BJpair() {};
		bool operator==(IN BJpair& other);
		bool operator<(IN BJpair& other);
		bool operator>(IN BJpair& other);
	};
private:
	myTree<BJpair> MBTree;
	vector<BJpair*> bjvec;                     //the bodies in the system, sort by their id. By default, the Joint id is the same as the outer body
	unsigned int nb;                           //the number of the bodies in the system
	BJpair* bjroot;                            //root body, it doesn't exit in the system in reality.
	VectorXi qdimvec;                          //record the number of the absolute coordinates of every Body.
	VectorXi ydimvec;                          //record the number of the generalized coordinates of every Body.
	VectorXi qdimcum;                          //record the beginning position in system coordinates of every Body.
	VectorXi ydimcum;                          //record the beginning position in system coordinates of every Body.
	VectorXd yprev;                            //generalized coordinates
	VectorXd dyprev;                           //generalized velocity
	VectorXd qprev;                            //absolute coordinates
	VectorXd dqprev;                           //absolute velocity
private:
	bool calytoq(IN double t,IN double* y);                   //calculate from y to q.
	bool caldytodq(IN double t,IN double* dy);                //calculate from dy to dq.
	bool update(IN double t, IN double* y,IN double* dy);               //update q dq Tij Ui betai
public:
	friend class MBSystem;
	friend class Solver;
	TreeSystem(IN vector<JointBase*>& jvec);
	bool calG0(IN double t,IN double* y, OUT MatrixXd& G0);
	bool calG(IN double t,IN double* y, OUT MatrixXd& G);
	bool calM(IN double t,IN double* y, OUT MatrixXd& M);
	bool calg(IN double t,IN double* y,IN double* dy, OUT MatrixXd& g);
	bool calf(IN double t,IN double* y,IN double* dy, OUT VectorXd& f);
	
	VectorXd rootacc(IN double t);
	unsigned int DOF()const;
	unsigned int NC()const;
	~TreeSystem();

};

class Solver;
class Equation;
class MBSystem
{
private:
	unsigned int dof = 0;                      //degree of freedom
	unsigned int nc = 0;                       //number of absolute coordinates.
	TreeSystem* mbtree = nullptr;
	Equation* peq = nullptr;
	Solver* psolver = nullptr;
protected:
	vector<JointBase*> jointvec;
	vector<int> driveid;
	VectorXd y0;
	VectorXd dy0;
	MatrixXd G;
	MatrixXd G0;
	MatrixXd M;
	MatrixXd g;
	MatrixXd Z;
	VectorXd z;
	VectorXd fey;
	VectorXd f;
	MatrixXd Zbar;
	VectorXd zbar;
private:
	MatrixXd& calG0(IN double t,IN double* y);
	MatrixXd& calG(IN double t,IN double* y);
	MatrixXd& calg(IN double t,IN double* y,IN double* dy);
	MatrixXd& calM(IN double t,IN double* y);
	VectorXd& calf(IN double t,IN double* y,IN double* dy);
	bool update(IN double t,IN double* y,IN double* dy);       //update G g M f G0
	bool initialize();                                         //call before calculation
protected:
	MatrixXd& calZ(IN double t,IN double* y);
	VectorXd& calz(IN double t,IN double* y,IN double* dy);
	VectorXd& calfey(IN double t,IN double* y,IN double* dy);       //remain completed.
	MatrixXd& calZbar(IN double t,IN double* y);
	VectorXd& calzbar(IN double t,IN double* y,IN double* dy);
public:
	friend class Equation;
	friend class Solver;
	friend class MBFileParser;
	MBSystem();
	~MBSystem();
	bool add(IN JointBase* j);
	bool del(IN JointBase* j);
	unsigned int DOF()const;                          //the number of the generalized coordinates of the system
	bool sety0(IN const VectorXd& y0);
	bool setdy0(IN const VectorXd& _dy0);
	bool setTimeInterval(IN double ti,IN double te, IN int N);
	bool setTolerance(IN double r = 1e-4, IN double a = 1e-3);
	bool setToInitial();                              //set postion and velocity to initial state.
	bool calculate();
	bool SaveAs(IN string fname, IN bool isbinary = false);
};


class Equation
{
private:
	MBSystem* pmbs;
protected:
	MatrixXd L;
	VectorXd R;
public:
	friend class Solver;
	Equation(IN MBSystem* p);
	~Equation();
	MatrixXd& Left(IN double t,IN VectorXd& y);                    //y is state space variable.
	VectorXd& Right(IN double t,IN VectorXd& y);
	VectorXd initialvalue()const;
	unsigned int DOF()const;
	void Initialize();
};


class Solver
{
private:
	Equation* pe;
	double t_ini;
	double t_end;
	int Nstep;
	double Rtol;
	double Atol;
	vector<double> tspan;
	list<VectorXd> Y;
	list<VectorXd> DY;
	list<VectorXd> Q;
	list<VectorXd> DQ;
public:
	friend class MBSystem;
	friend class MBFileParser;
	Solver(IN Equation* p);
	~Solver();
	bool setTimeInterval(IN double ti,IN double te,IN int N);
	bool setTolerance(IN double r=1e-4,IN double a=1e-3);
	bool calculate();
};

class MBFileParser
{
	MBSystem* pmbs = nullptr;
	vector<Body*> bodyvec;
	bool freememo = false;
protected:
	void clear();
	void CheckId(IN int id);
	void CheckMass(IN double m);
	void CheckJc(IN const Json::Value& Jc);
	void CheckRho(IN const Json::Value& rho);
	void CheckMat3d(IN const Json::Value& val);
	void CheckPos(IN const Json::Value& val, IN unsigned int k);
	void CheckVel(IN const Json::Value& val, IN unsigned int k);
	void GetJc(IN const Json::Value& Jc, OUT Matrix3d& Ic);
	void GetRho(IN const Json::Value& rho, OUT Vector3d& r);
	void GetMat3d(IN const Json::Value& val,OUT Matrix3d& M);
	void GetPos(IN const Json::Value& val,OUT VectorXd& p);
	void GetVel(IN const Json::Value& val,OUT VectorXd& v);
	void GetFlexibleBody(IN const string& fname,OUT Body* &p);
public:
	MBFileParser() :pmbs(nullptr) {};
	~MBFileParser();
	bool Read(IN const string& fname);  
	bool Write(IN const string& fname);  
	bool Simulate();
	bool SaveDataAs(IN const string& fname, IN bool isbinary);
};