#pragma once
#include "Body.h"
#include "Joint.h"
#include <vector>
#include <list>
#include "tree.h"
using namespace std;
using namespace Eigen;


class TreeSystem
{
public:
	struct BJpair
	{
		Body* bptr;
		JointBase* jptr;
		BJpair(Body* b = nullptr, JointBase* j = nullptr) { bptr = b; jptr = j; }
		~BJpair() {};
		bool operator==(BJpair& other);
		bool operator<(BJpair& other);
		bool operator>(BJpair& other);
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
	bool calytoq(double t,double* y);                   //calculate from y to q.
	bool caldytodq(double t,double* dy);                //calculate from dy to dq.
	bool update(double t, double* y, double* dy);               //update q dq Tij Ui betai
public:
	friend class MBSystem;
	TreeSystem(vector<JointBase*>& jvec);
	bool calG0(double t, double* y, OUT MatrixXd& G0);
	bool calG(double t, double* y, OUT MatrixXd& G);
	bool calM(double t, double* y, OUT MatrixXd& M);
	bool calg(double t, double* y, double* dy, OUT MatrixXd& g);
	bool calf(double t, double* y, double* dy, OUT VectorXd& f);
	
	VectorXd rootacc(double t);
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
	unsigned int nc = 0;                       //number of independent generalized coordinates.
	TreeSystem* mbtree = nullptr;
	Equation* peq = nullptr;
	Solver* psolver = nullptr;
protected:
	//vector<Body*> bodyvec;
	vector<JointBase*> jointvec;
	VectorXd y0;
	MatrixXd G;
	MatrixXd G0;
	MatrixXd M;
	MatrixXd g;
	MatrixXd Z;
	VectorXd z;
	VectorXd fey;
	VectorXd f;
private:
	MatrixXd& calG0(double t, double* y);
	MatrixXd& calG(double t,double* y);
	MatrixXd& calg(double t,double* y,double* dy);
	MatrixXd& calM(double t,double* y);
	VectorXd& calf(double t, double* y, double* dy);
	bool update(double t, double* y,double* dy);       //update G g M f G0
	bool initialize();                                //call before calculation
protected:
	MatrixXd& calZ(double t, double* y);
	VectorXd& calz(double t, double* y, double* dy);
	VectorXd& calfey(double t, double* y, double* dy);       //remain completed.
public:
	friend class Equation;
	MBSystem();
	~MBSystem();
	bool add(JointBase* j);
	bool del(JointBase* j);
	unsigned int DOF()const;                          //the number of the generalized coordinates of the system
	bool sety0(const VectorXd& y0);
	bool setTimeInterval(double ti, double te, int N);
	bool setTolerance(double r = 1e-4, double a = 1e-3);
	bool calculate();
	bool SaveAs(string fname, bool isbinary = false);
};


class Equation
{
private:
	MBSystem* pmbs;
protected:
	MatrixXd L;
	VectorXd R;
public:
	Equation(MBSystem* p);
	~Equation();
	MatrixXd& Left(double t, VectorXd& y);                    //y是状态空间变量
	VectorXd& Right(double t, VectorXd& y);
	VectorXd initialvalue()const;
	unsigned int DOF()const;
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
private:
	vector<double> tspan;
	list<VectorXd> Y;
	list<VectorXd> DY;
public:
	friend MBSystem;
	Solver(Equation* p);
	~Solver();
	bool setTimeInterval(double ti, double te, int N);
	bool setTolerance(double r=1e-4, double a=1e-3);
	bool calculate();
};