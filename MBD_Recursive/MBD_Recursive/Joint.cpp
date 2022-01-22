#include "Joint.h"

/*check once*/
MatR3CX Revolute::HhT() const
{
	return Vector3d(0, 0, 0);
}
/*check once*/
MatR3CX Revolute::HoT() const
{
	return Vector3d(1, 0, 0);
}
/*check once*/
Vector3d Revolute::Yita() const
{
	return Vector3d(0, 0, 0);
}
/*check once*/
Matrix3d& Revolute::calDih()
{
	// Assume rotate along x-axis.
	Dih.setIdentity();
	Dih(1, 1) = cos(yi[0]);
	Dih(1, 2) = -sin(yi[0]);
	Dih(2, 1) = -Dih(1, 2);
	Dih(2, 2) = Dih(1, 1);
	return Dih;
}

Revolute::Revolute(Body* Bi_ptr, Body* Bj_ptr, Vector3d& rho_i, Vector3d& rho_j):JointBase(Bi_ptr,Bj_ptr,rho_i,rho_j)
{
	unsigned int si = Bi_ptr->nMode();
	Ui.resize(6 + si, DOF() + si);
	Ui.setZero();
	if (si != 0)
	{
		Ui.block(6, DOF(), si, si).setIdentity();
	}
}

Revolute::~Revolute()
{
}

inline unsigned int Revolute::DOF() const
{
	return 1;
}

inline unsigned int Revolute::type() const
{
	return REVOLUTIONAL;
}
/*check once*/
Matrix3d& JointBase::calCiP()
{
	return CiP;
}
/*check once*/
Matrix3d& JointBase::calCjQ()
{
	return CjQ;
}
/*check once*/
Matrix3d& JointBase::calBiP()
{
	// TODO: insert return statement here
	if (Bi->type() != Body::FLEXIBLE)
		return BiP;
	/* important!!! ai should get from yi rather than qi*/
	Map<VectorXd> a(yi + DOF(), Bi->nMode());
	auto& PSI_iP = static_cast<FlexibleBody*>(Bi)->PSI[NP];
	Vector3d v = PSI_iP * a;
	Matrix3d til;
	AUX::tilde(v, til);
	BiP.setIdentity();
	BiP += til;
	return BiP;
}
/*check once*/
Matrix3d& JointBase::calBjQ()
{
	// TODO: insert return statement here
	if (Bj->type() != Body::FLEXIBLE)
		return BjQ;
	Map<VectorXd> a(Bj->pos + Body::NC, Bj->nMode());
	auto& PSI_jQ = static_cast<FlexibleBody*>(Bj)->PSI[NQ];
	Vector3d v = PSI_jQ * a;
	Matrix3d til;
	AUX::tilde(v, til);
	BjQ.setIdentity();
	BjQ += til;
	return BjQ;
}
/*check once*/
bool JointBase::calAh0(OUT Matrix3d& A)
{
	Bj->A(A);
	A *= calBjQ() * calCjQ();
	return true;
}
/*check once*/
Vector3d JointBase::Hi() const
{
	Map<VectorXd> y(yi, DOF());
	return HhT() * y;
}
/*check once*/
Vector3d JointBase::Wri() const
{
	Map<VectorXd> dy(dyi, DOF());
	return HoT() * dy;
}
/*check once*/
Vector3d JointBase::Vri() const
{
	Map<VectorXd> dy(dyi, DOF());
	return HhT() * dy;
}


/*check once*/
JointBase::JointBase(Body* Bi_ptr, Body* Bj_ptr, const Vector3d& rho_i, const Vector3d& rho_j)
{
	Bi = Bi_ptr;
	Bj = Bj_ptr;
	rhoi = rho_i;
	rhoj = rho_j;
	if (Bi->type() != Body::FLEXIBLE)
		NP = -1;
	else
	{
		FlexibleBody* pb = static_cast<FlexibleBody*>(Bi);
		(pb->rho - rhoi).colwise().squaredNorm().minCoeff(&NP);
	}
	if (Bj->type() != Body::FLEXIBLE)
		NQ = -1;
	else
	{
		FlexibleBody* pb = static_cast<FlexibleBody*>(Bj);
		(pb->rho - rhoj).colwise().squaredNorm().minCoeff(&NQ);
	}
	BiP.setIdentity();
	CiP.setIdentity();
	BjQ.setIdentity();
	CjQ.setIdentity();
	Dih.setIdentity();
	unsigned int si = Bi->nMode();
	unsigned int sj = Bj->nMode();
	betai.resize(6 + si);
	Tij.resize(6 + si, 6 + sj);
	betai.setZero();
	Tij.setZero();
	Tij.block<3, 3>(0, 0).setIdentity();
	Tij.block<3, 3>(3, 3).setIdentity();
}

JointBase::~JointBase()
{
}

inline unsigned int JointBase::DOF() const
{
	return 0;
}
/*check once*/
MatrixXd& JointBase::calUi(IN double t, IN double* y)
{
	Vector3d rhoip = rhoi;
	Vector3d ho, hh;
	Matrix3d Ai, Ah0;
	Bi->A(Ai);
	calAh0(Ah0);
	ho = Ah0 * HoT();
	hh = Ah0 * HhT();
	Matrix3d r;
	int deli = DOF();
	if (Bi->type() != Body::FLEXIBLE)
	{
		rhoip = Ai * rhoip;
		AUX::tilde(rhoip.data(), r);
		Ui.block(0, 0, 3, deli) = hh + r * ho;
		Ui.block(3, 0, 3, deli) = ho;
	}
	else
	{
		unsigned int si = Bi->nMode();
		Map<VectorXd> ai(yi + DOF(), si);
		auto& PHI_iP = static_cast<FlexibleBody*>(Bi)->PHI[NP];
		auto& PSI_iP = static_cast<FlexibleBody*>(Bi)->PSI[NP];
		rhoip += PHI_iP * ai;
		rhoip = Ai * rhoip;
		AUX::tilde(rhoip.data(), r);
		Matrix<double, 3, -1> Phi = Ai * PHI_iP;
		Matrix<double, 3, -1> Psi = Ai * PSI_iP;
		Ui.block(0, 0, 3, deli) = hh + r * ho;
		Ui.block(3, 0, 3, deli) = ho;
		Ui.block(0, deli, 3, si) = -Phi - r * Psi;
		Ui.block(3, deli, 3, si) = -Psi;
	}
	return Ui;
}
/*check once*/
VectorXd& JointBase::calbetai(IN double t, IN double* y, IN double* dy)
{
	// calculate beta2
	Vector3d wj = Bj->angularVel();
	Vector3d wi = Bi->angularVel();
	Vector3d wjQ = Bj->angularVel(NQ);
	Vector3d wiP = Bi->angularVel(NP);
	Vector3d wri = Wri();
	Vector3d beta2 = Yita();
	Matrix3d Ah0;
	calAh0(Ah0);
	wri = Ah0 * wri;
	beta2 = Ah0 * beta2;
	beta2 += wj.cross(wjQ - wj) + wjQ.cross(wri) - wi.cross(wiP - wi);
	betai.segment<3>(3) = beta2;
	//calculate beta1
	Vector3d rhojQ = Bj->rhoP(rhoj, NQ);
	Vector3d rhoiP = Bi->rhoP(rhoi, NP);
	Vector3d vriP = Bi->uiP(NP);
	Vector3d vrjQ = Bj->uiP(NQ);
	Vector3d vri = Ah0 * Vri();
	Vector3d hi = Ah0 * Hi();
	Vector3d beta1 = rhoiP.cross(beta2);
	beta1 += wj.cross(wj.cross(rhojQ) + 2 * vrjQ);
	beta1 += wjQ.cross(wjQ.cross(hi) + 2 * vri);
	beta1 -= wi.cross(wi.cross(rhoiP) + 2 * vriP);
	beta1 -= hi.cross(wj.cross(wjQ - wj));
	betai.segment<3>(0) = beta1;
	return betai;
}
/*check once*/
MatrixXd& JointBase::calTij(IN double t, IN double* y)
{
	// calculate Tij
	unsigned int si = Bi->nMode();
	unsigned int sj = Bj->nMode();
	Vector3d rhoip = rhoi;
	Vector3d rhojq = rhoj;
	Vector3d hi = Hi();
	Matrix3d Ai, Aj, Ah0;
	calAh0(Ah0);
	Bi->A(Ai);
	Bj->A(Aj);
	hi = Ah0 * hi;
	/*calculate rhoip and rhojq*/
	if (Bi->type() == Body::FLEXIBLE)
	{
		Map<VectorXd> ai(Bi->pos + Body::NC, si);
		auto& PHI_iP = static_cast<FlexibleBody*>(Bi)->PHI[NP];
		rhoip += PHI_iP * ai;
	}
	rhoip = Ai * rhoip;
	if (Bj->type() == Body::FLEXIBLE)
	{
		Map<VectorXd> aj(Bj->pos + Body::NC, sj);
		auto& PHI_jQ = static_cast<FlexibleBody*>(Bj)->PHI[NQ];
		rhojq += PHI_jQ * aj;
	}
	rhojq = Aj * rhojq;

	Matrix3d ri, rj, hitil;
	AUX::tilde(rhoip.data(), ri);
	AUX::tilde(rhojq.data(), rj);
	AUX::tilde(hi.data(), hitil);
	Tij.block<3, 3>(0, 3) = ri - rj - hitil;
	if (sj != 0)
	{
		auto& PHI_jQ = static_cast<FlexibleBody*>(Bj)->PHI[NQ];
		auto& PSI_jQ = static_cast<FlexibleBody*>(Bj)->PSI[NQ];
		MatR3CX Phi = Aj * PHI_jQ;
		MatR3CX Psi = Aj * PSI_jQ;
		Tij.block(0, 6, 3, sj) = Phi + ri * Psi - hitil * Psi;
		Tij.block(3, 6, 3, sj) = Psi;
	}
	return Tij;
}
/*check once*/
bool JointBase::calytoq(IN double t, IN double* y)
{
	Matrix3d Ai, Aj, Ah0;
	calAh0(Ah0);
	Ai=Ah0* calDih() * calCiP().transpose() * calBiP().transpose();
	Bj->A(Aj);
	//calculate rotation
	if (Body::m_s_rtype == RCORDS::EULERQUATERNION)
		AUX::AtoEQ(Ai, Bi->pos + 3);
	else if (Body::m_s_rtype == RCORDS::EULERANGLE)
		AUX::AtoEA(Ai, Bi->pos + 3);
	else if (Body::m_s_rtype == RCORDS::CARDANANGLE)
		AUX::AtoCA(Ai, Bi->pos + 3);
	//calculate translation
	Vector3d rhoip=rhoi, rhojq=rhoj;
	if (Bi->type() == Body::FLEXIBLE)
	{
		Map<VectorXd> ai(yi + DOF(), Bi->nMode());
		auto& PHI_iP = static_cast<FlexibleBody*>(Bi)->PHI[NP];
		rhoip += PHI_iP * ai;
	}
	rhoip = Ai * rhoip;
	if (Bj->type() == Body::FLEXIBLE)
	{
		Map<VectorXd> aj(Bj->pos + Body::NC, Bj->nMode());
		auto& PHI_jQ = static_cast<FlexibleBody*>(Bj)->PHI[NQ];
		rhojq += PHI_jQ * aj;
	}
	rhojq = Aj * rhojq;
	Map<Vector3d> rj(Bj->pos);
	Vector3d ri = rj - rhoip + rhojq + Ah0 * Hi();
	for (int i = 0; i < 3; ++i)
		Bi->pos[i] = ri(i);
	if (Bi->type() == Body::FLEXIBLE)
	{
		unsigned int nc = Body::NC;
		unsigned int del = DOF();
		for (unsigned int i = 0; i < Bi->nMode(); ++i)
		{
			Bi->pos[nc + i] = yi[del + i];
		}
	}
	return true;
}
/*check once*/
bool JointBase::caldytodq(IN double t, IN double* dy)
{
	/*ensuring that Tij and Ui have been already updated.*/
	unsigned int si = Bi->nMode();
	unsigned int sj = Bj->nMode();
	VectorXd vi;
	Map<VectorXd> y(dyi, DOF() + si);
	Map<VectorXd> vj(Bj->vel, 7 + sj);
	if (Body::m_s_rtype==RCORDS::EULERQUATERNION)
	{//Euler Quaternion
		MatrixXd Kj = MatrixXd::Zero(6 + sj, 7 + sj);
		Kj.block<3, 3>(0, 0).setIdentity();
		Kj.block(6, 7, sj, sj).setIdentity();
		MatR3C4 Rj, Ri;
		AUX::R(Bj->pos + 3, Rj);
		AUX::R(Bi->pos + 3, Ri);
		Kj.block<3, 4>(3, 3) = 2 * Rj;
		vi = Tij * Kj * vj + Ui * y;
		for (unsigned int i = 0; i < 3; ++i)
			Bi->vel[i] = vi(i);
		Vector4d lami = 0.5 * Ri.transpose() * vi.segment<3>(3);
		for (unsigned int i = 0; i < 4; ++i)
			Bi->vel[i + 3] = lami[i];
		for (unsigned int i = 0; i < si; ++i)
			Bi->vel[7 + i] = vi(6 + i);
	}
	else if (Body::m_s_rtype==RCORDS::EULERANGLE)
	{//Euler Angle
		MatrixXd Kj = MatrixXd::Zero(6 + sj, 6 + sj);
		Kj.block<3, 3>(0, 0).setIdentity();
		Kj.block(6, 6, sj, sj).setIdentity();
		MatR3CX Krj(3, 3), Kri(3, 3);
		AUX::Kr(Bj->pos + 3, RCORDS::EULERANGLE,Krj);
		AUX::Kr(Bi->pos + 3, RCORDS::EULERANGLE, Kri);
		Kj.block<3, 3>(3, 3) = Krj;
		vi = Tij * Kj * vj + Ui * y;
		for (unsigned int i = 0; i < 3; ++i)
			Bi->vel[i] = vi(i);
		Vector3d lami = Kri.partialPivLu().solve(vi.segment<3>(3));
		for (unsigned int i = 0; i < 3; ++i)
			Bi->vel[i + 3] = lami[i];
		for (unsigned int i = 0; i < si; ++i)
			Bi->vel[6 + i] = vi(6 + i);
	}
	else if (Body::m_s_rtype == RCORDS::CARDANANGLE)
	{//Cardan Angle
		MatrixXd Kj = MatrixXd::Zero(6 + sj, 6 + sj);
		Kj.block<3, 3>(0, 0).setIdentity();
		Kj.block(6, 6, sj, sj).setIdentity();
		MatR3CX Krj(3, 3), Kri(3, 3);
		AUX::Kr(Bj->pos + 3, RCORDS::CARDANANGLE, Krj);
		AUX::Kr(Bi->pos + 3, RCORDS::CARDANANGLE, Kri);
		Kj.block<3, 3>(3, 3) = Krj;
		vi = Tij * Kj * vj + Ui * y;
		for (unsigned int i = 0; i < 3; ++i)
			Bi->vel[i] = vi(i);
		Vector3d lami = Kri.partialPivLu().solve(vi.segment<3>(3));
		for (unsigned int i = 0; i < 3; ++i)
			Bi->vel[i + 3] = lami[i];
		for (unsigned int i = 0; i < si; ++i)
			Bi->vel[6 + i] = vi(6 + i);
	}
	return true;
}
/*check once*/
bool JointBase::setCiP(const Matrix3d& c)
{
	CiP = c;
	return true;
}
/*check once*/
bool JointBase::setCjQ(const Matrix3d& c)
{
	CjQ = c;
	return true;
}
/*check once*/
MatR3CX Universe::HhT() const
{
	return MatR3CX::Zero(3,2);
}
/*check once*/
MatR3CX Universe::HoT() const
{
	MatR3CX m = MatR3CX::Zero(3, 2);
	m(0, 0) = 1;
	m(1, 1) = cos(yi[0]);
	m(2, 1) = sin(yi[0]);
	return m;
}
/*check once*/
Vector3d Universe::Yita() const
{
	Vector3d v(0, 0, 0);
	v(1) = -sin(yi[0]) * dyi[0] * dyi[1];
	v(2) = cos(yi[0]) * dyi[0] * dyi[1];
	return v;
}
/*check once*/
Matrix3d& Universe::calDih()
{
	//assume rotate along x-axis and y-axis
	double s1 = sin(yi[0]), c1 = cos(yi[0]);
	double s2 = sin(yi[1]), c2 = cos(yi[1]);
	Dih(0, 0) = c2;
	Dih(0, 1) = 0;
	Dih(0, 2) = s2;
	Dih(1, 0) = s1 * s2;
	Dih(1, 1) = c1;
	Dih(1, 2) = -s1 * c2;
	Dih(2, 0) = -c1 * s2;
	Dih(2, 1) = s1;
	Dih(2, 2) = c1 * c2;
	return Dih;
}
/*check once*/
Universe::Universe(Body* Bi_ptr, Body* Bj_ptr, Vector3d& rho_i, Vector3d& rho_j) :JointBase(Bi_ptr, Bj_ptr, rho_i, rho_j)
{
	unsigned int si = Bi_ptr->nMode();
	Ui.resize(6 + si, DOF() + si);
	Ui.setZero();
	if (si != 0)
	{
		Ui.block(6, DOF(), si, si).setIdentity();
	}
}

Universe::~Universe()
{
}

inline unsigned int Universe::DOF() const
{
	return 2;
}
/*check once*/
inline unsigned int Universe::type() const
{
	return JointBase::UNIVERSAL;
}
/*check once*/
MatR3CX Virtual::HhT() const
{
	MatR3CX hht = MatR3CX::Zero(3, 6);
	hht.block<3, 3>(0, 0).setIdentity();
	return hht;
}
/*check once*/
MatR3CX Virtual::HoT() const
{
	MatR3CX hot = MatR3CX::Zero(3, 6);
	double s1 = sin(yi[3]), c1 = cos(yi[3]), s2 = sin(yi[4]), c2 = cos(yi[4]);
	hot(0, 3) = 1;
	hot(0, 5) = s2;
	hot(1, 4) = c1;
	hot(1, 5) = -c2 * s1;
	hot(2, 4) = s1;
	hot(2, 5) = c2 * c1;
	return hot;
}
/*check once*/
Vector3d Virtual::Yita() const
{
	Vector3d ans(0, 0, 0);
	double s1 = sin(yi[3]), c1 = cos(yi[3]);
	double s2 = sin(yi[4]), c2 = cos(yi[4]);
	double dq1 = dyi[3], dq2 = dyi[4], dq3 = dyi[5];
	ans(0) = c1 * dq2 * dq3;
	ans(1) = -s1 * dq1 * dq2 - c1 * c2 * dq1 * dq3 + s1 * s2 * dq2 * dq3;
	ans(2) = c1 * dq1 * dq2 - s1 * c2 * dq1 * dq3 - c1 * s2 * dq2 * dq3;
	return ans;
}
/*check once*/
Matrix3d& Virtual::calDih()
{
	double s1 = sin(yi[3]), c1 = cos(yi[3]);
	double s2 = sin(yi[4]), c2 = cos(yi[4]);
	double s3 = sin(yi[5]), c3 = cos(yi[5]);
	Dih(0, 0) = c2 * c3;
	Dih(0, 1) = -c2 * s3;
	Dih(0, 2) = s2;
	Dih(1, 0) = s1 * s2 * c3 + c1 * s3;
	Dih(1, 1) = -s1 * s2 * s3 + c1 * c3;
	Dih(1, 2) = -s1 * c2;
	Dih(2, 0) = -c1 * s2 * c3 + s1 * s3;
	Dih(2, 1) = c1 * s2 * s3 + s1 * c3;
	Dih(2, 2) = c1 * c2;
	return Dih;
}
/*check once*/
Virtual::Virtual(Body* Bi_ptr, Body* Bj_ptr, Vector3d& rho_i, Vector3d& rho_j) :JointBase(Bi_ptr, Bj_ptr, rho_i, rho_j)
{
	unsigned int si = Bi_ptr->nMode();
	Ui.resize(6 + si, DOF() + si);
	Ui.setZero();
	if (si != 0)
	{
		Ui.block(6, DOF(), si, si).setIdentity();
	}
}

Virtual::~Virtual()
{
}

inline unsigned int Virtual::DOF() const
{
	return 6;
}

inline unsigned int Virtual::type() const
{
	return JointBase::VIRTUAL;
}

MatR3CX Sphere::HhT() const
{
	return MatR3CX::Zero(3,3);
}
/*check once*/
MatR3CX Sphere::HoT() const
{
	MatR3CX hot = MatR3CX::Zero(3, 3);
	double s1 = sin(yi[0]), c1 = cos(yi[0]);
	double s2 = sin(yi[1]), c2 = cos(yi[1]);
	hot(0, 0) = 1;
	hot(0, 2) = s2;
	hot(1, 1) = c1;
	hot(1, 2) = -c2 * s1;
	hot(2, 1) = s1;
	hot(2, 2) = c2 * c1;
	return hot;
}
/*check once*/
Vector3d Sphere::Yita() const
{
	Vector3d ans(0, 0, 0);
	double s1 = sin(yi[0]), c1 = cos(yi[0]);
	double s2 = sin(yi[1]), c2 = cos(yi[1]);
	double dq1 = dyi[0], dq2 = dyi[1], dq3 = dyi[2];
	ans(0) = c1 * dq2 * dq3;
	ans(1) = -s1 * dq1 * dq2 - c1 * c2 * dq1 * dq3 + s1 * s2 * dq2 * dq3;
	ans(2) = c1 * dq1 * dq2 - s1 * c2 * dq1 * dq3 - c1 * s2 * dq2 * dq3;
	return ans;
}
/*check once*/
Matrix3d& Sphere::calDih()
{
	double s1 = sin(yi[0]), c1 = cos(yi[0]);
	double s2 = sin(yi[1]), c2 = cos(yi[1]);
	double s3 = sin(yi[2]), c3 = cos(yi[2]);
	Dih(0, 0) = c2 * c3;
	Dih(0, 1) = -c2 * s3;
	Dih(0, 2) = s2;
	Dih(1, 0) = s1 * s2 * c3 + c1 * s3;
	Dih(1, 1) = -s1 * s2 * s3 + c1 * c3;
	Dih(1, 2) = -s1 * c2;
	Dih(2, 0) = -c1 * s2 * c3 + s1 * s3;
	Dih(2, 1) = c1 * s2 * s3 + s1 * c3;
	Dih(2, 2) = c1 * c2;
	return Dih;
}
/*check once*/
Sphere::Sphere(Body* Bi_ptr, Body* Bj_ptr, Vector3d& rho_i, Vector3d& rho_j) :JointBase(Bi_ptr, Bj_ptr, rho_i, rho_j)
{
	unsigned int si = Bi_ptr->nMode();
	Ui.resize(6 + si, DOF() + si);
	Ui.setZero();
	if (si != 0)
	{
		Ui.block(6, DOF(), si, si).setIdentity();
	}
}

Sphere::~Sphere()
{
}

inline unsigned int Sphere::DOF() const
{
	return 3;
}

inline unsigned int Sphere::type() const
{
	return JointBase::SPHERICAL;
}
/*check once*/
Matrix<double, 3, -1> Prism::HhT() const
{
	return Vector3d(1, 0, 0);
}
/*check once*/
Matrix<double, 3, -1> Prism::HoT() const
{
	return Vector3d(0, 0, 0);
}
/*check once*/
Vector3d Prism::Yita() const
{
	return Vector3d(0,0,0);
}
/*check once*/
Matrix3d& Prism::calDih()
{
	Dih.setIdentity();
	return Dih;
}
/*check once*/
Prism::Prism(Body* Bi_ptr, Body* Bj_ptr, Vector3d& rho_i, Vector3d& rho_j) :JointBase(Bi_ptr, Bj_ptr, rho_i, rho_j)
{
	unsigned int si = Bi_ptr->nMode();
	Ui.resize(6 + si, DOF() + si);
	Ui.setZero();
	if (si != 0)
	{
		Ui.block(6, DOF(), si, si).setIdentity();
	}
}

Prism::~Prism()
{
}

inline unsigned int Prism::DOF() const
{
	return 1;
}

inline unsigned int Prism::type() const
{
	return JointBase::PRISMATIC;
}
/*check once*/
MatR3CX Cylinder::HhT() const
{
	MatR3CX hht = MatR3CX::Zero(3, 2);
	hht(0, 0) = 1;
	return hht;
}
/*check once*/
MatR3CX Cylinder::HoT() const
{
	MatR3CX hot = MatR3CX::Zero(3, 2);
	hot(0, 1) = 1;
	return hot;
}
/*check once*/
Vector3d Cylinder::Yita() const
{
	return Vector3d(0,0,0);
}
/*check once*/
Matrix3d& Cylinder::calDih()
{
	double s1 = sin(yi[1]), c1 = cos(yi[1]);
	Dih(1, 1) = c1;
	Dih(1, 2) = -s1;
	Dih(2, 1) = s1;
	Dih(2, 2) = c1;
	return Dih;
}
/*check once*/
Cylinder::Cylinder(Body* Bi_ptr, Body* Bj_ptr, Vector3d& rho_i, Vector3d& rho_j) :JointBase(Bi_ptr, Bj_ptr, rho_i, rho_j)
{
	unsigned int si = Bi_ptr->nMode();
	Ui.resize(6 + si, DOF() + si);
	Ui.setZero();
	if (si != 0)
	{
		Ui.block(6, DOF(), si, si).setIdentity();
	}
}

Cylinder::~Cylinder()
{
}

inline unsigned int Cylinder::DOF() const
{
	return 2;
}

inline unsigned int Cylinder::type() const
{
	return JointBase::CYLINDRICAL;
}
