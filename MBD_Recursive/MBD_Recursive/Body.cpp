#include "Body.h"

RCORDS Body::m_s_rtype = RCORDS::EULERQUATERNION;
unsigned int Body::NC = 7;

Body::Body(double m) noexcept
{
	tM = m;
	id = -1;             //  id not assigned
}

Body::Body() noexcept
{
	tM = 0.;
	id = -1;            //    id not assigned
}



bool Body::operator==(Body& other)const
{
	if (this->id == other.id)
		return true;
	else
		return false;
}

bool Body::operator<(Body& other)const
{
	if (this->id < other.id)
		return true;
	else
		return false;
}

bool Body::operator<=(Body& other)const
{
	return (*this) < other || (*this) == other;
}

bool Body::operator>(Body& other)const
{
	return !((*this) <= other);
}

bool Body::setID(int i)
{
	id = i;
	return true;
}

/*check once*/
bool Body::A(Matrix3d& M) const
{
	if (m_s_rtype == RCORDS::EULERQUATERNION)
		AUX::A(pos + 3, M);
	else if (m_s_rtype == RCORDS::EULERANGLE)
		AUX::EAtoA(pos + 3, M);
	else if (m_s_rtype == RCORDS::CARDANANGLE)
		AUX::CAtoA(pos + 3, M);
	return true;
}
/*check once*/
Vector3d Body::angularVel() const
{
	Vector3d w;
	if (m_s_rtype==RCORDS::EULERQUATERNION)
	{//Euler Quaternion
		Matrix<double, 3, 4> r;
		Map<Vector4d> lambda(vel + 3);
		AUX::R(pos + 3, r);
		w = 2 * r * lambda;
	}
	else if (m_s_rtype == RCORDS::CARDANANGLE)
	{
		MatR3CX Kr(3, 3);
		AUX::Kr(pos + 3, RCORDS::CARDANANGLE, Kr);
		Map<Vector3d> dq(vel + 3);
		w = Kr * dq;
	}
	else if (m_s_rtype == RCORDS::EULERANGLE)
	{
		MatR3CX Kr(3, 3);
		AUX::Kr(pos + 3, RCORDS::EULERANGLE, Kr);
		Map<Vector3d> dq(vel + 3);
		w = Kr * dq;
	}
	return w;
}
/*check once*/
Vector3d Body::angularVel(int EID) const
{
	return angularVel();
}
/*check once*/
Vector3d Body::rhoP(Vector3d& _rp, int EID) const
{
	Matrix3d T;
	A(T);
	return T * _rp;
}
/*check once*/
Vector3d Body::translationalVel(Vector3d& _rp) const
{
	Vector3d v(vel[0], vel[1], vel[2]);
	v += angularVel().cross(Body::rhoP(_rp, -1));
	return v;
}
/*check once*/
Vector3d Body::translationalVel(Vector3d& _rp, int EID) const
{
	return Body::translationalVel(_rp);
}
/*check once*/
VectorXd Body::inertiaForce()
{
	return VectorXd::Zero(6);
}
/*check once*/
Body::~Body()
{
}

bool Body::Write(Json::Value& body) const
{
	body["Id"] = Json::Value(id);
	return true;
}
/*check once*/
Vector3d Body::uiP(int EID) const
{
	return Vector3d(0,0,0);
}

Matrix3d& RigidBody::calM11()
{
	return M11;
}

Matrix3d& RigidBody::calM12()
{
	return M12;
}
/*check once*/
Matrix3d& RigidBody::calM22()
{
	Matrix3d A;
	Body::A(A);
	M22 = A * Jc * A.transpose();
	return M22;
}

RigidBody::RigidBody()
{
	M11.setZero();
	M12.setZero();
	InerF = VectorXd::Zero(6);
}

RigidBody::RigidBody(double m):Body(m)
{
	M11 = Matrix3d::Identity() * m;
	M12.setZero();
	InerF = VectorXd::Zero(6);
}

RigidBody::RigidBody(double m, Matrix3d& I):Body(m)
{
	M11 = Matrix3d::Identity() * m;
	M12.setZero();
	Jc = I;
	InerF = VectorXd::Zero(6);
}




unsigned int RigidBody::type() const
{
	return Body::RIGID;
}
/*check once*/
bool RigidBody::calMass(MatrixXd& M, int k)
{
	M.block<3, 3>(k, k) = calM11();
	M.block<3, 3>(k, k + 3) = calM12();
	M.block<3, 3>(k + 3, 3) = M12.transpose();
	M.block<3, 3>(k + 3, k + 3) = calM22();
	return true;
}

unsigned int RigidBody::nMode() const
{
	return 0;
}
/*check once*/
VectorXd RigidBody::inertiaForce()
{
	Vector3d w = angularVel();
	InerF.segment<3>(3) = -w.cross(M22 * w);
	return InerF;
}

RigidBody::~RigidBody()
{

}

bool RigidBody::Write(Json::Value& body) const
{
	Body::Write(body);
	body["Type"] = Json::Value("Rigid");
	body["Mass"] = Json::Value(tM);
	for (int i = 0; i < 3; ++i)
		for (int j = i; j < 3; ++j)
			body["Jc"].append(Jc(i, j));
	return true;
}

Matrix3d& FlexibleBody::calM11()
{
	// TODO: insert return statement here
	return M11;
}
/*check once*/
Matrix3d& FlexibleBody::calM12()
{
	// TODO: insert return statement here
	Matrix3d T;
	this->A(T);
	Matrix3d m;
	Vector3d v(0, 0, 0);
	Map<VectorXd> a(pos + NC, s);
	v = g2 * a + g1;
	AUX::tilde(v.data(), m);
	M12 = -T * m * T.transpose();
	return M12;
}

Matrix3d& FlexibleBody::calM22()
{
	Matrix3d T;
	this->A(T);
	calG5();
	Matrix3d m = G1;
	for (unsigned int i = 0; i < s; ++i)
	{
		m += G5[i] * pos[NC + i];
	}
	M22 = -T * m * T.transpose();
	return M22;
}
/*check once*/
Matrix<double,3,-1>& FlexibleBody::calM13()
{
	Matrix3d T;
	this->A(T);
	M13 = T * g2;
	return M13;
}
/*check once*/
Matrix<double,3,-1>& FlexibleBody::calM23()
{
	// calculate gama5 firstly.
	calg5();
	Matrix3d T;
	this->A(T);
	M23 = T * g5;
	return M23;
}
/*check once*/
MatrixXd& FlexibleBody::calM33()
{
	return M33;
}
/*check once*/
Vector3d FlexibleBody::uiP(int EID) const
{
	if (EID < 0)
		return Body::uiP(EID);
	Map<VectorXd> dai(vel + Body::NC, nMode());
	Matrix3d T;
	A(T);
	return T * PHI[EID] * dai;
}
/*check once*/
bool FlexibleBody::calg5()
{
	for (unsigned int j = 0; j < s; ++j)
	{
		g5.col(j) = g3.col(j);
		for (unsigned int i = 0; i < s; ++i)
		{
			g5.col(j) += g4[i][j] * pos[NC + i];
		}
	}
	return true;
}
/*check once*/
bool FlexibleBody::calG4()
{
	for (unsigned int j = 0; j < s; ++j)
	{
		G4[j] = G2[j];
		for (unsigned int i = 0; i < s; ++i)
		{
			G4[j] += G3[i][j] * pos[NC + i];
		}
	}
	return true;
}
/*check once*/
bool FlexibleBody::calG5()
{
	//calculate Gama4 firstly.
	calG4();
	for (unsigned int j = 0; j < s; ++j)
	{
		G5[j] = G2[j].transpose() + G4[j];
	}
	return true;
}
/*check once*/
FlexibleBody::FlexibleBody(int NE, int nmode)
{
	s = nmode;
	ne = NE;
	me.resize(ne);
	rho.resize(3, ne);
	PHI.resize(ne);
	PSI.resize(ne);
	Ka.resize(s, s);
	Ca.resize(s, s);
	M13.resize(3, s);
	M23.resize(3, s);
	M33.resize(s, s);
	g2.resize(3, s);
	g3.resize(3, s);
	g4.resize(s, vector<Vector3d>(s, Vector3d(0,0,0)));
	g5.resize(3, s);
	G2.resize(s, Matrix3d::Zero());
	G3.resize(s, vector<Matrix3d>(s, Matrix3d::Zero()));
	G4.resize(s, Matrix3d::Zero());
	G5.resize(s, Matrix3d::Zero());
	InerF.resize(6 + s);
}

bool FlexibleBody::setMe(VectorXd& Me)
{
	if (Me.rows() != ne)
		throw MBException(INI_FAILURE("FlexibleBody::setMe"));
	me = Me;
	tM = me.sum();
	return true;
}

/*
* me is set before rho is set
* when rho is set, then gama1 
* and Gama1 can be calculated
* at the same time.
*/
/*check once*/
bool FlexibleBody::setRho(MatrixXd& Rho)
{
	if (Rho.rows() != 3 || Rho.cols() != ne)
	{
		throw MBException(INI_FAILURE("FlexibleBody::setRho"));
	}
	rho = Rho;
	//calculate gama1
	g1 = rho * me;
	//calculate Gama1
	G1.setZero();
	Matrix3d R=Matrix3d::Zero();
	for (unsigned int k = 0; k < ne; ++k)
	{
		AUX::tilde(rho.col(k).data(), R);
		G1 += me[k] * R * R;
	}
	return true;
}

/*
* PHI is set after me and rho are set,
* when PHI is set, then gama2, gama3, gama4,
* Gama2, Gama3 can be calculated.
*/
/*check once*/
bool FlexibleBody::setPHI(vector<MatrixXd>& phi)
{
	if (me.rows() != ne || rho.rows() != 3 || rho.cols() != ne)
	{
		throw MBException(INI_FAILURE("FlexibleBody::setPHI"));
	}
	if (phi.size() != ne)
	{
		throw MBException(INI_FAILURE("FlexibleBody::setPHI"));
	}
	for (unsigned int i = 0; i < ne; ++i)
	{
		if (phi[i].cols() != s)
		{
			throw MBException(INI_FAILURE("FlexibleBody::setPHI"));
		}
		PHI[i] = phi[i];
	}
	//calculate gama2
	g2.setZero();
	for (unsigned int j = 0; j < s; ++j)
	{
		for (unsigned int k = 0; k < ne; ++k)
		{
			g2.col(j) += me(k) * PHI[k].col(j);
		}
	}
	//calculate gama3
	g3.setZero();
	for (unsigned int j = 0; j < s; ++j)
	{
		for (unsigned int k = 0; k < ne; ++k)
		{
			g3.col(j) += me(k) * rho.col(k).cross(PHI[k].col(j));
		}
	}
	//calculate gama4
	for (unsigned int i = 0; i < s; ++i)
	{
		for (unsigned int j = i; j < s; ++j)
		{
			g4[i][j].setZero();
			for (unsigned int k = 0; k < ne; ++k)
			{
				g4[i][j] += me(k) * PHI[k].col(i).cross(PHI[k].col(j));
			}
			g4[j][i] = -g4[i][j];
		}
	}
	//calculate Gama2
	Matrix3d rtil;
	Matrix3d ptil;
	for (unsigned int j = 0; j < s; ++j)
	{
		G2[j].setZero();
		for (unsigned int k = 0; k < ne; ++k)
		{
			AUX::tilde(rho.col(k).data(), rtil);
			AUX::tilde(PHI[k].col(j).data(), ptil);
			G2[j] += me(k) * rtil * ptil;
		}
	}
	//calculate Gama3
	for (unsigned int i = 0; i < s; ++i)
	{
		for (unsigned int j = i; j < s; ++j)
		{
			G3[i][j].setZero();
			for (unsigned int k = 0; k < ne; ++k)
			{
				AUX::tilde(PHI[k].col(i).data(), rtil);
				AUX::tilde(PHI[k].col(j).data(), ptil);
				G3[i][j] += me(k) * rtil * ptil;
			}
			if (i != j)
				G3[j][i] = G3[i][j].transpose();
		}
	}
	//calculate M33
	for (unsigned int i = 0; i < s; ++i)
	{
		for (unsigned int j = i; j < s; ++j)
		{
			M33(i, j) = 0;
			for (unsigned int k = 0; k < ne; ++k)
			{
				M33(i, j) += me(k) * PHI[k].col(i).dot(PHI[k].col(j));
			}
			M33(j, i) = M33(i, j);
		}
	}
	return true;
}
/*check once*/
bool FlexibleBody::setPSI(vector<MatrixXd>& psi)
{
	if (psi.size() != ne)
		throw MBException(INI_FAILURE("Flexible::setPSI"));
	for (unsigned int i = 0; i < ne; ++i)
	{
		if (psi[i].cols() != s)
			throw MBException(INI_FAILURE("FlexibleBody::setPSI"));
		PSI[i] = psi[i];
	}
	return true;
}
/*check once*/
bool FlexibleBody::setKa(MatrixXd& ka)
{
	if (ka.rows() != s || ka.cols() != s)
		throw MBException(INI_FAILURE("FlexibleBody::setKa"));
	Ka = ka;
	return true;
}
/*check once*/
bool FlexibleBody::setCa(MatrixXd& ca)
{
	if (ca.rows() != s || ca.cols() != s)
		throw MBException(INI_FAILURE("FlexibleBody::setCa"));
	Ca = ca;
	return true;
}


unsigned int FlexibleBody::type() const
{
	return Body::FLEXIBLE;
}
/*check once*/
bool FlexibleBody::calMass(MatrixXd& M, int k)
{
	M.block<3, 3>(k, k) = calM11();
	M.block<3, 3>(k + 3, k + 3) = calM22();
	M.block<3, 3>(k, k + 3) = calM12();
	M.block(k, k + 6, 3, s) = calM13();
	M.block(k + 3, k + 6, 3, s) = calM23();
	M.block(k + 6, k + 6, s, s) = calM33();
	M.block<3,3>(k + 3, k) = M12.transpose();
	M.block(k + 6, k, s, 3) = M13.transpose();
	M.block(k + 6, k + 3, s, 3) = M23.transpose();
	return true;
}

unsigned int FlexibleBody::nMode() const
{
	return s;
}
/*check once*/
Vector3d FlexibleBody::angularVel(int EID) const
{
	Vector3d wiP = Body::angularVel();
	if (EID < 0)
		return wiP;
	Map<VectorXd> dai(vel + Body::NC, nMode());
	Matrix3d Ai;
	Body::A(Ai);
	wiP += Ai * PSI[EID] * dai;
	return wiP;
}
/*check once*/
Vector3d FlexibleBody::rhoP(Vector3d& _rp, int EID) const
{
	Matrix3d T;
	A(T);
	if (EID < 0)
		return T * _rp;
	Map<VectorXd> a(pos + Body::NC, nMode());
	return T * (_rp + PHI[EID] * a);
}
/*check once*/
Vector3d FlexibleBody::translationalVel(Vector3d& _rp, int EID) const
{
	Vector3d v(vel[0], vel[1], vel[2]);
	v += Body::angularVel().cross(rhoP(_rp, EID)) + uiP(EID);
	return v;
}
/*check once*/
VectorXd FlexibleBody::inertiaForce()
{
	Matrix3d A;
	Body::A(A);
	Map<VectorXd> a(pos + Body::NC, s);
	Map<VectorXd> da(vel + Body::NC, s);
	Vector3d w1, w2;
	Vector3d w = Body::angularVel();
	//calculate w1
	w1 = w.cross(w.cross(A * (g1 + g2 * a)));
	w1 += 2 * w.cross(A * (g2 * da));
	//calculate w2
	w2 = w.cross(M22 * w);
	Matrix3d temp = Matrix3d::Zero();
	for (unsigned int i = 0; i < s; ++i)
	{
		temp += G4[i] * da(i);
	}
	w2 -= 2 * A * temp * A.transpose() * w;
	//calculate w3
	VectorXd w3 = VectorXd::Zero(s);
	Vector3d v1 = A.transpose() * w;
	Vector3d temp1;
	for (unsigned int i = 0; i < s; ++i)
	{
		temp1.setZero();
		for (unsigned int j = 0; j < s; ++j)
		{
			temp1 += g4[j][i] * a(j);
		}
		w3(i) = 2 * temp1.dot(v1) + v1.dot(G4[i] * v1);
	}
	w3 += Ka * a + Ca * da;
	//InerF
	InerF.head<3>() = -w1;
	InerF.segment<3>(3) = -w2;
	InerF.tail(s) = -w3;

	return InerF;
}

FlexibleBody::~FlexibleBody()
{

}
bool FlexibleBody::Write(Json::Value& body) const
{
	Body::Write(body);
	body["Type"] = Json::Value("Flexible");
	string fname = "flex_" + to_string(id) + string(".txt");
	body["ConfigFile"] = Json::Value(fname);
	fstream fout(fname, ios::out);
	if (!fout)
		throw MBException("Open Modal Data File Error!");
	fout << ne << endl;
	fout << s << endl;
	for (unsigned int i = 0; i < ne; ++i)
		fout << me(i) << endl;
	for (unsigned int i = 0; i < ne; ++i)
		fout << rho(0, i) << '\t' << rho(1, i) << '\t' << rho(2, i) << endl;
	for (unsigned int k = 0; k < ne; ++k)
	{
		for (unsigned int i = 0; i < s; ++i)
		{
			fout << PHI[k](0, i) << '\t' << PHI[k](1, i) << '\t' << PHI[k](2, i) << '\t';
			fout << PSI[k](0, i) << '\t' << PSI[k](1, i) << '\t' << PSI[k](2, i) << endl;
		}
	}
	for (unsigned int i = 0; i < s; ++i)
	{
		for (unsigned int j = 0; j < s; ++j)
			fout << Ka(i, j) << '\t';
		fout << endl;
	}
	for (unsigned int i = 0; i < s; ++i)
	{
		for (unsigned int j = 0; j < s; ++j)
			fout << Ca(i, j) << '\t';
		fout << endl;
	}
	fout.close();
	return true;
}
/*check once*/
BaseBody::BaseBody()
{
	tM = 0;
	id = 0;    //BaseBody id is definitely zero
	pos = new double[Body::NC];
	vel = new double[Body::NC];
	acc = new double[Body::NC];
	for (unsigned int i = 0; i < Body::NC; ++i)
	{
		pos[i] = 0;
		vel[i] = 0;
		acc[i] = 0;
	}
	if (Body::m_s_rtype == RCORDS::EULERQUATERNION)
		pos[3] = 1;
	else if (Body::m_s_rtype == RCORDS::EULERANGLE)
		pos[4] = PI;
}

BaseBody::BaseBody(VectorXd(*p)(double), VectorXd(*v)(double), VectorXd(*a)(double)):BaseBody()
{
	pFun = p;
	vFun = v;
	aFun = a;
}


BaseBody::~BaseBody()
{
	delete[] pos;
	delete[] vel;
	delete[] acc;
}

bool BaseBody::Write(Json::Value& body) const
{
	Body::Write(body);
	body["Type"] = Json::Value("Base");
	return true;
}

unsigned int BaseBody::type() const
{
	return BASE;
}

bool BaseBody::calMass(MatrixXd& M, int k)
{
	return true;
}

unsigned int BaseBody::nMode() const
{
	return 0;
}

VectorXd BaseBody::acceleration(double t)
{
	VectorXd a(6+nMode());
	if (Body::m_s_rtype==RCORDS::EULERQUATERNION)
	{
		Matrix<double, 3, 4> R;
		AUX::R(pos + 3, R);
		for (unsigned int i = 0; i < 3; ++i)
			a(i) = acc[i];
		Map<Vector4d> ddlam(acc + 3);
		a.segment(3,3) = 2 * R * ddlam;
		for (unsigned int i = 0; i < nMode(); ++i)
			a(6 + i) = acc[7 + i];
	}
	else if (Body::m_s_rtype == RCORDS::CARDANANGLE)
	{
		a.resize(6+nMode());
		for (unsigned int i = 0; i < 6; ++i)
			a(i) = acc[i];
	}
	return a;
}
/*check once*/
bool BaseBody::update(double t)
{
	if (pFun!=nullptr)
	{
		VectorXd p = pFun(t);
		for (unsigned int i = 0; i < Body::NC; ++i)
			pos[i] = p(i);
	}
	if (vFun != nullptr)
	{
		VectorXd v = vFun(t);
		for (unsigned int i = 0; i < Body::NC; ++i)
			vel[i] = v(i);
	}
	if (aFun != nullptr)
	{
		VectorXd a = aFun(t);
		for (unsigned int i = 0; i < Body::NC; ++i)
			acc[i] = a(i);
	}
	return true;
}
