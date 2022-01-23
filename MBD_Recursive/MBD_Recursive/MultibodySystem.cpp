#include "MultibodySystem.h"


/*
* transform y into q
* examine norm(ynew-yprev) 
* to avoid calculating repeatedly
* when ynew is equal to yprev then return true,
* else then calculate y->q and update yprev and qprev
* and assign the coordinates for the bodies in the system
*/
bool TreeSystem::calytoq(double t, double* y)
{
	double* yi = y;
	for (unsigned int i = 0; i < nb; ++i)
	{
		bjvec[i]->jptr->calytoq(t,yi);
		yi += ydimvec[i];
	}
	return true;
}

bool TreeSystem::caldytodq(double t, double* dy)
{
	double* dyi = dy;
	for (unsigned int i = 0; i < nb; ++i)
	{
		bjvec[i]->jptr->caldytodq(t, dyi);
		dyi += ydimvec[i];
	}
	return true;
}
/*check once*/
TreeSystem::TreeSystem(vector<JointBase*>& jvec)
{
	Body* baseptr = nullptr;
	for (auto& item : jvec)
	{
		if (item->Bj->type()==Body::BASE)
		{//find base body
			assert(item->Bj->id == 0);
			baseptr = item->Bj;
			break;
		}
	}
	MBTree.addRoot(BJpair(baseptr));
	queue<BJpair> bjque;
	bjque.push(BJpair(baseptr));
	auto NJ = jvec.size();
	bool* visited = new bool[NJ];
	for (unsigned int i = 0; i < NJ; ++i)
		visited[i] = false;
	while (!bjque.empty())
	{
		auto p = bjque.front();
		bjque.pop();
		for (unsigned int i = 0; i < NJ; ++i)
		{
			if (visited[i])
				continue;
			if (jvec[i]->Bj == p.bptr)
			{
				visited[i] = true;
				bjque.push(BJpair(jvec[i]->Bi, jvec[i]));
				MBTree.addNode(p, BJpair(jvec[i]->Bi, jvec[i]));
			}
		}
	}
	delete[]visited;
	MBTree.getNodevec(bjvec);
	bjroot = *bjvec.begin();
	bjvec.erase(bjvec.begin());
	nb = bjvec.size();
	qdimvec = VectorXi(nb, 0);
	ydimvec = VectorXi(nb, 0);
	qdimcum = VectorXi(nb, 0);
	ydimcum = VectorXi(nb, 0);
	unsigned int delta = 0;
	for (unsigned int i = 0; i < nb; ++i)
	{
		qdimvec[i] = bjvec[i]->bptr->nMode() + Body::NC;
		ydimvec[i] = bjvec[i]->jptr->DOF()+bjvec[i]->bptr->nMode();
		if (i != 0)
		{
			qdimcum[i] = qdimcum[i - 1] + qdimvec[i - 1];
			ydimcum[i] = ydimcum[i - 1] + ydimvec[i - 1];
		}
	}
	delta = ydimcum[nb - 1] + ydimvec[nb - 1];
	yprev.resize(delta);
	dyprev.resize(delta);
	unsigned int nc = qdimcum[nb - 1] + qdimvec[nb - 1];
	qprev.resize(nc);
	dqprev.resize(nc);
	//allocate address
	double* yi = yprev.data();
	double* dyi = dyprev.data();
	double* qi = qprev.data();
	double* dqi = dqprev.data();
	for (unsigned int i = 0; i < nb; ++i)
	{
		bjvec[i]->bptr->pos = qi;
		bjvec[i]->bptr->vel = dqi;
		bjvec[i]->jptr->yi = yi;
		bjvec[i]->jptr->dyi = dyi;
		qi += qdimvec[i];
		dqi += qdimvec[i];
		yi += ydimvec[i];
		dyi += ydimvec[i];
	}
}

bool TreeSystem::calG0(double t, double* y, OUT MatrixXd& G0)
{
	int I, J;
	Body * pbj;
	int brow, bcol = 0;
	int cdim = bjroot->bptr->nMode() + 6;
	int rdim;
	for (unsigned int i = 0; i < nb; ++i)
	{
		I = bjvec[i]->bptr->id - 1;
		pbj = bjvec[i]->jptr->Bj;
		brow = qdimcum[I] - (Body::NC - 6)*I;                //q����ŷ����Ԫ��������̬�� 7-1=6
		rdim = qdimvec[I] - (Body::NC - 6);
		if (pbj->id == 0)
		{
			//G_10=T_10*G_00
			//G0.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->calTij(t, y + ydimcum[I]);
			G0.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->Tij;
		}
		else
		{
			J = pbj->id - 1;
			int Jrow = qdimcum[J] - (Body::NC - 6)*J;
			int Jrdim = qdimvec[J] - (Body::NC - 6);
			//G_i0=T_ij*G_j0
			/*G0.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->calTij(t, y + ydimcum[I]) *
				G0.block(Jrow, bcol, Jrdim, cdim);*/
			G0.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->Tij *
				G0.block(Jrow, bcol, Jrdim, cdim);
		}
	}
	return true;
}

bool TreeSystem::calG(double t, double* y, OUT MatrixXd& G)
{
	Body* pbi = nullptr;
	Body* pbj = nullptr;
	int id;
	int brow, bcol, rdim, cdim;
	for (unsigned int i = 0; i < nb; ++i)
	{
		pbi = bjvec[i]->bptr;
		pbj = bjvec[i]->bptr;
		while (true)
		{
			if (pbj->id == 0)                    //pbj is pointer to the base body
				break;
			if (pbi->id == pbj->id)
			{
				id = pbi->id - 1;
				brow = qdimcum[id]- (Body::NC - 6) * id;
				bcol = ydimcum[id];
				rdim = qdimvec[id] - (Body::NC - 6);
				cdim = ydimvec[id];
				//G_ii=U_i
				G.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->Ui;
				//G.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->calUi(t, y + ydimcum[id]);
			}
			else
			{
				int I = pbi->id - 1;
				int K = pbj->id - 1;
				int J = bjvec[i]->jptr->Bj->id - 1;
				brow = qdimcum[I] - (Body::NC - 6) * I;
				bcol = ydimcum[K];
				rdim = qdimvec[I] - (Body::NC - 6);
				cdim = ydimvec[K];
				int Jrow = qdimcum[J] - (Body::NC - 6) * J;
				int Jrdim = qdimvec[J] - (Body::NC - 6);
				//G_ik=T_ij*G_jk
				/*G.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->calTij(t, y + ydimcum[I]) *
					G.block(Jrow, bcol, Jrdim, cdim);*/
				G.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->Tij *
					G.block(Jrow, bcol, Jrdim, cdim);
			}
			id = pbj->id - 1;
			pbj = bjvec[id]->jptr->Bj;
		}
	}
	return true;
}

bool TreeSystem::calM(double t, double* y, OUT MatrixXd& M)
{
	for (unsigned int i = 0; i < nb; ++i)
	{
		bjvec[i]->bptr->calMass(M, qdimcum[i] - i * (Body::NC - 6));
	}
	return true;
}

bool TreeSystem::calg(double t, double* y, double* dy, OUT MatrixXd& g)
{
	Body* pbi = nullptr;
	Body* pbj = nullptr;
	int id;
	int brow, bcol, rdim, cdim;
	for (unsigned int i = 0; i < nb; ++i)
	{
		pbi = bjvec[i]->bptr;
		pbj = bjvec[i]->bptr;
		while (true)
		{
			if (pbj->id == 0)                    //pbj is pointer to the base body
				break;
			if (pbi->id == pbj->id)
			{
				id = pbi->id - 1;
				brow = qdimcum[id] - (Body::NC - 6) * id;
				bcol = id;
				rdim = qdimvec[id] - (Body::NC - 6);
				cdim = 1;
				//g_ii=betai
				//g.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->calbetai(t, y + ydimcum[id],dy+ydimcum[id]);
				g.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->betai;
			}
			else
			{
				int I = pbi->id - 1;
				int K = pbj->id - 1;
				int J = bjvec[i]->jptr->Bj->id - 1;
				brow = qdimcum[I] - (Body::NC - 6) * I;
				bcol = K;
				rdim = qdimvec[I] - (Body::NC - 6);
				cdim = 1;
				int Jrow = qdimcum[J] - (Body::NC - 6) * J;
				int Jrdim = qdimvec[J] - (Body::NC - 6);
				//g_ik=T_ij*g_jk
				g.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->Tij *
					g.block(Jrow, bcol, Jrdim, cdim);
				/*g.block(brow, bcol, rdim, cdim) = bjvec[i]->jptr->calTij(t, y + ydimcum[I]) *
					g.block(Jrow, bcol, Jrdim, cdim);*/
			}
			id = pbj->id - 1;
			pbj = bjvec[id]->jptr->Bj;
		}
	}
	return true;
}

bool TreeSystem::calf(double t, double* y, double* dy, OUT VectorXd& f)
{
	for (unsigned int i = 0; i < nb; ++i)
	{
		f.segment(qdimcum[i] - (Body::NC - 6) * i, qdimvec[i] - (Body::NC - 6)) = bjvec[i]->bptr->inertiaForce();
	}
	return true;
}

bool TreeSystem::update(double t, double* y, double* dy)
{
	//update q dq Tij Ui betai
	if (AUX::isEqual(yprev.data(), y, yprev.rows()) && AUX::isEqual(dyprev.data(), dy, dyprev.rows()))
		return false;                        //already updated!
	if (!AUX::isEqual(yprev.data(), y, yprev.rows()))
	{
		for (int i = 0; i < yprev.rows(); ++i)
		{
			yprev[i] = y[i];
		}
		calytoq(t, y);
		static_cast<BaseBody*>(bjroot->bptr)->update(t);
		for (unsigned int i = 0; i < nb; ++i)
		{
			bjvec[i]->jptr->calTij(t, y + ydimcum[i]);
			bjvec[i]->jptr->calUi(t, y + ydimcum[i]);
		}
	}
	if (!AUX::isEqual(dyprev.data(), dy, dyprev.rows()))
	{
		for (int i = 0; i < dyprev.rows(); ++i)
		{
			dyprev[i] = dy[i];
		}
		/*when update dq, Tij and Ui have already been calculated.*/
		caldytodq(t, dy);
		double* yi = y;
		double* dyi = dy;
		for (unsigned int i = 0; i < nb; ++i)
		{
			bjvec[i]->jptr->calbetai(t, yi, dyi);
			yi += ydimvec[i];
		}
	}
	return true;
}

VectorXd TreeSystem::rootacc(double t)
{
	return static_cast<BaseBody*>(bjroot->bptr)->acceleration(t);
}

unsigned int TreeSystem::DOF() const
{
	return ydimcum[nb - 1] + ydimvec[nb - 1];
}

unsigned int TreeSystem::NC() const
{
	return qdimcum[nb - 1] + qdimvec[nb - 1] - nb * (Body::NC - 6);
}

TreeSystem::~TreeSystem()
{

}

bool TreeSystem::BJpair::operator==(BJpair& other)
{
	return (bptr == other.bptr) && (jptr == other.jptr);
}

bool TreeSystem::BJpair::operator<(BJpair& other)
{
	return (*bptr) < (*other.bptr);
}

bool TreeSystem::BJpair::operator>(BJpair& other)
{
	return (*bptr) > (*other.bptr);
}

MBSystem::MBSystem()
{
	if (peq != nullptr)
		delete peq;
	peq = new Equation(this);
	if (psolver != nullptr)
		delete psolver;
	psolver = new Solver(peq);
}

MBSystem::~MBSystem()
{
	if (!mbtree)
		delete mbtree;
	if (!psolver)
		delete psolver;

}

bool MBSystem::add(JointBase* j)
{
	if (std::find(jointvec.begin(), jointvec.end(), j) == jointvec.end())
	{
		jointvec.push_back(j);
		return true;
	}
	else
		return false;
}

bool MBSystem::del(JointBase* j)
{
	auto it = std::find(jointvec.begin(), jointvec.end(), j);
	if (it == jointvec.end())
		return false;
	jointvec.erase(it);
	return true;
}

unsigned int MBSystem::DOF() const
{
	return dof;
}

bool MBSystem::sety0(const VectorXd& _y0)
{
	y0 = _y0;
	return true;
}

bool MBSystem::setTimeInterval(double ti, double te, int N)
{
	return psolver->setTimeInterval(ti, te, N);
}

bool MBSystem::setTolerance(double r, double a)
{
	return psolver->setTolerance(r, a);
}

bool MBSystem::calculate()
{
	if (!initialize())
		return false;
	if (!psolver->calculate())
		return false;
	return true;
}
/* Save the data:
timespan
generalized coordinates velocity acceleration
absolute coordinates velocity acceleration*/
bool MBSystem::SaveAs(string fname, bool isbinary)
{
	ofstream fout;
	if (isbinary)
	{
		fname.append(".dat");
		fout.open(fname, ios::binary | ios::out);
	}
	else
	{
		fname.append(".txt");
		fout.open(fname, ios::out);
	}
	if (!fout)
	{
		cerr << "Open output file error!" << endl;
		return false;
	}
	unsigned int NT = psolver->tspan.size();
	
	if (isbinary)
	{
		fout.write(reinterpret_cast<char*>(&NT), sizeof(unsigned int));
		fout.write(reinterpret_cast<char*>(&dof), sizeof(unsigned int));
		fout.write(reinterpret_cast<char*>(&nc), sizeof(unsigned int));
		for (auto& item : psolver->tspan)
			fout.write(reinterpret_cast<char*>(&item), sizeof(double));
		for (auto& item : psolver->Y)
			fout.write(reinterpret_cast<char*>(&item), dof * sizeof(double));
		for (auto& item : psolver->DY)
			fout.write(reinterpret_cast<char*>(&item), dof * sizeof(double));
		for (auto& item : psolver->DY)
			fout.write(reinterpret_cast<char*>(&item + sizeof(double) * dof), dof * sizeof(double));
	}
	else
	{
		fout << NT << endl;
		fout << dof << endl;
		fout << nc << endl;
		for (auto& item : psolver->tspan)
			fout << item << endl;
		for (auto& item : psolver->Y)
			fout << item.head(dof).transpose() << endl;
		for (auto& item : psolver->DY)
			fout << item.head(dof).transpose() << endl;
		for (auto& item : psolver->DY)
			fout << item.tail(dof).transpose() << endl;
	}
	fout.close();
	return true;
}

Solver::Solver(Equation* p)
{
	pe = p;
	Atol = 1e-3;
	Rtol = 1e-4;
	t_ini = 0.;
	t_end = 1.;
	Nstep = 50;
}

Solver::~Solver()
{

}

bool Solver::setTimeInterval(double ti, double te, int N)
{
	t_ini = ti;
	t_end = te;
	Nstep = N;
	return true;
}

bool Solver::setTolerance(double r, double a)
{
	Rtol = r;
	Atol = a;
	return true;
}

bool Solver::calculate()
{
	double dt = (t_end - t_ini) / Nstep;
	VectorXd y0 = pe->initialvalue();
	VectorXd y = y0, dy;
	dy.resize(pe->DOF(), 0);
	dy.head(pe->DOF() / 2) = y.tail(pe->DOF() / 2);
	double t = t_ini;
	while (t<t_end)
	{
		tspan.push_back(t);
		Y.push_back(y);
		DY.push_back(dy);

		dy = pe->Left(t, y).partialPivLu().solve(pe->Right(t, y));
		y += dy * dt;
		t += dt;
	}
	return true;
}

Equation::Equation(MBSystem* p)
{

	pmbs = p;
	L.setZero(2 * pmbs->dof, 2 * pmbs->dof);
	R.setZero(2 * pmbs->dof);
	L.block(0, 0, pmbs->dof, pmbs->dof).setIdentity();
}

Equation::~Equation()
{

}

MatrixXd& Equation::Left(double t, VectorXd& y)
{
	// update first
	auto& n = pmbs->dof;
	pmbs->update(t, y.data(), y.data() + n);
	L.block(n, n, n, n) = pmbs->calZ(t, y.data());
	return L;
}

VectorXd& Equation::Right(double t, VectorXd& y)
{
	// update first
	pmbs->update(t, y.data(), y.data() + pmbs->dof);
	auto& n = pmbs->dof;
	R.segment(0, n) = y.segment(n, n);
	R.segment(n, n) = pmbs->calz(t, y.data(), y.data() + n) + pmbs->calfey(t, y.data(), y.data() + n);
	return R;
}

VectorXd Equation::initialvalue() const
{
	return pmbs->y0;
}

unsigned int Equation::DOF() const
{
	return pmbs->dof * 2;
}

MatrixXd& MBSystem::calG0(double t, double* y)
{
	// TODO: insert return statement here
	mbtree->calG0(t, y, G0);
	return G0;
}

MatrixXd& MBSystem::calG(double t, double* y)
{
	// TODO: insert return statement here
	mbtree->calG(t, y, G);
	return G;
}

MatrixXd& MBSystem::calg(double t, double* y, double* dy)
{
	// TODO: insert return statement here
	mbtree->calg(t, y, dy, g);
	return g;
}

MatrixXd& MBSystem::calM(double t, double* y)
{
	// TODO: insert return statement here
	mbtree->calM(t, y, M);
	return M;
}

VectorXd& MBSystem::calf(double t, double* y, double* dy)
{
	// TODO: insert return statement here
	mbtree->calf(t, y, dy, f);
	return f;
}

bool MBSystem::update(double t, double* y, double* dy)
{
	/*call mbtree to update first*/
	if (!mbtree->update(t, y, dy))
		return false;
	/*this function need to update some basic variable
	G g M f G0*/
	calM(t, y);
	calG0(t, y);
	calG(t, y);
	calg(t, y, dy);
	calf(t, y, dy);
	
	return true;
}


bool MBSystem::initialize()
{
	/*initialize for calculation*/
	if (mbtree != nullptr)
		delete mbtree;
	mbtree = new TreeSystem(jointvec);
	dof = mbtree->DOF();
	nc = mbtree->NC();
	//allocate space
	y0.resize(dof);
	G.resize(nc, dof);
	G0.resize(nc, 6+mbtree->bjroot->bptr->nMode());
	M.resize(nc, nc);
	g.resize(nc, mbtree->nb);
	Z.resize(dof, dof);
	z.resize(dof);
	fey.resize(dof);
	f.resize(nc);
	//set zeros
	y0.setZero();
	G.setZero();
	G0.setZero();
	M.setZero();
	g.setZero();
	Z.setZero();
	z.setZero();
	fey.setZero();
	f.setZero();
	return true;
}

MatrixXd& MBSystem::calZ(double t, double* y)
{
	// calculate generalized Mass matrix Z
	Z = G.transpose() * M * G;
	return Z;
}

VectorXd& MBSystem::calz(double t, double* y, double* dy)
{
	// TODO: insert return statement here
	z = G.transpose() * (f - M * g.rowwise().sum() - M * G0 * mbtree->rootacc(t));
	return z;
}

VectorXd& MBSystem::calfey(double t, double* y, double* dy)
{
	// TODO: insert return statement here
	return fey;
}