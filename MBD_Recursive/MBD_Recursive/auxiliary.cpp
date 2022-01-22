#include "auxiliary.h"

extern double AUX::Tol = 1e-10;

/*check once*/
bool AUX::A(IN double* q, OUT Matrix3d& M)
{
	M(0, 0) = 2 * (q[0] * q[0] + q[1] * q[1]) - 1;
	M(0, 1) = 2 * (q[1] * q[2] - q[0] * q[3]);
	M(0, 2) = 2 * (q[1] * q[3] + q[0] * q[2]);
	M(1, 0) = 2 * (q[1] * q[2] + q[0] * q[3]);
	M(1, 1) = 2 * (q[0] * q[0] + q[2] * q[2]) - 1;
	M(1, 2) = 2 * (q[2] * q[3] - q[0] * q[1]);
	M(2, 0) = 2 * (q[1] * q[3] - q[0] * q[2]);
	M(2, 1) = 2 * (q[2] * q[3] + q[0] * q[1]);
	M(2, 2) = 2 * (q[0] * q[0] + q[3] * q[3]) - 1;
	return true;
}
bool AUX::EQtoA(IN double* q, OUT Matrix3d& M)
{
	return A(q, M);
}
bool AUX::EAtoA(IN double* q, OUT Matrix3d& M)
{
	double s1 = sin(q[0]), c1 = cos(q[0]);
	double s2 = sin(q[1]), c2 = cos(q[1]);
	double s3 = sin(q[2]), c3 = cos(q[2]);
	M(0, 0) = c1 * c3 - s1 * c2 * s3;
	M(0, 1) = -c1 * s3 - s1 * c2 * c3;
	M(0, 2) = s1 * s2;
	M(1, 0) = s1 * c3 - c1 * c2 * s3;
	M(1, 1) = -s1 * s3 - c1 * c2 * c3;
	M(1, 2) = -c1 * s2;
	M(2, 0) = s2 * s3;
	M(2, 1) = s2 * c3;
	M(2, 2) = c2;
	return true;
}
bool AUX::CAtoA(IN double* q, OUT Matrix3d& M)
{
	double s1 = sin(q[0]), c1 = cos(q[0]);
	double s2 = sin(q[1]), c2 = cos(q[1]);
	double s3 = sin(q[2]), c3 = cos(q[2]);
	M(0, 0) = c2 * c3;
	M(0, 1) = -c2 * s3;
	M(0, 2) = s2;
	M(1, 0) = s1 * s2 * c3 + c1 * s3;
	M(1, 1) = -s1 * s2 * s3 + c1 * c3;
	M(1, 2) = -s1 * c2;
	M(2, 0) = -c1 * s2 * c3 + s1 * s3;
	M(2, 1) = c1 * s2 * s3 + s1 * c3;
	M(2, 2) = c1 * c2;
	return true;
}
/*check once*/
bool AUX::R(IN double* q, OUT Matrix<double, 3, 4>& M)
{
	M(0, 0) = -q[1];
	M(0, 1) = q[0];
	M(0, 2) = -q[3];
	M(0, 3) = q[2];
	M(1, 0) = -q[2];
	M(1, 1) = q[3];
	M(1, 2) = q[0];
	M(1, 3) = -q[1];
	M(2, 0) = -q[3];
	M(2, 1) = -q[2];
	M(2, 2) = q[1];
	M(2, 3) = q[0];
	return true;
}
/*check once*/
bool AUX::tilde(IN Vector3d& v, OUT Matrix3d& M)
{
	M(0, 0) = M(1, 1) = M(2, 2) = 0;
	M(0, 1) = -v(2);
	M(0, 2) = v(1);
	M(1, 2) = -v(0);
	M(1, 0) = -M(0, 1);
	M(2, 0) = -M(0, 2);
	M(2, 1) = -M(1, 2);
	return true;
}
/*check once*/
bool AUX::tilde(IN double* v, OUT Matrix3d& M)
{
	M(0, 0) = M(1, 1) = M(2, 2) = 0;
	M(0, 1) = -v[2];
	M(0, 2) = v[1];
	M(1, 2) = -v[0];
	M(1, 0) = -M(0, 1);
	M(2, 0) = -M(0, 2);
	M(2, 1) = -M(1, 2);
	return true;
}

/*check once*/
bool AUX::isEqual(const double* v, const double* w, const unsigned int& len)
{
	double e = 0.;
	for (unsigned int i = 0; i < len; ++i)
	{
		e += abs(v[i] - w[i]);
	}
	return (e < AUX::Tol);
}

bool AUX::AtoEQ(IN const Matrix3d& A, OUT double* q)
{
	q[0] = sqrt(1 + A.trace()) / 2;
	if (abs(q[0]) > AUX::Tol)
	{//lambda_0!=0
		q[3] = (A(1, 0) - A(0, 1)) / (4. * q[0]);
		q[1] = (A(2, 1) - A(1, 2)) / (4. * q[0]);
		q[2] = (A(0, 2) - A(2, 0)) / (4. * q[0]);
	}
	else
	{//lambda_0==0
		int m, n, l;
		A.diagonal().maxCoeff(&m);
		q[m + 1] = sqrt((1 + A(m, m)) / 2);
		n = (m + 1) % 3;
		l = (n + 1) % 3;
		q[n + 1] = (A(m, n) + A(n, m)) / (4.0 * q[m + 1]);
		q[l + 1] = (A(l, m) + A(m, l)) / (4.0 * q[m + 1]);
	}
	return true;
}
/*check once*/
bool AUX::AtoEA(IN const Matrix3d& A, OUT double* q)
{
	double s2 = sqrt(1 - A(2, 2) * A(2, 2));
	q[1] = acos(A(2, 2));
	q[0] = atan2(A(0, 2) / s2, -A(1, 2) / s2);
	q[2] = atan2(A(2, 0) / s2, A(2, 1) / s2);
	return true;
}
/*check once*/
bool AUX::AtoCA(IN const Matrix3d& A, OUT double* q)
{
	double c2 = sqrt(1 - A(0, 2) * A(0, 2));
	q[1] = asin(A(0, 2));
	q[2] = atan2(-A(0, 1) / c2, A(0, 0) / c2);
	q[0] = atan2(-A(1, 2) / c2, A(2, 2) / c2);
	return true;
}

bool AUX::Kr(IN double* q, IN RCORDS rtype, OUT MatR3CX& K)
{
	if (rtype == RCORDS::EULERANGLE)
	{
		K(0, 0) = K(1, 0) = K(2, 1) = 0;
		K(2, 0) = 1;
		K(0, 1) = cos(q[2]);
		K(1, 1) = sin(q[2]);
		K(0, 2) = sin(q[1]) * sin(q[0]);
		K(1, 2) = -sin(q[1]) * cos(q[0]);
		K(2, 2) = cos(q[1]);
	}
	else if (rtype == RCORDS::CARDANANGLE)
	{
		K(0, 0) = 1;
		K(0, 1) = K(1, 0) = K(2, 0) = 0;
		K(0, 2) = sin(q[1]);
		K(1, 1) = cos(q[0]);
		K(2, 1) = sin(q[0]);
		K(1, 2) = -cos(q[1]) * sin(q[0]);
		K(2, 2) = cos(q[1]) * cos(q[0]);
	}
	else if (rtype == RCORDS::EULERQUATERNION)
	{
		K(0, 0) = -q[1];
		K(0, 1) = q[0];
		K(0, 2) = -q[3];
		K(0, 3) = q[2];
		K(1, 0) = -q[2];
		K(1, 1) = q[3];
		K(1, 2) = q[0];
		K(1, 3) = -q[1];
		K(2, 0) = -q[3];
		K(2, 1) = -q[2];
		K(2, 2) = q[1];
		K(2, 3) = q[0];
		K *= 2;
	}
	return true;
}
