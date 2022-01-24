// MBD_Recursive.cpp : This file contains the 'main' function. Program execution begins and ends there.
// This is a Multibody Dynamics program using forward recursive formulation methods.
// 

#include <iostream>
#include <Dense>
#include "Body.h"
#include "auxiliary.h"
#include <ctime>
using namespace std;
using namespace Eigen;
int main()
{
	MatrixXd m(1000, 1000);
	MatrixXd n(1000, 1000);
	for (int i = 0; i < 1000; ++i)
		for (int j = 0; j < 1000; ++j)
		{
			m(i, j) = (i + j) % 2;
			n(i, j) = (i + j) % 2;
		}
	int N = 10;
	clock_t s;
	s = clock();
	for (int i = 0; i < N; ++i)
		m += n;
	cout << clock() - s << endl;
	s = clock();
	for (int i = 0; i < N; ++i)
		m += n;
	cout << clock() - s << endl;
	s = clock();
	for (int i = 0; i < N; ++i)
		m.noalias() += n;
	cout << clock() - s << endl;
    return 0;
}


