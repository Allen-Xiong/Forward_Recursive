// MBD_Recursive.cpp : This file contains the 'main' function. Program execution begins and ends there.
// This is a Multibody Dynamics program using forward recursive formulation methods.
// 

#include <iostream>
#include <Dense>
#include "Body.h"
#include "auxiliary.h"

using namespace std;
using namespace Eigen;
int main()
{
    Matrix3d A;
    A << 0.433, -0.25, 0.866,
        0.808, 0.5335, -0.25,
        -0.3995, 0.808, 0.433;
    double* q = new double[4];
    double* p = q;
    AUX::AtoCA(A, q);
    for (int i = 0; i < 4; ++i)
        cout << q[i] << endl;
    if (p == q)
        cout << "true" << endl;
    delete[]q;
    string s = "Body";
    string k = "Body";
    if (s == "Body")
        cout << "true" << endl;
    return 0;
}


