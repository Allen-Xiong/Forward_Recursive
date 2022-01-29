// MBD_Recursive.cpp : This file contains the 'main' function. Program execution begins and ends there.
// This is a Multibody Dynamics program using forward recursive formulation methods.
// 

#include <iostream>
#include "MultibodySystem.h"
#include <ctime>
using namespace std;
using namespace Eigen;
int main()
{
	clock_t start = clock();
    MBFileParser mbparser;
	try
	{
		mbparser.Read("solarpaneldeploy.json");
		mbparser.Simulate();
		mbparser.SaveDataAs("solarpaneldeploy", false);
	}
	catch (const MBException& mbexp)
	{
		cout << mbexp.what() << endl;
	}
	catch (const std::exception& stdexp)
	{
		cout << stdexp.what() << endl;
	}
	cout << "The computation costs " << (clock() - start) / CLOCKS_PER_SEC << " s." << endl;
    return 0;
}


