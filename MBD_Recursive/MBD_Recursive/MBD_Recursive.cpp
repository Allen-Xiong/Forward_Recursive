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
    MBFileParser mbparser;
	try
	{
		mbparser.Read("pendulum.json");
		mbparser.Simulate();
		mbparser.SaveDataAs("pendulumdata", false);
		mbparser.Write("pendulum.json");
	}
	catch (const MBException& mbexp)
	{
		cout << mbexp.what() << endl;
	}
	catch (const std::exception& stdexp)
	{
		cout << stdexp.what() << endl;
	}

    return 0;
}


