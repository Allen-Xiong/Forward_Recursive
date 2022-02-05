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
	char instr[200];
	cout << "Please Enter Input File Name(.json):" << endl;
	cin >> instr;
	string fname(instr);
	string outdatafile = fname;
	outdatafile.erase(outdatafile.find_last_of('.'), 5);
	clock_t start = clock();
    MBFileParser mbparser;
	try
	{
		mbparser.Read(fname);
		mbparser.Simulate();
		cout << "The computation costs " << (double)(clock() - start) / CLOCKS_PER_SEC << " s." << endl;
		mbparser.SaveDataAs(outdatafile, false);
		cout << "The Output Data File is " << outdatafile << ".txt" << endl;
	}
	catch (const MBException& mbexp)
	{
		cout << mbexp.what() << endl;
		cout << "Simulation Failed." << endl;
	}
	catch (const std::exception& stdexp)
	{
		cout << stdexp.what() << endl;
		cout << "Simulation Failed." << endl;
	}
	
    return 0;
}


