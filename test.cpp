#include "vypocetni.hpp"
#include<fstream>
#include<iostream>
#include<cmath> // jen kvuli sin a cos

using namespace std;

void refresh();

int main()
{
	refresh();
	double x[100], y[100];

	// nalinkovat osy
	for(int i = 0; i < 100; i++) {
		x[i] = -25.0 + i;
		y[i] = -50.0 + i; }

// 	// test gradh - OK
// 	for(int i = 0; i < 100; i++)
// 		for(int j = 0; j < 100; j++)
// 		{
// 			double gx,gy;
// 			gradh(x[i], y[j], &gx, &gy);
// 			cout << x[i] << "\t" << y[j] << "\t" << gx << "\t" << gy << "\n";
// 		}
// 	cout << "\n";

// 	// test v - OK
	ofstream rpole;
	rpole.open("rychlostni_pole.dat", ios_base::out);
	for(int i = 0; i < 100; i++)
		for(int j = 0; j < 100; j++)
		{
			double vx, vy;
			v(x[i], y[j], 3e-5, &vx, &vy);
			rpole << x[i] << "\t" << y[j] << "\t" << vx << "\t" << vy << "\n";
		}
	rpole << "\n";
	rpole.close();

	// test track_point()
	double krok = .05; // m
	vector<double> X,Y;
	ofstream lajny;
	lajny.open("lajny.dat", ios_base::out);
	
	for(double fi = 0; fi < 2*M_PI; fi+= M_PI/6) // po 30 stupnich
	{
		double x0 = (.1 + r[0]) * cos(fi); // startovni bod kousek od vrtu 2
		double y0 = (.1 + r[0]) * sin(fi);
		
		double T = track_point(x0, y0, 10.0, krok, &X, &Y); // zhotovit lajnu
	
		for(int i = 0; i < X.size(); i++) // zapsat lajnu do souboru
			lajna << X[i] << "\t" << Y[i] << "\n";
	}
	
	lajna.close();
	
	cerr << "T = " << T << " s" << endl;

return 0; }



void refresh()
{
	N = 5;

		// udaje o studnach:
	Q[0] = .005; // m3/s
	Q[1] =-.001; // m3/s
	r[0] = r[1] = .125; //m
	R[0] = 160;
	R[1] = 135; //m
	s[0] = 2.5;
	s[1] =-1.5;

	L = 50.0; // m

	// udaje o vrstvach:
	z[0] =  0.0; // teren
	z[1] =  2.0;
	z[2] =  5.0;
	z[3] =  8.0;
	z[4] = 12.0;
	z[5] = 15.0; // podlozi

	K[0] = 3e-5; // m/s
	K[1] = 3e-3; // m/s
	K[2] = 3e-6; // m/s
	K[3] = 3e-5; // m/s
	K[4] = 3e-6; // m/s

	T = 0;
	for(int i = 0; i < N; i++)
	{
		d[i] = z[i+1] - z[i];
		T+= K[i]*d[i];
	}

	H = 4.0; // hloubka hpv
}