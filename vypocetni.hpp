#ifndef VYPOCETNI_HPP
#define VYPOCETNI_HPP
#include<vector>

using namespace std;

extern int N;

extern double z[6], d[5], K[5];
extern double Q[2], r[2], R[2], s[2];
extern double L, T, H;


void gradh(double x, double y, double *dhdx, double *dhdy);
void v(double x, double y, double k, double *vx, double *vy);

void simple_track_point(double x0, double y0, double krok, vector<double> *X, vector<double> *Y);
void track_point(double x0, double y0, double z0, double krok, vector<double> *X, vector<double> *Y, double *T);
double track_point(double x0, double y0, double z0, double krok, vector<double> *X, vector<double> *Y);
double hydraulic_head(double x, double y);
void testovaci_input();
double wellYield(int idx);
double wellDrawdown(int idx);
#endif
