#include<cmath>
#include<iostream>
#include"vypocetni.hpp"

using namespace std;

int N;

double z[6], d[5], K[5];
double Q[2], r[2], R[2], s[2];
double L, T, H;

double getK(double Z);
int getLayer(double Z);
void gradG(double x, double y, double *dGdx, double *dGdy);
double IntegralOverOneLayer(double a, double b, double k, double h);
double Integral(double h);
// void solve_quadratic(double a, double b, double c, double *x1, double *x2);
double GirPot(double r1, double r2);

int sgn(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

// vypocitat cerpani ze snizeni
double wellYield(int idx)
{
    if((idx!=0) && (idx!=1))
    {
        cerr << "wellYield: bad parameter value!" << endl;
        return 0.0;
    }

    //pocitam cerpani z dane studny - ze snizeni v te studne a vydatnosti druhe studny
    int id2 = 1 - idx; // index druhe studny

    double ro[] = {fabs(idx*L - r[idx]), fabs(id2*L - r[idx])}; // vzdalenosti pro bod na plasti dane studny

    for(int k = 0; k < 2; k++)
        if(ro[k] > R[k]) ro[k] = R[k]; // zkratit vzdalenosti pro pripad, ze by depresni kuzel nektere studny nezasahoval az ke druhe

    // nechci hned prepsat Q-cko, proto NE "Q[idx] = ..."
    //Q[idx] = 1.0/log(r[idx]/R[idx]) * (2*M_PI * (Integral(H-s[idx]) - Integral(H)) - Q[id2]*log(r[id2]/R[id2]));
    return (2*M_PI * (Integral(H-s[idx]) - Integral(H)) - Q[id2]*log(r[id2]/R[id2])) / log(r[idx]/R[idx]);
}

// snizeni ve studnach, pokud jsou zadana Q[]
double wellDrawdown(int idx)
{
    if((idx!=0) && (idx!=1))
    {
        cerr << "wellDrawdown: bad parameter value!" << endl;
        return 0.0;
    }

    return (H - hydraulic_head(fabs(idx*L - r[idx]),0)); // mwahaha! ... :) vyuzivam toho, ze studna[0] ma vzdycky souradnice [0,0] a druha [L,0]
}

// pocitame hladinu v bode [x,y]
double hydraulic_head(double x, double y)
{
	double r1 = sqrt(x*x + y*y);
	double r2 = sqrt((x-L)*(x-L) + y*y);
	
    if(r1 < r[0])
        return hydraulic_head(r[0],0);

    if(r2 < r[1])
        return hydraulic_head(L+r[1],0);

	double I0 = Integral(0);
    double Gh = GirPot(r1, r2) - I0; // vyzaduje Q-cka
	double GH = Integral(H) - I0;
	
	double h_odhad = Gh/GH * H;
	double h_pom = H;
	
	double alfa = 0.6; // tlumic kroku
	
/*	cerr << "[ " << x << "; " << y << "] : " << endl;
	cerr << "GH = " << GH << endl;
	cerr << "Gh = " << Gh << endl;
	cerr << "h_odhad = " << h_odhad << endl;
	cerr << "Integral(0) = " << I0 << endl;
*/
	
	do {
		double G_odhad = Integral(h_odhad) - I0;
		h_pom = h_odhad;
		
        //cerr << "h_odhad = " << h_odhad << endl;
		
		h_odhad+= (Gh - G_odhad)/GH * H * alfa;
	} while(abs(h_odhad - h_pom)/alfa > 5e-4); // hyd. vyska s chybou .5 mm musi stacit kazdymu
	
	return h_odhad;
}

double GirPot(double r1, double r2)
{
    if((sgn(r1) < 1) || (sgn(r2) < 1))
    {
        cerr << "GirPot: zaporny argument!" << endl;
        return 0.0;
    }

	return (0.5/M_PI * (Q[0]*log(r1/R[0]) + Q[1]*log(r2/R[1])) + Integral(H));
}

// void solve_quadratic(double a, double b, double c, double *x1, double *x2)
// {
// 	double D = b*b - 4*a*c;
// 	
// 	if(D < 0)
// 	{
// 		cerr << "Kvadraticka rovnice nema koreny! :((" << endl;
// 		cerr << "a = " << a << endl;
// 		cerr << "b = " << b << endl;
// 		cerr << "c = " << c << endl;
// 		// *x1 = *x2 = 0.0;
// 		return;
// 	}
// 	
// 	*x1 = (-b + sqrt(D))/2.0/a;
// 	*x2 = (-b - sqrt(D))/2.0/a;
// }
//--------------------------------- H O T O V E -----------------------------------

double getK(double Z) // urci, jake K odpovida zadanemu z
{
	int idx = getLayer(Z);
	
	if(idx > -1)
		return K[idx];
	else
		return -1.0;
}

int getLayer(double Z)
{
	for(int i=0; i < N; i++)
		if( (z[i] <= Z) && (Z < z[i+1] )) //match
			return i;
	
	// jestli cyklus skoncil a my jsme dosli sem, ocividne se stala chyba - nesmyslne z
	return -1;
}


double track_point(double x0, double y0, double z0, double krok, vector<double> *X, vector<double> *Y)
{
	double T;
	
	track_point(x0, y0, z0, krok, X, Y, &T);
	
	return T;
}

void track_point(double x0, double y0, double z0, double krok, vector<double> *X, vector<double> *Y, double *T)
{
	X->clear();
	Y->clear();
	*T = 0;
	
	double x = x0;
	double y = y0;
	double K = getK(z0);
	double l[2];
	
	if (K < 0)
	{
		cerr << "Nepodarilo se najit odpovidajici K." << endl;
		return;
	}
	
	do {
		X->push_back(x);
		Y->push_back(y);
		
		double vx, vy;
		v(x, y, K, &vx, &vy); // slozky obj. hustoty toku
		double vel = sqrt(vx*vx+vy*vy);
		double dt = krok / vel; // skalovat vektor v na delku kroku
		double dx = vx * dt;
		double dy = vy * dt; // prevedeno na vektor delky "krok"
		
		x+= dx; // popojit o krok
		y+= dy;
		(*T)+= dt;
		
		l[0] = sqrt(x*x+y*y); // vzdalenost sledovane castice od prvni/druhe studny
		l[1] = sqrt((x-L)*(x-L)+y*y);
	
#ifdef DEBUG
		cout << "l[0] = " << l[0] << endl;
		cout << "l[1] = " << l[1] << endl;
#endif
		
	} while ( !( (l[0]<r[0]) || (l[1]<r[1]) || ( (l[0]>R[0]) && (l[0]>R[0]) )   ));
}

void simple_track_point(double x0, double y0, double krok, vector<double> *X, vector<double> *Y) // trajektorie konci, jakmile se dostaneme do vrtu
{
	X->clear();
	Y->clear();
	
	double x = x0;
	double y = y0;
	double dx, dy, koef;
	double l[2];
	
	do {
		X->push_back(x);
		Y->push_back(y);
		
		gradh(x, y, &dx, &dy); // pro anizotropii pouzit void v(...)!!!!!!!
		koef = krok / sqrt(dx*dx+dy*dy);
		dx*= koef;
		dy*= koef; // prevedeno na vektor delky "krok"
		
		x+= dx; // popojit o krok
		y+= dy;
		
		l[0] = sqrt(x*x+y*y); // vzdalenost sledovane castice od prvni/druhe studny
		l[1] = sqrt((x-L)*(x-L)+y*y);
	
#ifdef DEBUG
		cout << "l[0] = " << l[0] << endl;
		cout << "l[1] = " << l[1] << endl;
#endif
		
	} while ( !( (l[0]<r[0]) || (l[1]<r[1]) || ( (l[0]>R[0]) && (l[0]>R[0]) )   ));
}

double Integral(double h)
{
	bool konec = false;
	double Pot = 0.0;

	for(int i=0; i<N; i++)
	{
		double a = z[i];
		double b;

		if (h > z[i+1])
			b = z[i+1];
		else
		{
			b = h;
			konec = true;
		}

		Pot+= IntegralOverOneLayer(a,b,K[i],h);

		if (konec) break;
	}

	return Pot;
}

double IntegralOverOneLayer(double a, double b, double k, double h)
{
	return ( 0.5 * k * (b*b - a*a) - h * k * (b - a));
}

void gradG(double x, double y, double *dGdx, double *dGdy)
{
	double vzd2[2] = { (x*x+y*y), ((x-L)*(x-L)+y*y) }; // ctverce vzdalenosti od vrtu - hodi se i do vzorcu
	
	*dGdx = 0.0;
	*dGdy = 0.0;
	
	if(vzd2[0] < R[0]*R[0])
	{
		double pom = Q[0] / vzd2[0];
		*dGdx+= pom * x;
		*dGdy+= pom * y;
	}
	
	if(vzd2[1] < R[1]*R[1])
	{
		double pom = Q[1] / vzd2[1];
		*dGdx+= pom * (x-L);
		*dGdy+= pom * y;
	}
	
	double pom = 0.5/M_PI;
	*dGdx*= pom;
	*dGdy*= pom;
}

void gradh(double x, double y, double *dhdx, double *dhdy)
{
	gradG(x, y, dhdx, dhdy);
	*dhdx/= -T;
	*dhdy/= -T;
}

void v(double x, double y, double k, double *vx, double *vy)
{
	gradh(x,y,vx,vy);
	*vx*= -k; // tady by se zavedla anizotropie
	*vy*= -k;
}

void testovaci_input()
{
    N = 5;

        // udaje o studnach:
    Q[0] =-.005; // m3/s Q<0 cerpani, Q>0 zasakovani
    Q[1] = .003; // m3/s
    r[0] = r[1] = .125; //m
    R[0] = 40;
    R[1] = 20; //m
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

    // tenhle kus kodu pripravi z-souradnice vrstevnich rozhrani - obrati osu, polozi z=0 na podlozi
    double Z_base = z[N];

    for(int i = 0; i < N+1; i++)
    {
        z[i] = Z_base - z[i];
        //cerr << "z[" << i << "] = " << z[i] << endl;
    }

    for(int i = 0; i < (N+1)/2; i++)
    {
        double pom;
        pom = z[i];
        z[i] = z[N-i];
        z[N-i] = pom;
        //cerr << "z[" << i << "] = " << z[i] << endl;
    }

    H = Z_base - H;
}
