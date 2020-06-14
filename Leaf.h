#pragma once

#include "FindMin.h"

using namespace std;

// This functions contains functions to calculate leaf
// level photosynthesis

double GetAmax(double sn, double sc);

// Collatz 1992 equation A12
// temperature in celcius
double Q10Func(double K25, double Q10, double T);

// Collatz 1991 equation A3
double FGammastar(double Oi, double tau);

// eq a2 collatz 1991
double FJeC3(double a, double alpha, double Q0, double Ci, double gammastar);

double FJeC4(double Q0, double a, double alpharf);

//eq a5 collatz 1991
double FJcC3(double Vm, double Ci, double gammastar, double Kc, double Ko, double Oi);

// Collatz 1992 eq5
double FJcC4(double Vm);

// eq A6 Collatz 1991
double FJsC3(double Vm);

// collatz 1992 eq4
double FJpC4(double Ci, double P, double k);

// compute Vmax for C3 and C4 from Nt
// with a linear model
// follows Arora
// double FVmC3( double Nt );

// double FVmC4( double Nt );

// used for both c3 and c4.
// follows woodward Ci is 0.7 of atm CO2 for C3
// for C4 we assume Ci is 5 times atm CO2 ???
// based on Collatz et al. 1991 equations but woodward concept
// for calculating Vmax
// double FVm( double An, double Rd, double Ca, double CiFrac, double gammastar, 
// 		   double Kc, double Ko, double Oi, double T );
double FVm(double An, double T);

// double FVmWithCi( double An, double Rd, double Ci, double gammastar, 
// 				 double Kc, double Ko, double Oi, double T );

// day respiration
// rate should be 0.015 for C3 and 0.025 for C4
// many sources, not sure who was first, see Arora
double RmL(double Vm, double rate);

// takes gs in mol/m2/s gives Ad in micromols/m2/s
// re-arrangement for equation 6B from collatz 1992
double CollatzAd(double gs, double Ci, double Ca, double P);

// collatz 1991 equation 1
// hs rh at leaf surface
// gives conductance in mol/m2/s
double CollatzBallBerryLeaf(double A, double Ca, double P, double m, double b, double hs, double gb);

double GetAd(double m, double b, double Ab, double Rh, double P, double Ca, double Ci, double gb);

// gives aerodynamic conductance in mol/m2/s
// based on Jones 1992 empirical function
// coeff 0.271 assumes standard atm press and 25 deg C
double GetgbLEAF(double u, double cld);

double FindC3Ci(double Ci, double *aa);

double FindC4Ci(double Ci, double *aa);

//var lVm,lO2conc,ltau,lKc,lKo,laC4,lalpha,lQp,lCit,lP,lT,lws,lRh:real
double FindC3Min(double Vm, double Oi, double gammastar, double Kc, double Ko, double aC3, double mtree, double btree, double GammaC3, double alpha, double Q0,
		double P, double Ca, double T, double Rh, double gb);

//var lVm,lO2conc,ltau,lKc,lKo,laC4,lalpha,lQp,lCit,lP,lT,lws,lRh:real
double FindC4Min(double Vm, double aC4, double k, double mgrass, double bgrass, double GammaC4, double alpharf, double Q0, double P, double Ca, double T,
		double Rh, double gb);

// compute A0 and R for C3-plants
void GetC3A(double T, double sn, double sc, double Ca, double Oi, double Q0, double GammaC3, double aC3, double alpha, double Rh, double P, double wind,
		double CLDtree, double mtree, double btree, double tauC3K25, double tauC3Q10, double KoC3K25, double KoC3Q10, double KcC3K25, double KcC3Q10, double *A,
		double *RmLv);

// compute A0 and R for C4-plants
void GetC4A(double T, double sn, double sc, double Ca, double k, double Q0, double GammaC4, double aC4, double alpharf, double Rh, double P, double wind,
		double CLDgrass, double mgrass, double bgrass, double *A, double *RmLv);

//------------------------------------------------------------------------------------------------

double GetAmax(double sn, double sc)
{
	double am;

	am = 50. * pow(0.999927, sc);

	if (sn < 600.) am = am * 0.00166 * sn;

// 	return am*2.;
	return am;
}

double Q10Func(double K25, double Q10, double T)
{
	return K25 * pow(Q10, (T - 25.) / 10.);
}

double FGammastar(double Oi, double tau)
{
	return Oi / (2. * tau);
}

double FJeC3(double a, double alpha, double Q0, double Ci, double gammastar)
{
	return MyMax(0, a * alpha * Q0 * ((Ci - gammastar) / (Ci + 2. * gammastar)));
}

double FJeC4(double Q0, double a, double alpharf)
{
	return a * alpharf * Q0;
}

double FJcC3(double Vm, double Ci, double gammastar, double Kc, double Ko, double Oi)
{
	return Vm * (Ci - gammastar) / (Ci + Kc * (1. + Oi / Ko));
}

double FJcC4(double Vm)
{
	return Vm;
}

double FJsC3(double Vm)
{
	return Vm / 2.;
}

double FJpC4(double Ci, double P, double k)
{
	return k * Ci / P;
}

double FVm(double An, double T)
{
	double sol;

	sol = An * pow(2., 0.1 * (T - 25.)) / ((1. + exp(0.3 * (13. - T))) * (1. + exp(0.3 * (T - 36.))));

	return sol;
}

double RmL(double Vm, double rate)
{
	return rate * Vm;
}

double CollatzAd(double gs, double Ci, double Ca, double P)
{
	return (Ca - Ci) * gs / (1.6 * P);
}

double CollatzBallBerryLeaf(double A, double Ca, double P, double m, double b, double hs, double gb)
{
	double cs;

	cs = Ca - 1.4 * A * P / gb; 		//this line follows Arora

	double ret = MyMax(0, m * A * (hs / 100.) * P / cs + b); 			//ie for negative An we produce zero gs

	return ret;
}

double GetAd(double m, double b, double Ab, double Rh, double P, double Ca, double Ci, double gb)
{
	double gs;
	gs = CollatzBallBerryLeaf(Ab, Ca, P, m, b, Rh, gb);
	return CollatzAd(gs, Ci, Ca, P);
}

double GetAd_end(double m, double b, double Ab, double Rh, double P, double Ca, double Ci, double gb, double vt)
{
	double gs;
	gs = CollatzBallBerryLeaf(Ab, Ca, P, m, b, Rh, gb);

	if (vt == 3.) gs_C3_global += gs;
	if (vt == 4.) gs_C4_global += gs;
	return CollatzAd(gs, Ci, Ca, P);
}

double GetgbLEAF(double u, double cld)
{
	return 0.271 * 1e6 * pow(u / cld, 0.5);
}

double FindC3Ci(double Ci, double *aa)
{
	double Je;
	double Jc;
	double Js;
	double Ab;
	double Ad;
	double Jmin;
	double ret;

	Ab = 1e10;
	Ad = 1e5;

	if (Ci > 2.)
	{
		Je = FJeC3(aa[6], aa[7], aa[8], Ci, aa[3]);
		Jc = FJcC3(aa[1], Ci, aa[3], aa[4], aa[5], aa[2]);
		Js = FJsC3(aa[1]);
		Jmin = MyMin(MyMin(Je, Jc), Js);
		Ab = MyMax(0, Jmin - RmL(aa[1], aa[14])); //A_biochemical
		Ad = GetAd(aa[12], aa[13], Ab, aa[11], aa[16], aa[17], Ci, aa[15]); //A_diffusion
		ret = fabs(Ab - Ad);
	}
	else
		ret = 10000.;

	return ret;
}

double FindC4Ci(double Ci, double *aa)
{
	double Jc;
	double Je;
	double Jp;
	double Ab;
	double Ad;
	double Jmin;

	double ret = 0;

	Ab = 1e10;
	Ad = 1e5;

	if (Ci > .5)
	{
		Jc = FJcC4(aa[1]);
		Je = FJeC4(aa[6], aa[5], aa[3]);
		Jp = FJpC4(Ci, aa[7], aa[2]);
		Jmin = MyMin(MyMin(Je, Jc), Jp);
		Ab = MyMax(0, Jmin - RmL(aa[1], aa[4])); //A_biochemical
		Ad = GetAd(aa[10], aa[11], Ab, aa[9], aa[13], aa[14], Ci, aa[12]);  //A_diffusion
		ret = fabs(Ab - Ad);
	}
	else
		ret = 10000.;

	return ret;
}

double FindC3Min(double Vm, double Oi, double gammastar, double Kc, double Ko, double aC3, double mtree, double btree, double GammaC3, double alpha, double Q0,
		double P, double Ca, double T, double Rh, double gb)
{
	double ax;
	double bx;
	double cx;
	double fa;
	double fb;
	double fc;
	double xmin;
	double aa[18];

	aa[1] = Vm;
	aa[2] = Oi;
	aa[3] = gammastar;
	aa[4] = Kc;
	aa[5] = Ko;
	aa[6] = aC3;
	aa[7] = alpha;
	aa[8] = Q0;
	aa[9] = P;
	aa[10] = T;
	aa[11] = Rh;
	aa[12] = mtree;
	aa[13] = btree;
	aa[14] = GammaC3;
	aa[15] = gb;
	aa[16] = P;
	aa[17] = Ca;

	ax = .2;
	bx = 1000.;
	cx = 10000.;
	fa = FindC3Ci(ax, aa);
	fb = FindC3Ci(bx, aa);
	fc = FindC3Ci(cx, aa);

	mnbrak(FindC3Ci, aa, &ax, &bx, &cx, &fa, &fb, &fc);

	golden(FindC3Ci, aa, ax, bx, cx, 0.0000000001, &xmin);

	return xmin;
}

double FindC4Min(double Vm, double aC4, double k, double mgrass, double bgrass, double GammaC4, double alpharf, double Q0, double P, double Ca, double T,
		double Rh, double gb)
{
	double ax;
	double bx;
	double cx;
	double fa;
	double fb;
	double fc;
	double xmin;
	double aa[18];

	aa[1] = Vm;
	aa[2] = k;
	aa[3] = alpharf;
	aa[4] = GammaC4;
	aa[5] = aC4;
	aa[6] = Q0;
	aa[7] = P;
	aa[8] = T;
	aa[9] = Rh;
	aa[10] = mgrass;
	aa[11] = bgrass;
	aa[12] = gb;
	aa[13] = P;
	aa[14] = Ca;

	ax = .2;
	bx = 1000.;
	cx = 10000.;
	fa = FindC4Ci(ax, aa);
	fb = FindC4Ci(bx, aa);
	fc = FindC4Ci(cx, aa);

	mnbrak(FindC4Ci, aa, &ax, &bx, &cx, &fa, &fb, &fc);

	golden(FindC4Ci, aa, ax, bx, cx, 0.0000000001, &xmin);

	return xmin;
}

void GetC3A(double T, double sn, double sc, double Ca, double Oi, double Q0, double GammaC3, double aC3, double alpha, double Rh, double P, double wind,
		double CLDtree, double mtree, double btree, double tauC3K25, double tauC3Q10, double KoC3K25, double KoC3Q10, double KcC3K25, double KcC3Q10, double *A,
		double *RmLv)
{
	double Amax;
	double tau;
	double Ko;
	double Kc;
	double Vm;
	double gammastar;
	double Ci;
	double Je;
	double Jc;
	double Js;
	double Jmin;
	double gb_molar;

	Ci = 5.;
	Amax = GetAmax(sn, sc);
	tau = Q10Func(tauC3K25, tauC3Q10, T);
	Ko = Q10Func(KoC3K25, KoC3Q10, T);
	Kc = Q10Func(KcC3K25, KcC3Q10, T);
	gammastar = FGammastar(Oi, tau);
	Vm = FVm(Amax, T);
	gb_molar = GetgbLEAF(wind, CLDtree);

	Ci = FindC3Min(Vm, Oi, gammastar, Kc, Ko, aC3, mtree, btree, GammaC3, alpha, Q0, P, Ca, T, Rh, gb_molar);

	Je = FJeC3(aC3, alpha, Q0, Ci, gammastar);
	Jc = FJcC3(Vm, Ci, gammastar, Kc, Ko, Oi);
	Js = FJsC3(Vm);
	Jmin = MyMin(MyMin(Je, Jc), Js);
	*A = Jmin;
	*RmLv = RmL(Vm, GammaC3);

	GetAd_end(mtree, btree, *A, Rh, P, Ca, Ci, gb_molar, 3.); //A_diffusion

	return;
}

void GetC4A(double T, double sn, double sc, double Ca, double k, double Q0, double GammaC4, double aC4, double alpharf, double Rh, double P, double wind,
		double CLDgrass, double mgrass, double bgrass, double *A, double *RmLv)
{
	double Amax;
	double Vm;
	double Ci;
	double Je;
	double Jc;
	double Jp;
	double Jmin;
	double gb_molar;

	Ci = 5.;
	Amax = GetAmax(sn, sc);
	Vm = FVm(Amax, T);

	Vm /= 2.3;   // Arora 2.0, Collatz 90/39=2.3

	gb_molar = GetgbLEAF(wind, CLDgrass);

	Ci = FindC4Min(Vm, aC4, k, mgrass, bgrass, GammaC4, alpharf, Q0, P, Ca, T, Rh, gb_molar);

	Jc = FJcC4(Vm);                     //Vm is temp dep
	Je = FJeC4(Q0, aC4, alpharf);         //temp effect following collatz
	Jp = FJpC4(Ci, P, Q10Func(k, 2., T));   //no temp limit

	Jmin = MyMin(MyMin(Je, Jc), Jp);
	*A = Jmin;
	*RmLv = RmL(Vm, GammaC4);

	GetAd_end(mgrass, bgrass, *A, Rh, P, Ca, Ci, gb_molar, 4.);

	return;
}

