#pragma once

#include "MyMath.h"

using namespace std;

// Convert a double to a string
string doubleToString(double d)
{
	stringstream ss;
	ss << d;
	return (ss.str());
}

// gives gs in micro mol/m2/s from gs in m/s (derived from gNonMolar)
double gMolar(double T, double gs, double P);

// gives m/s from gs in micro mol/m2/s
// this comes from Arora submitted (following Sellers)
double gNonMolar(double T, double gs, double P);

// cox appendix 40
// double FQsumTree( double LAI );
// double FQsumGrass( double LAI );

// provides the depth of root given that it is shaped like a cylinder
double CylinderDepth(double Broot, double MinRootRad, double RootDensity, double Dmrd);

// gives aerodynamic conductance in m/s
// follows sellers 1996, eq. (14)
// follows Jones 1992
// 0.1681 = 0.41^2 = (von Karman constant)^2
// note need to fix this up - zwind-zd can be negative if vegetation gets very tall
double GetgbCANOPY(double Uz, double hbar);

// compute the wind at a certain height
double WindAtHeight(double height, double wind_ref);

// root and stem maintenance respiration following Arora Agr and For Met
// kg/day
// double MaintResp( double Biomass, double BetaN, double upsilon, double T );

double MaintRespFast(double Biomass, double BetaN, double upsilon, double T_fac);

double MaintRespTmpFac(double T);

//-----------------------------------------------------------------------------------------------

// woodward after shugart; height (m) from LAI
// double E25h( double Lai );

// compute available light of unshaded tree
// double FQsumShade( double LAIupper, double LAIlower );

// compute available light of a big tree (which shades an other one)
// double FQsumShadeBigPlant( double LAIbig, double LAIsmall );

// compute available light of a small tree (which is shaded)
// double FQsumShadeSmallPlant( double LAIbig, double LAIsmall );

// light competition index based on Fshade reasoning
// double FQi( double LAIcomp, double LAItarget );

// collatz 1992 equation 1
// hs rh at leaf surface
// gives conductance in micro mol/m2/s
// double CollatzBallBerryCanopy( double A, double resp, double Ca, double P, double m, double b, double hs, double gb, double LAI );

//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

double gMolar(double T, double gs, double P)
{
	return gs * 1e6 * 44.64286 * (273.16 / (T + 273.15)) * (1.013e5 / P);
}

//-----------------------------------------------------------------------------------------------

double gNonMolar(double T, double gs, double P)
{
	return gs * 1e-6 * 0.0224 * ((T + 273.15) / 273.16) * (P / 1.013e5);
}

//-----------------------------------------------------------------------------------------------

// double FQsumTree( double LAI )
// {
// 	return (1.-exp(-K_CAN_EXT_TREE*(LAI+0.4)))*K_CAN_EXT_TREE_INVERSE;     // NOTE: inserted 0.4 for callibration (Qsum=1 for LAI=1)
// }

// double FQsumGrass( double LAI )
// {
// 	return (1.-exp(-K_CAN_EXT_GRASS*(LAI+0.4)))*K_CAN_EXT_GRASS_INVERSE;     // NOTE: inserted 0.4 for callibration (Qsum=1 for LAI=1)
// }

//-----------------------------------------------------------------------------------------------

double CylinderDepth(double Broot, double MinRootRad, double RootDensity, double Dmrd)
{
	return MyMin(Dmrd, Broot / (RootDensity * M_PI * MinRootRad * MinRootRad));
// 	return Dmrd + Broot/(RootDensity*M_PI*MinRootRad*MinRootRad);
}

//-----------------------------------------------------------------------------------------------

double GetgbCANOPY(double Uz, double hbar)
{
	double Zd;
	double Z0;

	if (hbar > 11.) hbar = 11.; // necessary to avoid negative arguments in log

	Zd = ZD_CONST * hbar;
	Z0 = Z0_CONST * hbar;

	double log_tmp = 1. / log((REF_HEIGHT_Z - Zd) / Z0);

// 	cout << setw(15) << hbar << setw(15) << Uz << setw(15) << KARMAN_CONST_QUAD*Uz*log_tmp*log_tmp << setw(15) << log_tmp << endl;	                    

// 	return KARMAN_CONST_QUAD*Uz/pow( log((10.-Zd)/Z0),2. );
// 	return KARMAN_CONST_QUAD*Uz/pow( log((REF_HEIGHT_Z-Zd)/Z0),2. );
	return KARMAN_CONST_QUAD * Uz * log_tmp * log_tmp;
}

//-----------------------------------------------------------------------------------------------
//  NOTE: Old code
// double MaintResp( double Biomass, double BetaN, double upsilon, double T )
// {
// 	double tmp;
// 	
// 	T = MyMax(T,20.);
// 	tmp = pow( MAIN_RESP_TRANS-MAIN_RESP_FAC*T, (T-MAIN_RESP_TEMP_TRANS)/MAIN_RESP_TEMP_FAC )
// 			            * BetaN*Biomass*MAIN_RESP_CARB_PROPORT/upsilon;
// // 	return MyMax(1., pow( MAIN_RESP_TRANS-MAIN_RESP_FAC*T, (T-MAIN_RESP_TEMP_TRANS)/MAIN_RESP_TEMP_FAC )
// // 			            * BetaN*Biomass*MAIN_RESP_CARB_PROPORT/upsilon );
// 	
// // 	cout << tmp << endl;
// 	
// 	return tmp;
// }

double MaintRespFast(double Biomass, double BetaN, double upsilon, double T_fac)
{
	double tmp;

	tmp = T_fac * BetaN * Biomass * MAIN_RESP_CARB_PROPORT / upsilon;

// 	cout << "             " << tmp << endl;

	return tmp;
}

double MaintRespTmpFac(double T)
{
	T = MyMax(T, 20.);

	return pow(MAIN_RESP_TRANS - MAIN_RESP_FAC * T, (T - MAIN_RESP_TEMP_TRANS) / MAIN_RESP_TEMP_FAC);
}

//-----------------------------------------------------------------------------------------------

double WindAtHeight(double height, double ref_wind)
{
	double ws = 0.1;

	if (height > 1.3)
	{
// 		ws = ref_wind * log( (height-DISPLACEMENT_HEIGHT)/ROUGHNESS_LENGTH ) / 
// 					  log( (REF_HEIGHT_Z-DISPLACEMENT_HEIGHT)/ROUGHNESS_LENGTH );
		ws = ref_wind * log((height - DISPLACEMENT_HEIGHT) / ROUGHNESS_LENGTH) * WIND_AT_HEIGHT_HELPER;
	}

	return ws;
}

//-----------------------------------------------------------------------------------------------

