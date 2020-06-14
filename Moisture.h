#pragma once

#include <vector>

using namespace std;

// This file contains functions to calculate rainfall sequneces
// and to update the soil bucket model.

// for each soil layer work out Gtheta; see Arora et al. 2003 for equations
// Output : GetST
void GetSoilTheta(double *diff_fc_wp, double *Theta, double *ThetaWP, double *GetST, const size_t soil_layers);

// Simulate rain sequendce for year
// Input  : days, rdo, pwet, ralpha, rbeta
// Output : RYear
void RainFallYear(int days, arry12 rdo, arry12 pwet, arry12 ralpha, arry12 rbeta, arryYear RYear);

// Read values of daily rainfall from file
// void RainFallYearFromFile( char *filename, arryYear RYear );

// Update water content in soil layers
// Input  : drain, theta, ThetaFC, depth
// Output : BIn
void BucketIn(double drain, double *Theta, double *ThetaFC, double *thickness, const size_t soil_layers);

// Input  : Et, Theta, ThetaWP, depth
// Output : BOut
void BucketOut(double Et, double *Theta, double *ThetaWP, double *thickness, const size_t soil_layers);

// ----------------------------------------------------------------------------------------------------------

void GetSoilTheta(double *diff_fc_wp, double *Theta, double *ThetaWP, std::vector<double> &GetST, const size_t soil_layers)
{
	for (size_t i = 0; i < soil_layers; i++)
	{
		double beta = MyMax(0., MyMin(1., (Theta[i] - ThetaWP[i]) * diff_fc_wp[i]));
		GetST[i] = 2. * beta - beta * beta;
	}
}

// This procedure generates a vector RYear, which contains the rain/day.
void RainFallYear(int days, arry12 rdo, arry12 pwet, arry12 ralpha, arry12 rbeta, arryYear RYear)
{
	/**
	 * This procedure is not used!
	 * Model data is read in directly from input data.
	 */
	for (int d = 1; d <= days; d++)
		if (MyRand() <= pwet[(int) floor(((double) d) / 30.42)]) MyRand();

	return;
}

void BucketIn(double drain, double *Theta, double *ThetaFC, double *thickness, const size_t soil_layers)
{
	double soilin = drain * 1e-3;			   	//rain input (m/day)

	if (soilin > 0)
	{
		for (size_t i = 0; i < soil_layers; i++)
		{
			const double SDi = (ThetaFC[i] - Theta[i]) * thickness[i]; 			//how much water can be put to soil layer (m)
			if (soilin > SDi)
			{
				Theta[i] = ThetaFC[i];
				soilin -= SDi;
				soilin = MyMax(0., soilin);
			}
			else
			{
				Theta[i] = Theta[i] + soilin / thickness[i];
				soilin = 0.;
				break;
			}
		}
	}
	return;
}

void BucketOut(double Et, double *Theta, double *ThetaWP, double *thickness, const size_t soil_layers)
{

	double soilout = Et * 1e-3; 							//convert mm/day into m/day

	if (soilout > 0)
	{
		for (size_t i = 0; i < soil_layers; i++)
		{
			const double SLi = (Theta[i] - ThetaWP[i]) * thickness[i]; 			//how much water can be lost from soil layer (m)
			if (soilout > SLi)
			{
				Theta[i] = ThetaWP[i];
				soilout -= SLi;
				soilout = MyMax(0., soilout);
			}
			else
			{
				Theta[i] = Theta[i] - soilout / thickness[i];
				soilout = 0.;
				break;
			}
		}
	}

	return;
}

