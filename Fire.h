#pragma once

using namespace ::std;

// This file contains routines to calculate fire intensity, patchiness,
// combustion and emissions.

double FireIntensity(double fuel, double fuel_moisture, double wind_speed);

double Patchiness(double intensity);

double ScorchHeight(double intensity);

double CombComplFine(double scorch);

double CombComplCoarse(double scorch);

double CombComplHeavy(double scorch);

double CombComplTopkill(double scorch, double height);

double CombComplTopkillHelper(double scorch);

double CombComplTopkillTree(double cc_topkill_helper, double height);

double N2Ograss(double bm);

double N2Oleaf(double bm);

double N2Ocoarse(double bm);

double N2Oheavy(double bm);

double N2Oshrub(double bm);

double CH4grass(double bm);

double CH4leaf(double bm);

double CH4coarse(double bm);

double CH4heavy(double bm);

double CH4shrub(double bm);

double LightFire(double dead_fuel, double live_fuel, double dead_fuel_moisture, double live_fuel_moisture, double wind_speed, int day);

void GetIgnitions(arryYear ignitions);

// -----------------------------------------------------------------------------------------------------------------
// Compute fire intensity for given biomass (fuel), moisture and wind speed

double FireIntensity(double fuel, double fuel_moisture, double wind_speed)
{
// 	return ( INT_FLUX*HEAT_YIELD_GRASS*fuel*( 1.+pow( wind_speed, INT_WIND_SPEED ) ) )/
// 		   ( INT_IGNITION_MOISTURE*fuel_moisture-INT_IGNITION_VEGETATION*( 1.- fuel_moisture ) );
// 	return ( INT_FLUX*HEAT_YIELD_GRASS*fuel*( wind_speed/(INT_WIND_SPEED+wind_speed) ) )/
// 		   ( INT_IGNITION_MOISTURE*fuel_moisture-INT_IGNITION_VEGETATION*( 1.- fuel_moisture ) );
	return (FIRE_H * 1000. * fuel * atan(wind_speed) * FIRE_C * 1000. * fuel / (1000. * fuel + FIRE_AW))
			/ (FIRE_QM * fuel_moisture + FIRE_QV * (1. - fuel_moisture));
}

double Patchiness(double intensity)
{
	return 0.989 - 0.595 * exp(-0.000837 * intensity);
}

double ScorchHeight(double intensity)
{
	return 20.1 * (1. - exp(-0.41 * intensity / 1000));
}

double CombComplFine(double scorch)
{
	return 0.884 / (1. + exp(0.01534 - 1.2833 * log(scorch)));
}

double CombComplCoarse(double scorch)
{
	return 0.884 / (1. + exp(2.3504 - 0.9838 * log(scorch))) - 0.00884;
}

double CombComplHeavy(double scorch)
{
	return 0.884 / (1. + exp(2.7184 - 1.0601 * log(scorch))) - 0.00884;
}

double CombComplTopkill(double scorch, double height)
{
	return 0.884 / (1. + exp(1.5669 - 0.7246 * log(scorch))) * (1.) / (1. + exp((height - 3.) / 1.));
}

double CombComplTopkillHelper(double scorch)
{
	return 0.884 / (1. + exp(1.5669 - 0.7246 * log(scorch)));
}

double CombComplTopkillTree(double cc_topkill_helper, double height)
{
	return cc_topkill_helper * 1. / (1. + exp((height - 3.)));
}

double N2Ograss(double bm)
{
	return bm * 0.44 * 0.0090 * 0.0076 * 1.57;
}
double N2Oleaf(double bm)
{
	return bm * 0.49 * 0.0127 * 0.0076 * 1.57;
}
double N2Ocoarse(double bm)
{
	return bm * 0.5 * 0.0074 * 0.0076 * 1.57;
}
double N2Oheavy(double bm)
{
	return bm * 0.5 * 0.0074 * 0.0076 * 1.57;
}
double N2Oshrub(double bm)
{
	return bm * 0.5 * 0.0093 * 0.0076 * 1.57;
}

double CH4grass(double bm)
{
	return bm * 0.44 * 0.0035 * 1.33;
}
double CH4leaf(double bm)
{
	return bm * 0.49 * 0.0035 * 1.33;
}
double CH4coarse(double bm)
{
	return bm * 0.5 * 0.0035 * 1.33;
}
double CH4heavy(double bm)
{
	return bm * 0.5 * 0.0035 * 1.33;
}
double CH4shrub(double bm)
{
	return bm * 0.5 * 0.0035 * 1.33;
}

// -----------------------------------------------------------------------------------------------------------------
// Function which "decides", if a fire can can break out. The function computes the amount of fuel,
// the resulting intensity and updates the biomasses, if the fire really spreads.

double LightFire(double dead_fuel, double live_fuel, double dead_fuel_moisture, double live_fuel_moisture, double wind_speed, int day)
{
	double ret = 0;

	double fuel = dead_fuel + live_fuel;
	double fuel_moisture = (live_fuel * live_fuel_moisture + dead_fuel * dead_fuel_moisture) / (live_fuel + dead_fuel);
	double fire_intensity = FireIntensity(fuel, fuel_moisture, wind_speed);

	if (MyRand() < IGNITION_PROB)
	{
		if (fire_intensity > IGNITION_MIN_INT)
		{
			ret = fire_intensity;
		}
	}

	return ret;
}

void GetIgnitions(arryYear ignitions, double param)
{
	for (int i = 0; i < 365; i++)
	{
		if (MyRand() < IGNITION_PAR_1 * param + IGNITION_PAR_2) ignitions[i] = 1.;
		else
			ignitions[i] = 0.;
	}

}

