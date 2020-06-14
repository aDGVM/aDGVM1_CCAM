#pragma once
#include <random>

#define MODMULT(a, b, c, m, s) q = s/a; s = b*(s-a*q)-c*q; if (s < 0) s += m; 

std::linear_congruential_engine<uint64_t, 25214903917, 11, 281474976710656> gen;

inline double MyRand()
{
	return std::generate_canonical<double, 48>(gen);
}
// Minimum of two values
double MyMin(double x, double y);

// Maximum of two values
double MyMax(double x, double y);

// Transform degree to rad
double DegToRad(double lat);



// -----------------------------------------------------------------------------

double MyMin(double x, double y)
{
	return x < y ? x : y;
}

double MyMax(double x, double y)
{
	return x > y ? x : y;
}

// -----------------------------------------------------------------------------

double DegToRad(double lat)
{
	return lat * M_PI / 180.0;
}

// -----------------------------------------------------------------------------
