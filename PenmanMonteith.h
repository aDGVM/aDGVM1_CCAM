#pragma once



// compute evapotranspiration according to the
// Penman-Monteith-equations.
// Input: P, T, Z, SGC, Tmax, Tmin, RHmax, RHmin

//-----------------------------------------------------------------------------

// psychromatic constant
double FOAgama(double P);

// slope of curve, expressing the relation between vapor saturation and temp.
double FOAs(double T);

// air density
double FOArho(double P, double T);

double FOAe0(double T);

double FOAeS(double Tmax, double Tmin);

double FOAVPD(double eA, double eS);

// Compute the evapotranspiration
double FOAPenMan(double Rn, double ga, double gn, double s, double cp, double gama, double rho, double VPD);

// Compute the reference evapotranspiration
double FOAPenManRef(double temp_actual, double temp_last, double delta, double radiation, double gamma, double wind, double e_mean, double e_actual);

//-----------------------------------------------------------------------------

double FOAgama(double P)
{
	return P * 0.000665;
}

double FOAs(double T)
{
	return 2504.0 * exp((17.27 * T) / (T + 237.2)) / pow(T + 237.3, 2);
}

double FOArho(double P, double T)
{
	return (P * 1e-3) / (SGC * (1.01 * (T + 273.3)));
}

double FOAe0(double T)
{
	return 0.6108 * exp((17.27 * T) / (T + 237.3));
}

double FOAeS(double Tmax, double Tmin)
{
	return (FOAe0(Tmax) + FOAe0(Tmin)) / 2.0;
}

double FOAVPD(double eA, double eS)
{
	return (eS - eA);
}

double FOAPenMan(double Rn, double ga, double gn, double s, double cp, double gama, double rho, double VPD)
{
// 	cout << setw(14) << Rn << setw(14) << ga << setw(14) << gn << setw(14) << s << setw(14) << cp << setw(14) << gama << setw(14) << rho << setw(14) <<  VPD << endl;

	if (gn <= 0.000000000001) return 0.0;
	else
		return (s * Rn + 86400. * rho * cp * VPD * ga) / (LAMBDA * (s + gama * (1.0 + ga / gn)));
}

double FOAPenManRef(double temp_actual, double temp_last, double delta, double radiation, double gamma, double wind, double e_mean, double e_actual)
{
	double temp_grad;

	temp_grad = 0.14 * (temp_actual - temp_last) * 0.15;

	return (0.408 * delta * (radiation - temp_grad) + gamma * (900. / (temp_actual + 273.)) * wind * (e_mean - e_actual)) / (delta + gamma * (1. + 0.34 * wind));
}

