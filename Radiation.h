#pragma once

// ----------------------------------------------------------------------------

// Calculate inverse sun-earth distance for a given day in the year
double FOAdr(int DayofYear);

// Calculate declination angle of sun for a given day in year,
// the declination is the angle between the sun and the equator 
double FOAdelta(int DayofYear);

// Calculate angle at sunset as a function of latitude psi and
// declination delta
double FOAOmega(double psi, double delta);

// Calculate sun radiation as a function of solar constant gsc=0.082, the
// inverse sun-earth distance dr, the angle at sunset omega, the
// latitude psi and the declination delta
double FOARa(double dr, double omega, double psi, double delta);

// Calculate sun hours per day as a function of sun angle omega
double FOADaylighthours(double omega);

// Calculate solar radiation
double FOARs(double Ra, double n, double dlh);

// Calculate solar radiation at clear sky, no consideration of 
// sunshine duration per day
double FOARso(double Ra);


// Calculate net solar radiation given the radiation and albedo
double FOARns(double Rs);

// Calculate long wave radiation given radiation Rs, radiation at clear sky Rso,
// min/max temperature, vapor pressure and density of air
double FOARnl(double Rs, double Rso, double Tmax, double Tmin, double eA);

// Calculate net radiation, this value is used in the Model
void GetNetRadiation(double latitude, double psunshine, double Tmax, double Tmin, double eA, int day, double *Rn, double *SunHrs);

// Calculate photosynthetically active radiation (PAR)
double GetPARRadiation(double latitude, double psunshine, int day);

//-----------------------------------------------------------------------------

double FOAdr(int DayofYear)
{
	return 1.0 + 0.033 * cos(DayofYear * 2.0 * M_PI / 365.0);
}

double FOAdelta(int DayofYear)
{
	return 0.409 * sin(DayofYear * 2.0 * M_PI / 365.0 - 1.39);
}

double FOAOmega(double psi, double delta)
{
	return acos(-tan(psi) * tan(delta));
}

double FOARa(double dr, double omega, double psi, double delta)
{
	return (24.0 * 60.0 / M_PI) * GSC * dr * (omega * sin(psi) * sin(delta) + cos(psi) * cos(delta) * sin(omega));
}

double FOADaylighthours(double omega)
{
	return omega * 24.0 / M_PI;
}

double FOARs(double Ra, double n, double dlh)
{
	return (ANGSTRONG_A + ANGSTRONG_B * n / dlh) * Ra;
}

double FOARso(double Ra)
{
	return (ANGSTRONG_A + ANGSTRONG_B) * Ra;
}


double FOARns(double Rs)
{
	return Rs * (1.0 - ALBEDO);
}

double FOARnl(double Rs, double Rso, double Tmax, double Tmin, double eA)
{
	double TmaxK;
	double TminK;
	double RsRso;

	TmaxK = 273.16 + Tmax;
	TminK = 273.16 + Tmin;
	RsRso = MyMin(Rs / Rso, 1.0);

	return SBC * ((pow(TmaxK, 4.0) + pow(TminK, 4.0)) / 2.0) * (0.35 - 0.14 * sqrt(eA)) * (1.35 * RsRso - 0.35);
}

void GetNetRadiation(double latitude, double psunshine, double Tmax, double Tmin, double eA, int day, double *Rn, double *SunHrs)
{
	double psi;
	double dr;
	double delta;
	double omega;
	double Ra;
	double hrs;
	double Rs;
	double Rso;
	double Rnl;
	double Rns;

	psi = DegToRad(latitude);
	dr = FOAdr(day);
	delta = FOAdelta(day);
	omega = FOAOmega(psi, delta);
	Ra = FOARa(dr, omega, psi, delta);
	hrs = FOADaylighthours(omega);
	*SunHrs = psunshine * hrs;
	Rs = FOARs(Ra, *SunHrs, hrs);
	Rso = FOARso(Ra);
	Rnl = FOARnl(Rs, Rso, Tmax, Tmin, eA);
	Rns = FOARns(Rs);
	*Rn = (Rns - Rnl);                     // in MJ/m^2/day

	return;
}

double GetPARRadiation(double latitude, double psunshine, int day)
{
	double psi;
	double dr;
	double delta;
	double omega;
	double Ra;
	double hrs;
	double Rs;
	double sunhrs;

	psi = DegToRad(latitude);
	dr = FOAdr(day);
	delta = FOAdelta(day);
	omega = FOAOmega(psi, delta);
	Ra = FOARa(dr, omega, psi, delta);
	hrs = FOADaylighthours(omega);
	sunhrs = psunshine * hrs;
	Rs = FOARs(Ra, sunhrs, hrs);   // MJ/m^2/day

	// this converts from MJ/m^2/day to mmol/m^2/s
	return 5. * 1e6 / (sunhrs * 3600.) * Rs; //assumes 1 J/m2/s = 5 micro mol/m2/s photons

}

