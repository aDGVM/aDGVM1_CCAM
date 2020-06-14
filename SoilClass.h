#pragma once

// This class implements the Yasso soil model
// following Liski et al 2005 Ecological Modelling

#include "SoilClassGlobals.h"

using namespace std;

class clSoil
{
private:
	double x_fwl_;      // fine woody litter
	double x_cwl_;      // coarse woody litter
	double x_ext_;      // extractives
	double x_cel_;      // celluloses
	double x_lig_;      // lignin-like compounds
	double x_hum1_;     // humus
	double x_hum2_;     // more recalcitrant humus

	double r_ext_;      // CO2 release from extractives
	double r_cel_;      // CO2 release from cellulose
	double r_lig_;      // CO2 release from lignin-like
	double r_hum1_;     // CO2 release from humus1
	double r_hum2_;     // CO2 release from humus2

	double u_nwl_;      // non-woody litter input
	double u_fwl_;      // fine woodly litter input
	double u_cwl_;      // coarse woody litter input

	double sens_hum_;   // tmp and drought sensitivity of k for humus
	double sens_oth_;   // tmp and drought sensitivity of k and a for other compartments

public:
	clSoil();
	clSoil(double soil_carbon);
	~clSoil()
	{
	}
	;

	double UpdateCarbonPools(double u_nwl, double u_fwl, double u_cwl, double T, double D);
	double GetCarbonInput();
	double GetCarbonStored();
	double GetCarbonRelease();
	double GetCarbonInHumus();
	void PrintCarbonPools();

	double GetXfwl()
	{
		return x_fwl_;
	}
	double GetXcwl()
	{
		return x_cwl_;
	}
	double GetXext()
	{
		return x_ext_;
	}
	double GetXcel()
	{
		return x_cel_;
	}
	double GetXlig()
	{
		return x_lig_;
	}
	double GetXhu1()
	{
		return x_hum1_;
	}
	double GetXhu2()
	{
		return x_hum2_;
	}
	double GetUnwl()
	{
		return u_nwl_;
	}
	double GetUfwl()
	{
		return u_fwl_;
	}
	double GetUcwl()
	{
		return u_cwl_;
	}
	double GetRext()
	{
		return r_ext_;
	}
	double GetRcel()
	{
		return r_cel_;
	}
	double GetRlig()
	{
		return r_lig_;
	}
	double GetRhu1()
	{
		return r_hum1_;
	}
	double GetRhu2()
	{
		return r_hum2_;
	}

private:
	void getSensHum(double T, double D);
	void getSensOth(double T, double D);
};

clSoil::clSoil()
{
	sens_hum_ = 0.;
	sens_oth_ = 0.;
	x_fwl_ = 0.;    // fine woody litter
	x_cwl_ = 0.;    // coarse woody litter
	x_ext_ = 0.;    // extractives
	x_cel_ = 0.;    // celluloses
	x_lig_ = 0.;    // lignin-like compounds
	x_hum1_ = 1.;    // humus
	x_hum2_ = 1.;    // more recalcitrant humus

	r_ext_ = 0.;     // CO2 release from extractives
	r_cel_ = 0.;     // CO2 release from cellulose
	r_lig_ = 0.;     // CO2 release from lignin-like
	r_hum1_ = 0.;     // CO2 release from humus1
	r_hum2_ = 0.;     // CO2 release from humus2

	u_nwl_ = 0.;     // non-woody litter input
	u_fwl_ = 0.;     // fine woodly litter input
	u_cwl_ = 0.;     // coarse woody litter input
}

clSoil::clSoil(double soil_carbon)
{
	sens_hum_ = 0.;
	sens_oth_ = 0.;
	x_fwl_ = soil_carbon * 0.01152994;    // fine woody litter
	x_cwl_ = soil_carbon * 0.06949435;    // coarse woody litter
	x_ext_ = soil_carbon * 0.00851789;    // extractives
	x_cel_ = soil_carbon * 0.06624854;    // celluloses
	x_lig_ = soil_carbon * 0.06360159;    // lignin-like compounds
	x_hum1_ = soil_carbon * 0.2377162;     // humus
	x_hum2_ = soil_carbon * 0.5428915;     // more recalcitrant humus

	r_ext_ = 0.;     // CO2 release from extractives
	r_cel_ = 0.;     // CO2 release from cellulose
	r_lig_ = 0.;     // CO2 release from lignin-like
	r_hum1_ = 0.;     // CO2 release from humus1
	r_hum2_ = 0.;     // CO2 release from humus2

	u_nwl_ = 0.;     // non-woody litter input
	u_fwl_ = 0.;     // fine woodly litter input
	u_cwl_ = 0.;     // coarse woody litter input
}

// ---------------------------------------------------------------
// update carbon pools (run yasso)
// function takes input values in biomass (kg/m^2) and
// transforms it to C (kg/m^2)
double clSoil::UpdateCarbonPools(double u_nwl, double u_fwl, double u_cwl, double T, double D)
{
	getSensHum(T, D);
	getSensOth(T, D);

	u_nwl_ = 0.44 * u_nwl;     // non-woody litter input
	u_fwl_ = 0.44 * u_fwl;     // fine woodly litter input
	u_cwl_ = 0.44 * u_cwl;     // coarse woody litter input

	r_ext_ = (1. - YASSO_P_EXT) * YASSO_K_EXT * x_ext_;
	r_cel_ = (1. - YASSO_P_CEL) * YASSO_K_CEL * x_cel_;
	r_lig_ = (1. - YASSO_P_LIG) * YASSO_K_LIG * x_lig_;
	r_hum1_ = (1. - YASSO_P_HUM1) * YASSO_K_HUM1 * x_hum1_;
	r_hum2_ = YASSO_K_HUM2 * x_hum2_;

	x_fwl_ += (u_fwl_ - sens_oth_ * YASSO_A_FWL * x_fwl_);

	x_cwl_ += (u_cwl_ - sens_oth_ * YASSO_A_CWL * x_cwl_);

	x_ext_ += (u_nwl_ * YASSO_C_NWL_EXT + YASSO_C_FWL_EXT * sens_oth_ * YASSO_A_FWL * x_fwl_ + YASSO_C_CWL_EXT * sens_oth_ * YASSO_A_CWL * x_cwl_
			- sens_oth_ * YASSO_K_EXT * x_ext_);

	x_cel_ += (u_nwl_ * YASSO_C_NWL_CEL + YASSO_C_FWL_CEL * sens_oth_ * YASSO_A_FWL * x_fwl_ + YASSO_C_CWL_CEL * sens_oth_ * YASSO_A_CWL * x_cwl_
			- sens_oth_ * YASSO_K_CEL * x_cel_);

	x_lig_ += (u_nwl_ * YASSO_C_NWL_LIG + YASSO_C_FWL_LIG * sens_oth_ * YASSO_A_FWL * x_fwl_ + YASSO_C_CWL_LIG * sens_oth_ * YASSO_A_CWL * x_cwl_
			- sens_oth_ * YASSO_K_LIG * x_lig_ + YASSO_P_EXT * sens_oth_ * YASSO_K_EXT * x_ext_ + YASSO_P_CEL * sens_oth_ * YASSO_K_CEL * x_cel_);

	x_hum1_ += (YASSO_P_LIG * sens_oth_ * YASSO_K_LIG * x_lig_ - sens_hum_ * YASSO_K_HUM1 * x_hum1_);

	x_hum2_ += (YASSO_P_HUM1 * sens_hum_ * YASSO_K_HUM1 * x_hum1_ - sens_hum_ * YASSO_K_HUM2 * x_hum2_);

	return GetCarbonRelease();
}

double clSoil::GetCarbonInput()
{
	return u_nwl_ + u_fwl_ + u_cwl_;
}

double clSoil::GetCarbonStored()
{
	return x_fwl_ + x_cwl_ + x_ext_ + x_cel_ + x_lig_ + x_hum1_ + x_hum2_;
}

double clSoil::GetCarbonRelease()
{
	return r_ext_ + r_cel_ + r_lig_ + r_hum1_ + r_hum2_;
}

// ------------------------------------------------------------------------
// returns soil carbon stroed in humus, in kgC/m^2
double clSoil::GetCarbonInHumus()
{
	return r_hum1_ + r_hum2_;
}

void clSoil::PrintCarbonPools()
{
	cout << setw(14) << x_fwl_ << setw(14) << x_cwl_ << setw(14) << x_ext_ << setw(14) << x_cel_ << setw(14) << x_lig_ << setw(14) << x_hum1_ << setw(14)
			<< x_hum2_ << setw(14) << u_nwl_ << setw(14) << u_fwl_ << setw(14) << u_cwl_ << setw(14) << r_ext_ << setw(14) << r_cel_ << setw(14) << r_lig_
			<< setw(14) << r_hum1_ << setw(14) << r_hum2_ << endl;

	return;
}

void clSoil::getSensHum(double T, double D)
{
	sens_hum_ = 1. + YASSO_S_HUM * YASSO_BETA * (T - YASSO_T0) + YASSO_GAMMA * D;
}

void clSoil::getSensOth(double T, double D)
{
	sens_oth_ = 1. + YASSO_BETA * (T - YASSO_T0) + YASSO_GAMMA * D;
}

