#pragma once

#include "Radiation.h"
#include "PenmanMonteith.h"
#include "Leaf.h"
#include "Fxn.h"
#include "GlobalTypes.h"

using namespace std;

// This class stores all site specific input data

class clInData
{
public:
	clInData();
	~clInData()
	{
	}
	void calcAtmospheric(double latitude);
	void calcLeafPhoto(double ca_par_preassure, double soil_C);

	const size_t soil_layers_ = N_SOIL_LAYERS;

	double elev_;
	double atm_press_;
	double latitude_;
	double longitude_;

	double thickness_[N_SOIL_LAYERS];     // thickness of soil layers
	double depth_[N_SOIL_LAYERS];         // depth of soil layers

	double soil_N_[N_SOIL_LAYERS];
	double soil_C_[N_SOIL_LAYERS];
	double theta_[N_SOIL_LAYERS];
	double theta_wp_[N_SOIL_LAYERS];
	double theta_fc_[N_SOIL_LAYERS];
	double bulk_dens_[N_SOIL_LAYERS];

	arry12 tmp_;
	arry12 tmp_day_;
	arry12 tmp_min_;
	arry12 tmp_max_;
	arry12 ralpha_;
	arry12 rbeta_;
	arry12 reh_;
	arry12 pwet_;
	arry12 sun_;
	arry12 wnd_;
	arry12 rdo_;
	arry12 frost_;

	arry12 s12_;						// slope of vapour pressure curve kPa/C
	arry12 gama12_;						// psychrometric constant [kPa/C]
	arry12 eA12_;						// actual vapour pressure [kPa]
	arry12 eS12_;						// saturation vapour pressure [kPa]
	arry12 rho12_;						// density of air g/m3
	arry12 VPD12_;						// saturation vapour pressure deficit [kPa]
	arry12 Q012_;						// PAR
	arry12 A012C3_;
	arry12 RmL12C3_;
	arry12 A012C4_;
	arry12 RmL12C4_;
};

clInData::clInData()
{
	for (size_t i = 0; i < N_SOIL_LAYERS; i++)
	{
		bulk_dens_[i] = 0.;
		depth_[i] = 0.;
		soil_C_[i] = 0.;
		soil_N_[i] = 0.;
		theta_[i] = 0.;
		theta_fc_[i] = 0.;
		theta_wp_[i] = 0.;
		thickness_[i] = 0.;
	}

	elev_ = 0;
	atm_press_ = 0;
	latitude_ = 0;
	longitude_ = 0;

	for (size_t i = 0; i < 12; i++)
	{
		tmp_[i] = 0;
		tmp_day_[i] = 0;
		tmp_min_[i] = 0;
		tmp_max_[i] = 0;
		ralpha_[i] = 0;
		rbeta_[i] = 0;
		reh_[i] = 0;
		pwet_[i] = 0;
		sun_[i] = 0;
		wnd_[i] = 0;
		rdo_[i] = 0;
		frost_[i] = 0;
		s12_[i] = 0;
		gama12_[i] = 0;
		eA12_[i] = 0;
		eS12_[i] = 0;
		rho12_[i] = 0;
		VPD12_[i] = 0;
		Q012_[i] = 0;
		A012C3_[i] = 0;
		RmL12C3_[i] = 0;
		A012C4_[i] = 0;
		RmL12C4_[i] = 0;
	}
}

void clInData::calcAtmospheric(double latitude)
{
	int midmonthday;

	// compute some atmospheric values of the current coordinates
	for (int i = 0; i < 12; i++)
	{
		midmonthday = (int) floor((double) (i + 1.) * 30.42 - 15.21);
		// average saturation vapor pressure, KPa
		eS12_[i] = FOAeS(tmp_max_[i], tmp_min_[i]);
		// FOAe0 = saturation vapor pressure, reh = rel. hum in % => reh/100 in KPa
		// actual actual vapor pressure
		eA12_[i] = reh_[i] / 100. * eS12_[i];  // ((FOAe0(tmax[i])+FOAe0(tmin[i]))/2.);
		// Radiation in MJ/m^2/day
		s12_[i] = FOAs(tmp_[i]); 							// KPa/degC
		// psychochromatic constant 
		gama12_[i] = FOAgama(atm_press_ / 1000.);  				// KPa/degC
		// vapor pressure difference
		VPD12_[i] = FOAVPD(eA12_[i], eS12_[i]); 				// KPa
		// density of air
		rho12_[i] = FOArho(atm_press_, tmp_[i]);				// g/m3 //takes atm_press in Pascals!
		// PAR
		Q012_[i] = GetPARRadiation(latitude, sun_[i], midmonthday);
// 		cout << VPD12_[i] << endl;
// 		cout << "------- " << Q012_[i] << endl;
	}

}

void clInData::calcLeafPhoto(double ca_par_preassure, double soil_C)
{
	double A0;
	double RmL;

	// get the photosynthetic rate A0 and the leaf maintainance respiration
	// for the given atmospheric values and for C3/C4-plants.
	for (int i = 0; i < 12; i++)
	{
		GetC3A(tmp_day_[i], soil_N_[0], soil_C, ca_par_preassure, OI_PAR_PREASSURE, Q012_[i], R_MAINT_RESP_TR[0], ABS_PHOTONS_C3, ALPHA_C3, reh_[i], atm_press_,
				WindAtHeight(VEGETATION_HEIGHT, wnd_[i]), CLD_TREE, M_TREE[0], B_TREE[0], TAU_K25_C3, TAU_Q10_C3, KO_K25_C3, KO_Q10_C3, KC_K25_C3, KC_Q10_C3,
				&A0, &RmL);
		A012C3_[i] = A0;
		RmL12C3_[i] = RmL;

		GetC4A(tmp_day_[i], soil_N_[0], soil_C, ca_par_preassure, KAPPA_C4, Q012_[i], R_MAINT_RESP_GR[0], ABS_PHOTONS_C4, ALPHAR_F_C4, reh_[i], atm_press_,
				WindAtHeight(VEGETATION_HEIGHT, wnd_[i]), CLD_GRASS, M_GRASS[0], B_GRASS[0], &A0, &RmL);

		A012C4_[i] = A0;
		RmL12C4_[i] = RmL;
	}

	return;
}
