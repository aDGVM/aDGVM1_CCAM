#pragma once

#include "Globals.h"
#include "Fxn.h"

using namespace std;

class clGrass
{
private:
	int Dpos_; 								// dormancy counter (days where nCGT>0)
	int Dneg_;								// standby counter (days where nCGT=0)
	int Dormant_;							// dormant yes=1, no=0
	int grass_type_;						// C4: open grass (1), sav canopy grass (0) or for. can. grass (2)
											// C3: open grass (4), sav canopy grass (3) or for. can. grass (5)
	int active_days_;						// counts active days per year
	double Bl_;    							// leaf biomass
	double Br_;    							// root biomass
	double Bs_;    							// stem biomass

	double Bld_;								// dead leaf biomass, kg per m^2
	double Bsd_;								// dead stem biomass, kg per m^2
	double Brd_;								// dead root biomass, kg per m^2

	double leaf_combustion_;
	double stem_combustion_;

	double gpp_;
	double Rma_;
	double Rgr_;

	double Qsum_;
	double Qi_;
	double Ci_;
	double alloc_denom_;
	double gc_;    								// canopy stomatal conductance (umol/m2/s)
	double gb_;									// canopy boundary layer conductanced (m/s)    non-molar
	double fuel_moisture_;

	double Droot_;   								// rooting depth m
	double nCGT_;
	double Gw_;
	double Aindex_;
	double height_;								// height m
	double canopy_radius_;
	double canopy_area_;							// canopy area m^2
	double LAI_;   	 							// leaf area index of grass

	double A_factor_;		// factor to scale A0 if the system is grazed

	double Bl_year_[365];
	double Br_year_[365];

private:
	void calDroot();
	void calPlantHeight();
	void calCanopyRadius();
	void calCanopyArea();
	void calPlantLAI();
	void calSoilMoisture(double *G_theta, double *thickness);
	void calGc(double A, double resp, double Ca, double P, double hs, double T);
	void calGb(double wind);
	void AllocateCorbon();
	void UpdateAllometry();
	void MassTurnover();
	void Respirate(double rRoot, double rStem, double rLeaf);
	double AllocLeaf();
	double AllocRoot();
	double calQSum();

public:
	clGrass(double init_mass, int grass_type);
	~clGrass();

	void RunPhysiology(double A0_C4, double A0_C3, double RmL_C4, double RmL_C3, double mmsTOkgd, double T, double wind, double *G_theta, double Ca, double P,
			double rh, double mean_sav_height, double mean_for_height, double dead_bm, double T_fac, int frost, double *thickness, int day);
	double getEt(double T, double P, double radnet_day, double s12, double SP_HEAT, double gama12, double rho12, double VPD12);
	void setStateAfterFire(double patchiness, double cc_fine);

	int getDormant()
	{
		return Dormant_;
	}
	int getActiveDays()
	{
		return active_days_;
	}
	double getBl()
	{
		return Bl_;
	}
	double getBr()
	{
		return Br_;
	}
	double getBs()
	{
		return Bs_;
	}
	double getBld()
	{
		double tmp = Bld_;
		Bld_ = 0;
		return tmp;
	} // get dead biomass and set to pool to zero
	double getBsd()
	{
		double tmp = Bsd_;
		Bsd_ = 0;
		return tmp;
	} // get dead biomass and set to pool to zero
	double getBrd()
	{
		double tmp = Brd_;
		Brd_ = 0;
		return tmp;
	} // get dead biomass and set to pool to zero
	double getLeafCombustion()
	{
		return leaf_combustion_;
	}
	double getStemCombustion()
	{
		return stem_combustion_;
	}
	double getGPP()
	{
		return gpp_;
	}
	double getRma()
	{
		return Rma_;
	}
	double getRgr()
	{
		return Rgr_;
	}
	double getQsum()
	{
		return Qsum_;
	}
	double getQi()
	{
		return Qi_;
	}
	double getGc()
	{
		return gc_;
	}
	double getGb()
	{
		return gb_;
	}
	double getFuelMoisture()
	{
		return fuel_moisture_;
	}
	double getDroot()
	{
		return Droot_;
	}
	double getnCGT()
	{
		return nCGT_;
	}
	double getGw()
	{
		return Gw_;
	}
	double getAindex()
	{
		return Aindex_;
	}
	double getHeight()
	{
		return height_;
	}
	double getCanopyRadius()
	{
		return canopy_radius_;
	}
	double getCanopyArea()
	{
		return canopy_area_;
	}
	double getLAI()
	{
		return LAI_;
	}
	double getBlYearMax();
	double getBrYearMax();

};

// -------------------------------------------------------------------------------------------------------------------
clGrass::clGrass(double init_mass, int grass_type)
{

	Dpos_ = 0;
	Dneg_ = 0;
	Dormant_ = 1;
	grass_type_ = grass_type;
	active_days_ = 0;

	Bl_ = init_mass * A0_LEAF_GRASS[grass_type_];
	Br_ = init_mass * A0_ROOT_GRASS[grass_type_];
	Bs_ = 0.001;
	Bld_ = 0.;
	Brd_ = 0.;
	Bsd_ = 0.;
	leaf_combustion_ = 0.;
	stem_combustion_ = 0.;
	gpp_ = 0.;
	Rma_ = 0.;
	Rgr_ = 0.;
	Qsum_ = 1.;
	Qi_ = 1.;
	gc_ = 0.;
	gb_ = 0.;
	fuel_moisture_ = 1.;
	nCGT_ = 0.;
	Gw_ = 1.;
	Aindex_ = 1.;

	A_factor_ = 1.;

	for (int i = 0; i < 365; i++)
	{
		Bl_year_[i] = 0;
		Br_year_[i] = 0;
	}

	UpdateAllometry();
}

// -------------------------------------------------------------------------------------------------------------------
clGrass::~clGrass()
{
}

// -------------------------------------------------------------------------------------------------------------------
void clGrass::calDroot()
{
	Droot_ = MyMin(MAX_ROOT_DEP_GRASS, CylinderDepth(Br_, MIN_ROOT_RAD_GRASS, ROOT_DENSITY_GRASS, MAX_ROOT_DEP_GRASS));
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// Standard power function, based on fits to Niklas&Enquist data, estimates the plant height from the stem biomass. Arora
void clGrass::calPlantHeight()
{
	height_ = MIN_HEIGHT + HEIGHT_C1_GRASS[grass_type_] * pow(Bl_, HEIGHT_C2_GRASS[grass_type_]);
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// Compute radius of the canopy
void clGrass::calCanopyRadius()
{
	canopy_radius_ = height_ * GAMMA_CANOPY_GRASS[grass_type_];
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// Compute the canopy area of a Grass,
// assumes canopy radius is constant proportion of height
void clGrass::calCanopyArea()
{
	canopy_area_ = M_PI * canopy_radius_ * canopy_radius_;
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// compute leaf area index
void clGrass::calPlantLAI()
{
	LAI_ = Bl_ * SLA_GRASS[grass_type_] + MIN_LAI;

	return;
}

// -------------------------------------------------------------------------------------------------------------------
// compute the allometric values of the Grass
void clGrass::UpdateAllometry()
{
	calDroot();
	calPlantHeight();
	calCanopyRadius();
	calCanopyArea();
	calPlantLAI();

	return;
}

// -------------------------------------------------------------------------------------------------------------------
// allocation function of net carbon gain to the roots
double clGrass::AllocRoot()
{
	return (1. + A0_ROOT_GRASS[grass_type_] - Gw_) * alloc_denom_;
}

// -------------------------------------------------------------------------------------------------------------------
// allocation function of net carbon gain to the leaf
double clGrass::AllocLeaf()
{
	return (1. - Ci_) * alloc_denom_;
}

// -------------------------------------------------------------------------------------------------------------------
// alloaction functions for leaf, stem and root, and add it to the Grass biomass
void clGrass::AllocateCorbon()
{
	if (Dormant_ == 0)
	{
		alloc_denom_ = 1. / (2. + A0_ROOT_GRASS[grass_type_] - Gw_ - Ci_);
		Bl_ += (gpp_ - Rgr_) * AllocLeaf();
		Br_ += (gpp_ - Rgr_) * AllocRoot();
	}

	return;
}

// ----------------------------------------------------------------------------------------------------------------------
// run the grass physiology
void clGrass::RunPhysiology(double A0_C4, double A0_C3, double RmL_C4, double RmL_C3, double mmsTOkgd, double T, double wind, double *G_theta, double Ca,
		double P, double rh, double mean_sav_height, double mean_for_height, double dead_bm, double T_fac, int frost, double *thickness, int day)
{
	double A0;
	double Acs;		// canopy and canopy-stressed photosythesis
	double RmS;
	double RmR;		// maint respiration total, stem, root
	double RmLc;		// canopy leaf maint resp
	double RmLg;		// leaf maint resp in kg/day/plant
	double RmL;

	if (day == 0)
	{
		// 		if (number_==0 )cout << setw(14) << active_days_ << endl;
		active_days_ = 0;
	}

	if (grass_type_ <= 2)
	{
		A0 = A0_C4;
		RmL = RmL_C4;
	}
	else
	{
		A0 = A0_C3;
		RmL = RmL_C3;
	}

	gpp_ = 0.;
	Rma_ = 0.;
	Rgr_ = 0.;

	// water availability
	calSoilMoisture(G_theta, thickness);

	// light availability
	Qi_ = 1.;
	if (dead_bm > Bl_) Qi_ *= LICMP_1[grass_type_ + 2][grass_type_ + 2] / dead_bm * Bl_ + LICMP_2[grass_type_ + 2][grass_type_ + 2]; // dead grass shading effects

	if ((grass_type_ == GR_C4_SAV || grass_type_ == GR_C3_SAV) && mean_sav_height > height_) // savanna tree sub-canopy grass
	Qi_ *= MyMax(0.01, LICMP_1[TR_SAV][grass_type_ + 2] / mean_sav_height * height_ + LICMP_2[TR_SAV][grass_type_ + 2]);

	if ((grass_type_ == GR_C4_FOR || grass_type_ == GR_C3_FOR) && mean_for_height > height_) // forest tree sub-canopy grass
	Qi_ *= MyMax(0.01, LICMP_1[TR_FOR][grass_type_ + 2] / mean_for_height * height_ + LICMP_2[TR_FOR][grass_type_ + 2]);

	Qsum_ = Qi_ * calQSum();

	Ci_ = MyMin(1., Bl_ / (Bl_ + Br_) / A0_LEAF_GRASS[grass_type_]);

	// conopy photosynthesis
	A0 *= A_factor_;   // A_factor only for grazing
	Acs = A0 * Qsum_ * Gw_;			// water stressed canopy photosynthesis

	// carbon gain
	gpp_ = Acs * mmsTOkgd * (double) (1 - Dormant_);   //water and light stressed carbon gain, kg/day/m^2, gpp=0 if plant is dormant

	// components of maintenance respiration

	RmLc = R_MAINT_RESP_GR[grass_type_] * Acs; //maint resp leaf is now function of water stress and canopy photosynthesis
	RmLg = RmLc * mmsTOkgd;       //leaf maint resp in kg/day/m^2
	RmS = 0;
	RmR = RmLg + MaintRespFast(BETA_ROOT_GRASS * Br_, BETA_N, UPSILON_ROOT_GRASS[grass_type_], T_fac);

	Rgr_ = gpp_ * SIGMA_GROW_RESP_GRASS[grass_type_];   //growth resp, kg/day/m^2
	nCGT_ = gpp_ - Rgr_ - RmR - RmS - RmLg;

	//canopy conductance
	calGb(wind);  // m/s
	calGc(Acs, RmLc, Ca, P, rh, T);  // (mmol/m^2/s)

	double Ti;
	if (tmp_min_month > TI_CONST) Ti = 0.;
	else
		Ti = 2.0 * (-1. + tmp_min_month / TI_CONST);

	Aindex_ = A0 * ((WATER_IND_GRASS[grass_type_] * G_theta[3] + WATER_IND_GRASS_1[grass_type_] * Gw_) + Ti) - RmL;

	// if no carbon gain increment days without carbon gain
	if (Dormant_ == 0 && Aindex_ <= STRESS_INDEX_GRASS[grass_type_]) Dneg_++;
	// if carbon gain then reset days without carbon gain counter
	if (Dormant_ == 0 && Aindex_ > STRESS_INDEX_GRASS[grass_type_]) Dneg_ = 0;
	// if dormant and carbon gain then increment dormancy break counter
	if (Dormant_ == 1 && Aindex_ > STRESS_INDEX_GRASS[grass_type_]) Dpos_++;
	// if dormant and carbon gain = 0 then reset dormancy break counter
	if (Dormant_ == 1 && Aindex_ <= STRESS_INDEX_GRASS[grass_type_]) Dpos_ = 0;
	// days with frost increase number of neg. days
	if (frost == 1)
	{
		Dneg_++;
		Dpos_ = 0;
	}
	// if in growth stage then dpos (dormancy break counter) remains zero.
	if (Dormant_ == 0) Dpos_ = 0;
	if (Dormant_ == 1) Dneg_ = 0;

	//if days of negative nCGT > threshold then assign to dormant state
	if (Dneg_ >= D_NEG_GRASS[grass_type_])
	{
		fuel_moisture_ = 1.;
		Dormant_ = 1;
		Dneg_ = 0;
		Bld_ += (1. - REM_BM_GRASS) * Bl_;
		Bl_ = REM_BM_GRASS * Bl_;
	}

	//if days of positive nCGT > threshold
	if (Dpos_ >= D_POS_GRASS[grass_type_])
	{
		Dormant_ = 0; //assign to growth state
		Dpos_ = 0;
		A_factor_ = 1.;
	}

	active_days_ += (1 - Dormant_);
	// 	if ( number_==0 ) cout << setw(14) << active_days_ << setw(14) << Dneg_ << setw(14) << Dpos_ << setw(14) << Dormant_ << endl;

	// maintainance respiration
	Respirate(RmR, RmS, RmLg);

	// biomass turnover
	MassTurnover();

	// allocate carbon gain
	AllocateCorbon();

	// update allometry
	UpdateAllometry();

#	ifdef SITE_GRAZING
	if ( GLOB_YEAR <= 187 )
	{
		Bl_ *= SITE_GRAZING_RATE_1/100.;
		A_factor_ *= SITE_GRAZING_RATE_1/100.*0.01+0.99;
	}
	else
	{
		Bl_ *= SITE_GRAZING_RATE_2/100.;
		A_factor_ *= SITE_GRAZING_RATE_2/100.*0.01+0.99;
	}
#	endif

	for (int i = 0; i < 364; i++)
	{
		Bl_year_[i] = Bl_year_[i + 1];
		Br_year_[i] = Br_year_[i + 1];

		Bl_year_[364] = Bl_;
		Br_year_[364] = Br_;
	}

	fuel_moisture_ *= DESIC_COEFF[grass_type_];
// 	cout << fuel_moisture_ << endl;

	return;
}

// -------------------------------------------------------------------------------------------------------------------
// compute the water aviability of a Grass. The value is between 0 and 1 and is influenced by the
// soil moisture and the rooting depth
void clGrass::calSoilMoisture(double *G_theta, double *thickness)
{
	int layers = 1;    // counts the number of layers, in which the plant has roots
	double depth = thickness[0];
	Gw_ = G_theta[0] * thickness[0];

	while (Droot_ > depth)
	{
		Gw_ += (G_theta[layers] * thickness[layers]);
		depth += thickness[layers];
		layers++;
	}

	Gw_ /= depth;

	return;
}

// --------------------------------------------------------------------------------------------------------------
// Turnover from live leaf and root biomass to dead standing leaf biomass
// and dead root bimoass.
void clGrass::MassTurnover()
{
	if (Dormant_ == 0)
	{
		Bld_ += HELPER_BLS_GRASS[grass_type_] * Bl_;
		Bl_ *= HELPER_BL__GRASS[grass_type_];
	}

// 	Bsd_ += HELPER_BSD_GRASS*Bs_;
// 	Bs_  *= HELPER_BS__GRASS;

	Brd_ += HELPER_BRD_GRASS * Br_;
	Br_ *= HELPER_BR__GRASS;

	return;
}

// --------------------------------------------------------------------------------------------------------------

void clGrass::Respirate(double rRoot, double rStem, double rLeaf)
{
	double tmp_bm = Br_ + Bl_;

	Br_ = MyMax(0.0001, Br_ - rRoot);
//     Bs_  = 0.001;
	Bl_ = MyMax(0.0001, Bl_ - rLeaf);

	Rma_ = tmp_bm - (Br_ + Bl_);

	return;
}

// calculate gc_ in molar units
void clGrass::calGc(double A, double resp, double Ca, double P, double hs, double T)
{
	double cs;
	double gb = gMolar(T, gb_, P);

	if (gb > 0)
	{
		cs = Ca - CS_FACTOR * A * P / gb; 							//this line follows Arora
		gc_ = MyMax(0., M_GRASS[grass_type_] * (A - resp) * hs / 100. * P / cs + B_GRASS[grass_type_]);	//ie for negative An we produce zero gs
	}
	else
		gc_ = 0.;

// 	cout << grass_type_ << " " << hs << " " << gc_ << endl;
	return;
}

void clGrass::calGb(double wind)
{
	gb_ = GetgbCANOPY(wind, height_);

	return;
}

void clGrass::setStateAfterFire(double patchiness, double cc_fine)
{
	leaf_combustion_ = Bl_ * patchiness * cc_fine;
	Bl_ -= leaf_combustion_;

	Dneg_ = 0;
	Dpos_ = 0;
	UpdateAllometry();
	Dormant_ = 1;

	return;
}

double clGrass::getEt(double T, double P, double radnet_day, double s12, double SP_HEAT, double gama12, double rho12, double VPD12)
{
	if (Dormant_ == 1) return 0;

	double gc_nmolar = gNonMolar(T, gc_, P);  //gn

	return FOAPenMan(radnet_day, gb_, gc_nmolar, s12, SP_HEAT, gama12, rho12, VPD12);  //mm/day
}

double clGrass::getBlYearMax()
{
	double max = 0;
	for (int i = 0; i < 365; i++)
	{
		if (Bl_year_[i] > max) max = Bl_year_[i];
	}
	return max;
}

double clGrass::getBrYearMax()
{
	double max = 0;
	for (int i = 0; i < 365; i++)
	{
		if (Br_year_[i] > max) max = Br_year_[i];
	}
	return max;
}

double clGrass::calQSum()
{
	// NOTE: inserted 0.4 for callibration (Qsum=1 for LAI=1)
	return (1. - exp(-K_CAN_EXT_GRASS[grass_type_] * (LAI_ + 0.4))) * K_CAN_EXT_GRASS_INVERSE[grass_type_];
}

