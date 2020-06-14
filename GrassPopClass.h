#pragma once

#include "GrassClass.h"
#include <vector>

class clGrassPop
{
private:
	int pop_size_;					// counter
	std::vector<clGrass> Grasses;					// array of plants
	double root_bm_live_;				// kg per m^2
	double stem_bm_live_;				// kg per m^2
	double leaf_bm_live_;				// kg per m^2
	double root_bm_dead_;				// kg per m^2
	double stem_bm_dead_st_;			// kg per m^2
	double stem_bm_dead_ly_;			// kg per m^2
	double leaf_bm_dead_st_;			// kg per m^2
	double leaf_bm_dead_ly_;			// kg per m^2
	double leaf_live_comb_;
	double leaf_dead_st_comb_;
	double leaf_dead_ly_comb_;
	double soil_nwl_;					// soil non-woody litter (foilage, fine roots)
	double gpp_;
	double Rma_;
	double Rgr_;
	double C34_ratio_;                 // ratio between C3 and C4 grasses
	double cover_[6];                  // cover of different grass types

	double C34_ratio_year_[3650];
	double active_days_;

public:
	clGrassPop();
	~clGrassPop();
	void addGrass(double init_mass, int grass_type);
	void setStateAfterFire(double patchiness, double cc_fine);
	void emptyCombustionPools();
	void RunPhysiology(double A0_C4, double A0_C3, double RmL_C4, double RmL_C3, double mmsTOkgd, double T, double wind, double *G_theta, double Ca, double P,
			double rh, double mean_sav_height, double mean_for_height, double T_fac, int frost, double p_can_sav, double p_can_for, double *thickness, int day);
	double getEt(double T, double P, double radnet_day, double s12, double SP_HEAT, double gama12, double rho12, double VPD12, double p_canopy);
	double getWetBiomassForFire();
	double getDryBiomassForFire();

	double getLeafBmLive()
	{
		return leaf_bm_live_;
	}			// kg/m^2
	double getStemBmLive()
	{
		return stem_bm_live_;
	}			// kg/m^2
	double getRootBmLive()
	{
		return root_bm_live_;
	}			// kg/m^2

	double getLeafBmDeadSt()
	{
		return leaf_bm_dead_st_;
	}		// kg/m^2
	double getLeafBmDeadLy()
	{
		return leaf_bm_dead_ly_;
	}		// kg/m^2
	double getStemBmDeadSt()
	{
		return stem_bm_dead_st_;
	}		// kg/m^2
	double getStemBmDeadLy()
	{
		return stem_bm_dead_ly_;
	}		// kg/m^2
	double getRootBmDead()
	{
		return root_bm_dead_;
	}			// kg/m^2
	double getLeafLiveCombustion()
	{
		return leaf_live_comb_;
	}			// kg/m^2
	double getLeafDeadStCombustion()
	{
		return leaf_dead_st_comb_;
	}		// kg/m^2
	double getLeafDeadLyCombustion()
	{
		return leaf_dead_ly_comb_;
	}		// kg/m^2
	double getGPP()
	{
		return gpp_;
	}					// kg/m^2
	double getRma()
	{
		return Rma_;
	}					// kg/m^2
	double getRgr()
	{
		return Rgr_;
	}					// kg/m^2
	double getSoilNWL()
	{
		double tmp = soil_nwl_;
		soil_nwl_ = 0.;
		return tmp;
	}	// in kg/m^2
	double getC34Ratio()
	{
		return C34_ratio_;
	}
	double getHeight(int i)
	{
		return Grasses[i].getHeight();
	}
	double getActiveDays()
	{
		return active_days_;
	}

	double getFuelMoisture();
	int getAllDormant();
	double getBlMaxYearC4();			// kg/m^2
	double getBrMaxYearC4();			// kg/m^2
	double getBlMaxYearC3();			// kg/m^2
	double getBrMaxYearC3();			// kg/m^2


	double getC34RatioYearMean();

};

// -----------------------------------------------------------------------------------------------------------

clGrassPop::clGrassPop() :
		Grasses(0, clGrass(0.0, 0))
{
	pop_size_ = 0.;
	root_bm_live_ = 0.;
	stem_bm_live_ = 0.;
	leaf_bm_live_ = 0.;

	root_bm_dead_ = 0.;
	stem_bm_dead_st_ = 0.;
	stem_bm_dead_ly_ = 0.;
	leaf_bm_dead_st_ = 0.;
	leaf_bm_dead_ly_ = 0.;
	leaf_live_comb_ = 0.;
	leaf_dead_st_comb_ = 0.;
	leaf_dead_ly_comb_ = 0.;
	soil_nwl_ = 0.;
	active_days_ = 0.;

	gpp_ = 0.;
	Rma_ = 0.;
	Rgr_ = 0.;

	C34_ratio_ = 0.5;

	cover_[0] = 0.;
	cover_[1] = 0.;
	cover_[2] = 0.;
	cover_[3] = 0.;
	cover_[4] = 0.;
	cover_[5] = 0.;

	for (int i = 0; i < 3650; i++)
	{
		C34_ratio_year_[i] = 0;
	}

}

clGrassPop::~clGrassPop()
{
}

// ------------------------------------------------------------------------------------------------------------
// This method adds grass to the population. We can set the initial size of the 
// trees by the arguments.

void clGrassPop::addGrass(double init_mass, int grass_type)
{
	clGrass Grass(init_mass, grass_type);
	Grasses.push_back(Grass);
	pop_size_++;

	return;
}

// ------------------------------------------------------------------------------------------------------------
// Run physiology for the grass population, calls the physiology for each grass

void clGrassPop::RunPhysiology(double A0_C4, double A0_C3, double RmL_C4, double RmL_C3, double mmsTOkgd, double T, double wind, double *G_theta, double Ca,
		double P, double rh, double mean_sav_height, double mean_for_height, double T_fac, int frost, double p_can_sav, double p_can_for, double *thickness,
		int day)
{
	double dead_bm = 0.;    // stores standing dead grass biomass for the two individuals, needed for ligth competition
	leaf_bm_live_ = 0.;
	stem_bm_live_ = 0.;
	root_bm_live_ = 0.;
	gpp_ = 0.;
	Rma_ = 0.;
	Rgr_ = 0.;

	for (int count = 0; count < pop_size_; count++)
	{
		dead_bm = leaf_bm_dead_st_ * cover_[count];  // passed to physiology, required for light competition

		Grasses[count].RunPhysiology(A0_C4, A0_C3, RmL_C4, RmL_C3, mmsTOkgd, T, wind, G_theta, Ca, P, rh, mean_sav_height, mean_for_height, dead_bm, T_fac,
				frost, thickness, day);

		leaf_bm_live_ += Grasses[count].getBl() * cover_[count];
		root_bm_live_ += Grasses[count].getBr() * cover_[count];

		leaf_bm_dead_st_ += Grasses[count].getBld() * cover_[count];
		root_bm_dead_ += Grasses[count].getBrd() * cover_[count];

		gpp_ += Grasses[count].getGPP() * cover_[count];
		Rma_ += Grasses[count].getRma() * cover_[count];
		Rgr_ += Grasses[count].getRgr() * cover_[count];
	}

	// Transition from standing leaf biomass to lying biomass
	leaf_bm_dead_ly_ += HELPER_BLTL_GRASS * leaf_bm_dead_st_;
	leaf_bm_dead_st_ *= HELPER_BLT__GRASS;

	// Transition of dead biomass into Yasso soil pools
	double tr_leaf_d_nwl = HELPER_LD_NWL_GRASS * leaf_bm_dead_ly_;		// dead leaves to non woody litter
	double tr_root_f_nwl = HELPER_RF_NWL_GRASS * root_bm_dead_;			// roots to non woody litter

	soil_nwl_ += (tr_leaf_d_nwl + tr_root_f_nwl);

	leaf_bm_dead_ly_ -= tr_leaf_d_nwl;
	root_bm_dead_ -= tr_root_f_nwl;

//	double helper_c3 = Grasses[GR_C3_OPN].getGPP();
	double helper_c4 = Grasses[GR_C4_OPN].getGPP();

	if (helper_c4 > 0)
	{
		C34_ratio_ += 0.001 * (Grasses[GR_C3_OPN].getBl() / Grasses[GR_C4_OPN].getBl() - 1.) * (1. - C34_ratio_);
		C34_ratio_ = MyMax(0.01, MyMin(0.99, C34_ratio_));
	}

	for (int i = 0; i < 3649; i++)
	{
		C34_ratio_year_[i] = C34_ratio_year_[i + 1];
	}
	C34_ratio_year_[3649] = C34_ratio_;

	cover_[GR_C4_SAV] = (1. - C34_ratio_) * p_can_sav;
	cover_[GR_C4_OPN] = (1. - C34_ratio_) * (1. - MyMin(1., p_can_sav + p_can_for));
	cover_[GR_C4_FOR] = (1. - C34_ratio_) * p_can_for;
	cover_[GR_C3_SAV] = (C34_ratio_) * p_can_sav;
	cover_[GR_C3_OPN] = (C34_ratio_) * (1. - MyMin(1., p_can_sav + p_can_for));
	cover_[GR_C3_FOR] = (C34_ratio_) * p_can_for;

	if (day == 364)
	{
		//	active_days_ == 0;
		for (int count = 0; count < pop_size_; count++)
			active_days_ += Grasses[count].getActiveDays();
		active_days_ /= (double) pop_size_;
// 		cout << "GSL" << setw(24) << active_days_ << endl;
	}

	return;
}

// ---------------------------------------------------------------------------------------------------------------------------
// This method updates the biomass of the population after a fire. This happens by reducing the leaf and stem
// biomass by fixed part.

void clGrassPop::setStateAfterFire(double patchiness, double cc_fine)
{
	leaf_bm_live_ = 0; //

	leaf_live_comb_ = 0; // this is done in the loop for the individuals
	leaf_dead_st_comb_ = leaf_bm_dead_st_ * patchiness * cc_fine; //
	leaf_dead_ly_comb_ = leaf_bm_dead_ly_ * patchiness * cc_fine; //

	leaf_bm_dead_st_ -= leaf_dead_st_comb_;
	leaf_bm_dead_ly_ -= leaf_dead_ly_comb_;

	combustion_global += leaf_dead_st_comb_ + leaf_dead_ly_comb_;

	// simulate fire effects for grass super-individuals
	for (int count_plant = 0; count_plant < pop_size_; count_plant++)
	{
		Grasses[count_plant].setStateAfterFire(patchiness, cc_fine);

		leaf_bm_live_ += cover_[count_plant] * Grasses[count_plant].getBl();
		leaf_live_comb_ += cover_[count_plant] * Grasses[count_plant].getLeafCombustion();
		combustion_global += leaf_live_comb_;
	}

	return;
}

void clGrassPop::emptyCombustionPools()
{
	leaf_live_comb_ = 0.;
	leaf_dead_st_comb_ = 0.;
	leaf_dead_ly_comb_ = 0.;
	return;
}

// ---------------------------------------------------------------------------------------------------------------------------
// Computes the available biomass for fire by evaluating the mean value
// of sub- and between-canopy-grass. We distinguish between wet (live)
// and dry (dead) biomass

double clGrassPop::getWetBiomassForFire()
{
	return leaf_bm_live_;			// kg/m^2
}

// ---------------------------------------------------------------------------------------------------------------------------

double clGrassPop::getDryBiomassForFire()
{
	return leaf_bm_dead_ly_ + leaf_bm_dead_st_;			// kg/m^2
}

double clGrassPop::getFuelMoisture()
{
	double tmp_moist = 0.;
	for (int i = 0; i < pop_size_; i++)
// 		tmp_moist += (Grasses[i].getFuelMoisture());
		tmp_moist += (Grasses[i].getFuelMoisture() * cover_[i]);

	return tmp_moist;
// 	return tmp_moist/(double)pop_size_;
}

int clGrassPop::getAllDormant()
{
	int tmp = 0;

	for (int i = 0; i < pop_size_; i++)
		tmp += Grasses[i].getDormant();

	return tmp;
}

double clGrassPop::getEt(double T, double P, double radnet_day, double s12, double SP_HEAT, double gama12, double rho12, double VPD12, double p_canopy)
{
	double Et_grasses = 0;

	Et_grasses = Grasses[0].getEt(T, P, radnet_day, s12, SP_HEAT, gama12, rho12, VPD12) * p_canopy
			+ Grasses[1].getEt(T, P, radnet_day, s12, SP_HEAT, gama12, rho12, VPD12) * (1. - p_canopy);

	return Et_grasses;

}

double clGrassPop::getBlMaxYearC4()
{
	return Grasses[0].getBlYearMax() + Grasses[1].getBlYearMax() + Grasses[2].getBlYearMax();	// kg/m^2
}

double clGrassPop::getBrMaxYearC4()
{
	return Grasses[0].getBrYearMax() + Grasses[1].getBrYearMax() + Grasses[2].getBrYearMax();	// kg/m^2
}

double clGrassPop::getBlMaxYearC3()
{
	return Grasses[3].getBlYearMax() + Grasses[4].getBlYearMax() + Grasses[5].getBlYearMax();	// kg/m^2
}

double clGrassPop::getBrMaxYearC3()
{
	return Grasses[3].getBrYearMax() + Grasses[4].getBrYearMax() + Grasses[5].getBrYearMax();	// kg/m^2
}

double clGrassPop::getC34RatioYearMean()
{
	double mean = 0;
	for (int i = 0; i < 3650; i++)
		mean += C34_ratio_year_[i];
	return mean / 3650.;
}

