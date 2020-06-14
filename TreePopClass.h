#pragma once

#include "TreeClass.h"

class clTreePop
{
private:
	int pop_size_;
	int dead_trees_;						// number of trees which died during year
	int new_born_;
	int num_of_seeds_[2];
	int num_of_seeds_collect_[2];
	int count_trees_;
	int wet_days_;
	double root_bm_live_;
	double stem_bm_live_;
	double leaf_bm_live_;
	double root_bm_dead_;
	double stem_bm_dead_st_;
	double stem_bm_dead_ly_;
	double leaf_bm_dead_st_;
	double leaf_bm_dead_ly_;
	double stem_live_comb_;
	double leaf_live_comb_;
	double stem_dead_st_comb_;
	double stem_dead_ly_comb_;
	double leaf_dead_st_comb_;
	double leaf_dead_ly_comb_;
	double gpp_;
	double Rma_;
	double Rgr_;
	double soil_nwl_;							// soil non-woody litter (foilage, fine roots)
	double soil_fwl_;							// soil fine woody litter (branches, coarse roots)
	double soil_cwl_;							// soil coarse woody litter (stems)
	double max_root_depth_;
	double active_days_;

	double canopy_area_;
	double gc_weighted_;
	double mean_height_;
	double basal_area_;
	double max_height_;
	double mean_lai_;
	double evapotransp_;

	double p_canopy_sav_;
	double p_canopy_for_;
	double p_canopy_tot_;

	double Bl_year_[3650];						// stores leaf biomasses of the 10 last years
	double Br_year_[3650];						// stores leaf biomasses of the 10 last years
	double Bs_year_[3650];						// stores leaf biomasses of the 10 last years
	double men_height_year_[3650];				// stores mean heights of the 10 last years
	double max_height_year_[3650];				// stores max. heights of the 10 last years
	double basal_area_year_[3650];				// stores basal areas of the 10 last years
	double p_canopy_sav_year_[3650];			// stores tree cover of the 10 last years, forest+savanna trees
	double p_canopy_for_year_[3650];			// stores tree cover of the 10 last years, forest+savanna trees

	std::vector<clTree> Trees;								// array of trees

	void calpCanopy();

public:
	clTreePop(double max_root_depth);
	~clTreePop();
	void delTree(int number);
	void RunDeathProcess(int frost);
	double getDryBiomassForFire();
	void addTree(double init_mass, int tree_type);
	void addTreeFromSeedbank(int month, double moisture, double field_capacity, double wilt_point, double temperature, int frost);
	void RunPhysiology(double A0, double RmL, double mmsTOkgd, double T, double wind, double *G_theta, double Ca, double P, double rh, double C4_grass_height,
			double C3_grass_height, int day, int month, double soil_moisture, double field_capacity, double wilt_point, double T_fac, int frost,
			double C34_ratio, double *thickness, int soil_layers);
	double getEt(double T, double P, double radnet_day, double s12, double SP_HEAT, double gama12, double rho12, double VPD12);
	void setStateAfterFire(double intensity, double patchiness, double cc_fine, double cc_coarse, double cc_heavy, double cc_tk_helper);
	void emptyCombustionPools();
	void setBornToZero();
	void setDeadToZero();

public:
	int getPopSize()
	{
		return pop_size_;
	}
	int getNewBorn()
	{
		return new_born_;
	}
	int getDeadTrees()
	{
		return dead_trees_;
	}
	int getNumber()
	{
		return count_trees_;
	}
	int getCompTrees();
	int getNumOfSeeds(int i)
	{
		return num_of_seeds_[i];
	}
	int getSmallTrees();
	int getTallTrees();

	double getActiveDays()
	{
		return active_days_;
	}
	double getLeafBmLive()
	{
		return leaf_bm_live_ / GRID_SIZE;
	}			// in kg/m^2
	double getStemBmLive()
	{
		return stem_bm_live_ / GRID_SIZE;
	}			// in kg/m^2
	double getRootBmLive()
	{
		return root_bm_live_ / GRID_SIZE;
	}			// in kg/m^2

	double getLeafBmDeadSt()
	{
		return leaf_bm_dead_st_ / GRID_SIZE;
	}		// in kg/m^2
	double getLeafBmDeadLy()
	{
		return leaf_bm_dead_ly_ / GRID_SIZE;
	}		// in kg/m^2
	double getStemBmDeadSt()
	{
		return stem_bm_dead_st_ / GRID_SIZE;
	}		// in kg/m^2
	double getStemBmDeadLy()
	{
		return stem_bm_dead_ly_ / GRID_SIZE;
	}		// in kg/m^2
	double getRootBmDead()
	{
		return root_bm_dead_ / GRID_SIZE;
	}			// in kg/m^2
	double getStemLiveCombustion()
	{
		return stem_live_comb_ / GRID_SIZE;
	}		// in kg/m^2
	double getLeafLiveCombustion()
	{
		return leaf_live_comb_ / GRID_SIZE;
	}		// in kg/m^2
	double getStemDeadStCombustion()
	{
		return stem_dead_st_comb_ / GRID_SIZE;
	}	// in kg/m^2
	double getStemDeadLyCombustion()
	{
		return stem_dead_ly_comb_ / GRID_SIZE;
	}	// in kg/m^2
	double getLeafDeadStCombustion()
	{
		return leaf_dead_st_comb_ / GRID_SIZE;
	}	// in kg/m^2
	double getLeafDeadLyCombustion()
	{
		return leaf_dead_ly_comb_ / GRID_SIZE;
	}	// in kg/m^2
	double getGPP()
	{
		return gpp_ / GRID_SIZE;
	}					// in kg/m^2
	double getRma()
	{
		return Rma_ / GRID_SIZE;
	}					// in kg/m^2
	double getRgr()
	{
		return Rgr_ / GRID_SIZE;
	}					// in kg/m^2
	double getSoilNWL()
	{
		double tmp = soil_nwl_;
		soil_nwl_ = 0.;
		return tmp / GRID_SIZE;
	}	// in kg/m^2
	double getSoilFWL()
	{
		double tmp = soil_fwl_;
		soil_fwl_ = 0.;
		return tmp / GRID_SIZE;
	}	// in kg/m^2
	double getSoilCWL()
	{
		double tmp = soil_cwl_;
		soil_cwl_ = 0.;
		return tmp / GRID_SIZE;
	}	// in kg/m^2

	double getCanopyArea()
	{
		return canopy_area_;
	}
	double getgcWeighted()
	{
		return gc_weighted_;
	}
	double getMeanHeight()
	{
		return mean_height_ / (double) pop_size_;
	}
	double getBasalArea()
	{
		return basal_area_;
	}
	double getMaxHeight()
	{
		return max_height_;
	}
	double getMeanLai()
	{
		return mean_lai_;
	}
	double getEvapotransp()
	{
		return evapotransp_;
	}
	double getTreeAindex(int number);
	double getTreeHeight(int number);
	double getTreeGc(int number);
	double getTreeGb(int number);
	double getTreeBiomass(int number);
	double getTreeCanopyArea(int number);
	double getTreeDroot(int number);
	double getTreeBl(int number);
	double getTreeBr(int number);
	double getTreeBs(int number);
	double getTreeLAI(int number);
	double getTreenCGT(int number);
	double getTreeGw(int number);
	double getTreeQi(int number);
	double getDormant(int number);
	double getBlYearMax();
	double getBsYearMax();
	double getBrYearMax();
	double getMeanHeightYearMean();
	double getMaxHeightYearMean();
	double getBasalAreaYearMean();
	double getMaxBasalAreaYearMean();
// 		double getpCanopyYearMax();
	double getpCanopySavYearMean();
	double getpCanopyForYearMean();
	double getpCanopy()
	{
		return p_canopy_tot_;
	}
	double getpCanopySav()
	{
		return p_canopy_sav_;
	}
	double getpCanopyFor()
	{
		return p_canopy_for_;
	}

	int getSavTreeNum();
	int getForTreeNum();
	double getMeanSavHeight();
	double getMeanForHeight();

};

// ---------------------------------------------------------------------------------------------------------------------------
clTreePop::~clTreePop()
{
}

// ---------------------------------------------------------------------------------------------------------------------------
// Constructor for a TreePop, initializes some values.
clTreePop::clTreePop(double max_root_depth) :
		Trees()
{
	p_canopy_tot_ = 0.;

	wet_days_ = 0.;
	pop_size_ = 0.;
	dead_trees_ = 0.;
	new_born_ = 0.;
	num_of_seeds_[TR_SAV] = 0;
	num_of_seeds_[TR_FOR] = 0;
	num_of_seeds_collect_[TR_SAV] = 0;
	num_of_seeds_collect_[TR_FOR] = 0;
	count_trees_ = 0.;
	root_bm_live_ = 0.;
	stem_bm_live_ = 0.;
	leaf_bm_live_ = 0.;
	max_root_depth_ = max_root_depth;
	active_days_ = 0.;

	root_bm_dead_ = 0.;
	stem_bm_dead_st_ = 0.;
	stem_bm_dead_ly_ = 0.;
	leaf_bm_dead_st_ = 0.;
	leaf_bm_dead_ly_ = 0.;
	stem_live_comb_ = 0.;
	leaf_live_comb_ = 0.;
	stem_dead_st_comb_ = 0.;
	stem_dead_ly_comb_ = 0.;
	leaf_dead_st_comb_ = 0.;
	leaf_dead_ly_comb_ = 0.;
	soil_nwl_ = 0.;
	soil_fwl_ = 0.;
	soil_cwl_ = 0.;

	gpp_ = 0.;
	Rma_ = 0.;
	Rgr_ = 0.;
	canopy_area_ = 0.;
	gc_weighted_ = 0.;
	mean_height_ = 0.;
	basal_area_ = 0.;
	max_height_ = 0.;
	mean_lai_ = 0.;
	p_canopy_sav_ = 0.;
	p_canopy_for_ = 0.;

	evapotransp_ = 0.;

	for (int i = 0; i < 3650; i++)
	{
		Bl_year_[i] = 0;
		Bs_year_[i] = 0;
		Br_year_[i] = 0;
		men_height_year_[i] = 0;
		max_height_year_[i] = 0;
		basal_area_year_[i] = 0;
		p_canopy_sav_year_[i] = 0;
		p_canopy_for_year_[i] = 0;
	}
}

// ---------------------------------------------------------------------------------------------------------------------------
// This method adds a tree to the population. We can set the initial size of the
// trees by the arguments.
void clTreePop::addTree(double init_mass, int tree_type)
{
	double canopy = -1.;
	if (Trees.size() > 0) canopy = p_canopy_sav_ + p_canopy_for_;
	clTree Tree(init_mass, canopy, pop_size_, count_trees_, tree_type, max_root_depth_);

	Trees.push_back(Tree);
	pop_size_ = Trees.size();
	new_born_++;
	count_trees_++;

	root_bm_live_ += Tree.getBr();
	stem_bm_live_ += Tree.getBs();
	leaf_bm_live_ += Tree.getBl();
	canopy_area_ += Tree.getCanopyArea();
	gc_weighted_ += Tree.getGc();	//*Tree.getCanopyArea();
	basal_area_ += Tree.getStemArea();
	max_height_ = MyMax(max_height_, Tree.getHeight());

	return;
}

// ---------------------------------------------------------------------------------------------------------------------------
// This method delets a tree from the population. First, we add the stem and leaf
// biomass to the dead biomass, then we copy the trees 0..number-1 and
// number+1..pop_size_ to the Trees-array. Then we must update the
// biomass of the population and reduce pop_size_ by one.
void clTreePop::delTree(int number)
{

	leaf_bm_live_ -= Trees[number].getBl();
	stem_bm_live_ -= Trees[number].getBs();
	root_bm_live_ -= Trees[number].getBr();

	if (MyRand() < 0.5)
	{
		leaf_bm_dead_st_ += (Trees[number].getBl() + Trees[number].getBld());
		stem_bm_dead_st_ += (Trees[number].getBs() + Trees[number].getBsd());
	}
	else
	{
		leaf_bm_dead_ly_ += (Trees[number].getBl() + Trees[number].getBld());
		stem_bm_dead_ly_ += (Trees[number].getBs() + Trees[number].getBsd());
	}

	root_bm_dead_ += Trees[number].getBr() + Trees[number].getBrd();

	canopy_area_ -= Trees[number].getCanopyArea();
	basal_area_ -= Trees[number].getStemArea();
	gc_weighted_ -= Trees[number].getGc();	//*Trees[number].getCanopyArea();
	mean_height_ -= Trees[number].getHeight();
	mean_lai_ = ((mean_lai_ * (double) pop_size_) - Trees[number].getLAI()) / (double) (pop_size_ - 1);

	Trees.erase(Trees.begin() + number);

	dead_trees_++;
	pop_size_ = Trees.size();

	return;
}

// ---------------------------------------------------------------------------------------------------------------------------
// Run physiology for the tree population, calls the physiology for each tree
void clTreePop::RunPhysiology(double A0, double RmL, double mmsTOkgd, double T, double wind, double *G_theta, double Ca, double P, double rh,
		double C4_grass_height, double C3_grass_height, int day, int month, double soil_moisture, double field_capacity, double wilt_point, double T_fac,
		int frost, double C34_ratio, double *thickness, int soil_layers)
{
	int comp;
	int comp_type;

	double comp_height;

	double G_theta_weighted[soil_layers];

	for (int layer = 0; layer < soil_layers; layer++)
		G_theta_weighted[layer] = G_theta[layer] * thickness[layer];

	// Empty variable that stores seed number of this year
	// This happens 1/2 year after day where seeds are produced
	if ((int) (A0_max_C3_day + 182) % 365 == day)
	{
		num_of_seeds_[TR_SAV] = num_of_seeds_collect_[TR_SAV];
		num_of_seeds_collect_[TR_SAV] = 0;
		num_of_seeds_[TR_FOR] = num_of_seeds_collect_[TR_FOR];
		num_of_seeds_collect_[TR_FOR] = 0;
		wet_days_ = 0;
	}

	addTreeFromSeedbank(month, soil_moisture, field_capacity, wilt_point, T, frost);

	leaf_bm_live_ = 0.;
	stem_bm_live_ = 0.;
	root_bm_live_ = 0.;
	basal_area_ = 0.;
	canopy_area_ = 0.;
	gc_weighted_ = 0.;
	mean_height_ = 0.;
	max_height_ = 0.;
	mean_lai_ = 0.;
	evapotransp_ = 0.;
	gpp_ = 0.;
	Rma_ = 0.;
	Rgr_ = 0.;

	for (int count = 0; count < pop_size_; count++)
	{
		comp = Trees[count].getCompetitor();   // check if tree has a competitor (comp is index of competitor tree)
		if (comp >= 0 && comp < pop_size_)       // if yes get height and type of competitor tree
		{
			comp_height = Trees[comp].getHeight();
			comp_type = Trees[comp].getTreeType();
		}
		else
		{
			comp_height = 0;
			comp_type = 0;
		}

		Trees[count].RunPhysiology(A0, RmL, mmsTOkgd, T, wind, G_theta_weighted, G_theta[3], Ca, P, rh, comp_height, C4_grass_height, C3_grass_height, day,
				T_fac, frost, comp_type, C34_ratio, thickness);

		leaf_bm_live_ += Trees[count].getBl();
		stem_bm_live_ += Trees[count].getBs();
		root_bm_live_ += Trees[count].getBr();

		leaf_bm_dead_st_ += Trees[count].getBld();
		stem_bm_dead_st_ += Trees[count].getBsd();
		root_bm_dead_ += Trees[count].getBrd();

		gpp_ += Trees[count].getGPP();
		Rma_ += Trees[count].getRma();
		Rgr_ += Trees[count].getRgr();

		basal_area_ += Trees[count].getStemArea();
		canopy_area_ += Trees[count].getCanopyArea();
		gc_weighted_ += Trees[count].getGc();       //*Trees[count].getCanopyArea();
		mean_height_ += Trees[count].getHeight();
		max_height_ = MyMax(max_height_, Trees[count].getHeight());
		mean_lai_ += Trees[count].getLAI();

		if (day == A0_max_C3_day) num_of_seeds_collect_[Trees[count].getTreeType()] += Trees[count].getSeedProduction();

		if (Trees[count].getRootSucker() == 1 && day > (A0_max_C3_day - 1) && day < (A0_max_C3_day + 1)) addTree(1., Trees[count].getTreeType());
	}

	num_of_seeds_collect_[TR_SAV] = (int) ceil(SEED_DECAY_RATE[TR_SAV] * (double) num_of_seeds_collect_[TR_SAV]);
	num_of_seeds_collect_[TR_FOR] = (int) ceil(SEED_DECAY_RATE[TR_FOR] * (double) num_of_seeds_collect_[TR_FOR]);

	num_of_seeds_[TR_SAV] = (int) ceil(SEED_DECAY_RATE[TR_SAV] * (double) num_of_seeds_[TR_SAV]);
	num_of_seeds_[TR_FOR] = (int) ceil(SEED_DECAY_RATE[TR_FOR] * (double) num_of_seeds_[TR_FOR]);

	mean_lai_ /= (double) pop_size_;

	calpCanopy();

	for (int i = 0; i < 3649; i++)
	{
		Bl_year_[i] = Bl_year_[i + 1];
		Bs_year_[i] = Bs_year_[i + 1];
		Br_year_[i] = Br_year_[i + 1];
		men_height_year_[i] = men_height_year_[i + 1];
		max_height_year_[i] = max_height_year_[i + 1];
		basal_area_year_[i] = basal_area_year_[i + 1];
		p_canopy_sav_year_[i] = p_canopy_sav_year_[i + 1];
		p_canopy_for_year_[i] = p_canopy_for_year_[i + 1];
	}

	Bl_year_[3649] = leaf_bm_live_;
	Bs_year_[3649] = stem_bm_live_;
	Br_year_[3649] = root_bm_live_;
	men_height_year_[3649] = mean_height_ / (int) pop_size_;
	max_height_year_[3649] = max_height_;
	basal_area_year_[3649] = basal_area_;
	p_canopy_sav_year_[3649] = p_canopy_sav_;
	p_canopy_for_year_[3649] = p_canopy_for_;

	// Transition from hanging/standing leaf/stem biomass to lying biomass
	leaf_bm_dead_ly_ += HELPER_BLTL_TREE * leaf_bm_dead_st_;
	leaf_bm_dead_st_ *= HELPER_BLT__TREE;

	stem_bm_dead_ly_ += HELPER_BSTL_TREE * stem_bm_dead_st_;
	stem_bm_dead_st_ *= HELPER_BST__TREE;

	// Transition from lying biomass pools into soil pools
	double bl_bs_ratio = leaf_bm_live_ / stem_bm_live_;
	double tr_leaf_d_nwl = HELPER_LD_NWL_TREE * leaf_bm_dead_ly_;					// dead leaves to non woody litter
	double tr_root_f_nwl = HELPER_RF_NWL_TREE * root_bm_dead_ * bl_bs_ratio;	// fine roots to non woody litter
	double tr_root_c_fwl = HELPER_RC_FWL_TREE * root_bm_dead_ * (1. - bl_bs_ratio);	// coarse roots to fine woody litter
	double tr_stem_f_fwl = HELPER_SF_FWL_TREE * stem_bm_dead_ly_ * SFRAC_FINE;		// fine stem biomass to fine woody litter
	double tr_stem_c_cwl = HELPER_SC_CWL_TREE * stem_bm_dead_ly_ * (1. - SFRAC_FINE);	// coarse/heavy stem biomas to coarse woody litter

	soil_nwl_ += (tr_leaf_d_nwl + tr_root_f_nwl);
	soil_fwl_ += (tr_root_c_fwl + tr_stem_f_fwl);
	soil_cwl_ += (tr_stem_c_cwl);

	leaf_bm_dead_ly_ -= tr_leaf_d_nwl;
	root_bm_dead_ -= tr_root_f_nwl;
	root_bm_dead_ -= tr_root_c_fwl;
	stem_bm_dead_ly_ -= tr_stem_f_fwl;
	stem_bm_dead_ly_ -= tr_stem_c_cwl;


	if (day == 364)
	{
//		active_days_ == 0;  //We are looking at cumulative active days
		for (int count = 0; count < pop_size_; count++)
			active_days_ += Trees[count].getActiveDays();
		active_days_ /= (double) pop_size_;
// 		cout << "GSL" << setw(14) << active_days_ << endl;
	}

	return;
}

// ---------------------------------------------------------------------------------------------------------------------------
// Compute survival probability for each tree and delete tree from the population if it is too small.
// Tree 0 cannot be deleted
void clTreePop::RunDeathProcess(int frost)
{
	int count = pop_size_ - 1;

	while (count > 1)
	{
		if (Trees[count].WillIDie(frost) == 1)
		{
			delTree(count);
		}
		count--;
	}
	//tree 0+1 is special case as it is immortal and cannot be removed from population
	Trees[0].WillIDie(frost);
	Trees[1].WillIDie(frost);
	calpCanopy();

	return;
}

// ---------------------------------------------------------------------------------------------------------------------------
void clTreePop::setStateAfterFire(double intensity, double patchiness, double cc_fine, double cc_coarse, double cc_heavy, double cc_tk_helper)
{

	canopy_area_ = 0; //
	basal_area_ = 0; //
	gc_weighted_ = 0; //

	leaf_bm_live_ = 0; //
	stem_bm_live_ = 0; //

	stem_live_comb_ = 0; //
	leaf_live_comb_ = 0; //
	stem_dead_st_comb_ = 0; // NOTE what happens with this pool?
	stem_dead_ly_comb_ = stem_bm_dead_ly_ * patchiness * (cc_fine * SFRAC_FINE + cc_coarse * SFRAC_COARSE + cc_heavy * SFRAC_HEAVY); //
	leaf_dead_st_comb_ = leaf_bm_dead_st_ * patchiness * cc_fine; //
	leaf_dead_ly_comb_ = leaf_bm_dead_ly_ * patchiness * cc_fine; //

	stem_bm_dead_st_ -= stem_dead_st_comb_;
	stem_bm_dead_ly_ -= stem_dead_ly_comb_;
	leaf_bm_dead_st_ -= leaf_dead_st_comb_;
	leaf_bm_dead_ly_ -= leaf_dead_ly_comb_;

	combustion_global += stem_dead_st_comb_ + stem_dead_ly_comb_ + leaf_dead_st_comb_ + leaf_dead_ly_comb_;

	// simulate fire effects for single trees and update live biomass pools of population
	for (int count_plant = 0; count_plant < pop_size_; count_plant++)
	{
		Trees[count_plant].setStateAfterFire(intensity, patchiness, cc_fine, cc_tk_helper);

		stem_bm_live_ += Trees[count_plant].getBs();
		leaf_bm_live_ += Trees[count_plant].getBl();
		stem_bm_dead_st_ += Trees[count_plant].getBsd();   // stems from topkilled trees move to dead standing stem pool
		leaf_bm_dead_st_ += Trees[count_plant].getBld();   // nothing happens to this pool
		stem_live_comb_ += Trees[count_plant].getStemCombustion();
		leaf_live_comb_ += Trees[count_plant].getLeafCombustion();
		basal_area_ += Trees[count_plant].getStemArea();
		canopy_area_ += Trees[count_plant].getCanopyArea();
		gc_weighted_ += Trees[count_plant].getGc();   //*Trees[count_plant].getCanopyArea();
	}

	int count = pop_size_ - 1;

	while (count > 1)
	{
		if (Trees[count].WillIDieAfterFire() == 1)
		{
			delTree(count);
		}
		count--;
	}

	// tree 0+1 is special case as it is immortal and cannot be removed from population
	Trees[0].WillIDieAfterFire();
	Trees[1].WillIDieAfterFire();
	calpCanopy();

	return;
}

void clTreePop::emptyCombustionPools()
{
	stem_live_comb_ = 0.;
	leaf_live_comb_ = 0.;
	stem_dead_st_comb_ = 0.;
	stem_dead_ly_comb_ = 0.;
	leaf_dead_st_comb_ = 0.;
	leaf_dead_ly_comb_ = 0.;
	return;
}

// ---------------------------------------------------------------------------------------------------------------------------
// The intensity of fire is dependent on dry tree biomass per m^2
double clTreePop::getDryBiomassForFire()
{  // NOTE: needs update
	return leaf_bm_dead_ly_ / GRID_SIZE + stem_bm_dead_ly_ * SFRAC_FINE / GRID_SIZE;   // NOTE /GRID_SIZE for kg/m^2
}

// ---------------------------------------------------------------------------------------------------------------------------
void clTreePop::setBornToZero()
{
	new_born_ = 0;
}

void clTreePop::setDeadToZero()
{
	dead_trees_ = 0;
}

// ---------------------------------------------------------------------------------------------------------------------------
void clTreePop::addTreeFromSeedbank(int month, double moisture, double field_capacity, double wilt_point, double temperature, int frost)
{
	if (moisture >= field_capacity && frost == 0) wet_days_++;
	if (moisture < field_capacity || frost == 1) wet_days_--;
	if (wet_days_ < 0) wet_days_ = 0;

	if (wet_days_ >= 3)
	{
		wet_days_ -= 1;

		for (int type = TR_SAV; type <= TR_FOR; type++)
		{
			if (num_of_seeds_[type] > 0)
			{
				int seeds = (int) ceil((double) num_of_seeds_[type] * SEED_FRAC_DAY[type]);
				for (int count = 0; count < seeds; count++)
				{
					if (MyRand() < SEED_GERM_PROB[type])
					{
						addTree(INIT_MASS_TREE, type);
					}
					num_of_seeds_[type]--;
				}
			}
		}
	}

	return;
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeHeight(int number)
{
	return Trees[number].getHeight();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeAindex(int number)
{
	return Trees[number].getAindex();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeGc(int number)
{
	return Trees[number].getGc();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeGb(int number)
{
	return Trees[number].getGb();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeBiomass(int number)
{
	return Trees[number].getBl() + Trees[number].getBr() + Trees[number].getBs();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeCanopyArea(int number)
{
	return Trees[number].getCanopyArea();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeDroot(int number)
{
	return Trees[number].getDroot();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeBl(int number)
{
	return Trees[number].getBl();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeBr(int number)
{
	return Trees[number].getBr();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeBs(int number)
{
	return Trees[number].getBs();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeLAI(int number)
{
	return Trees[number].getLAI();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreenCGT(int number)
{
	return Trees[number].getnCGT();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeGw(int number)
{
	return Trees[number].getGw();
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getTreeQi(int number)
{
	return Trees[number].getQi();
}

double clTreePop::getDormant(int number)
{
	return Trees[number].getDormant();
}

// ---------------------------------------------------------------------------------------------------------------------------
int clTreePop::getCompTrees()
{
	int wt = 0;

	for (int i = 0; i < pop_size_; i++)
		if (Trees[i].getCompetitor() == -1) wt++;

	return pop_size_ - wt;
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getEt(double T, double P, double radnet_day, double s12, double SP_HEAT, double gama12, double rho12, double VPD12)
{
	double Et_trees = 0;

	for (int cnt = 0; cnt < pop_size_; cnt++)
	{
		Et_trees += Trees[cnt].getEt(T, P, radnet_day, s12, SP_HEAT, gama12, rho12, VPD12);
	}

	return Et_trees / GRID_SIZE;  // in m/s
}

// ---------------------------------------------------------------------------------------------------------------------------
int clTreePop::getSmallTrees()
{
	int small_trees = 0;

	for (int cnt = 0; cnt < pop_size_; cnt++)
		if (Trees[cnt].getHeight() <= 1.) small_trees++;

	return small_trees;
}

// ---------------------------------------------------------------------------------------------------------------------------
int clTreePop::getTallTrees()
{
	int tall_trees = 0;

	for (int cnt = 0; cnt < pop_size_; cnt++)
		if (Trees[cnt].getHeight() >= 5.) tall_trees++;

	return tall_trees;
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getBlYearMax()
{
	double max = 0;
	for (int i = 0; i < 3650; i++)
	{
		if (Bl_year_[i] > max) max = Bl_year_[i];
	}

	return max / GRID_SIZE;  // in kg/m^2
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getBsYearMax()
{
	double max = 0;
	for (int i = 0; i < 3650; i++)
	{
		if (Bs_year_[i] > max) max = Bs_year_[i];
	}

	return max / GRID_SIZE;  // in kg/m^2
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getBrYearMax()
{
	double max = 0;
	for (int i = 0; i < 3650; i++)
	{
		if (Br_year_[i] > max) max = Br_year_[i];
	}

	return max / GRID_SIZE;  // in kg/m^2
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getMeanHeightYearMean()
{
	double mean = 0;
	for (int i = 0; i < 3650; i++)
		mean += men_height_year_[i];
	return mean / 3650.;
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getMaxHeightYearMean()
{
	double mean = 0;
	for (int i = 0; i < 3650; i++)
		mean += max_height_year_[i];
	return mean / 3650.;
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getBasalAreaYearMean()
{
	double mean = 0;
	for (int i = 0; i < 3650; i++)
		mean += basal_area_year_[i];
	return mean / 3650.;
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getMaxBasalAreaYearMean()
{
	double max = 0;
	for (int i = 0; i < 3650; i++)
		max = MyMax(basal_area_year_[i], max);
	return max;
}

// ---------------------------------------------------------------------------------------------------------------------------
double clTreePop::getpCanopySavYearMean()
{
	double mean = 0;
	for (int i = 0; i < 3650; i++)
		mean += p_canopy_sav_year_[i];
	return mean / 3650.;
}

double clTreePop::getpCanopyForYearMean()
{
	double mean = 0;
	for (int i = 0; i < 3650; i++)
		mean += p_canopy_for_year_[i];
	return mean / 3650.;
}

// ---------------------------------------------------------------------------------------------------------------------------
// new functions for savanna and forest trees
int clTreePop::getSavTreeNum()
{
	int cnt = 0;
	for (int i = 0; i < pop_size_; i++)
		cnt += Trees[i].getTreeType();

	return pop_size_ - cnt;
}

int clTreePop::getForTreeNum()
{
	int cnt = 0;
	for (int i = 0; i < pop_size_; i++)
		cnt += Trees[i].getTreeType();

	return cnt;
}

double clTreePop::getMeanSavHeight()
{
	int cnt = 0;
	int sum = 0;
	for (int i = 0; i < pop_size_; i++)
	{
		if (Trees[i].getTreeType() == TR_SAV)
		{
			cnt++;
			sum += Trees[i].getHeight();
		}
	}

	return sum / (double) cnt;
}

double clTreePop::getMeanForHeight()
{
	int cnt = 0;
	int sum = 0;
	for (int i = 0; i < pop_size_; i++)
	{
		if (Trees[i].getTreeType() == TR_FOR)
		{
			cnt++;
			sum += Trees[i].getHeight();
		}
	}

	return sum / (double) cnt;
}

// ---------------------------------------------------------------------------------------------------------------------------
// new functions for savanna and forest trees

void clTreePop::calpCanopy()
{
	int i;
	double caSumTot = 0;
	double caSumSav = 0;
	double caSumFor = 0;

	for (i = 0; i < pop_size_; i++)
	{
		if (Trees[i].getCompetitor() == -1 && Trees[i].getHeight() > 0.5)
		{
			caSumTot += Trees[i].getCanopyArea();
			if (Trees[i].getTreeType() == TR_SAV) caSumSav += Trees[i].getCanopyArea();
			if (Trees[i].getTreeType() == TR_FOR) caSumFor += Trees[i].getCanopyArea();
		}
	}

	caSumTot /= (double) GRID_SIZE;
	caSumSav /= (double) GRID_SIZE;
	caSumFor /= (double) GRID_SIZE;

	p_canopy_tot_ = MyMin(1., caSumTot);
	p_canopy_sav_ = MyMin(1., caSumSav);
	p_canopy_for_ = MyMin(1., caSumFor);

	if (p_canopy_tot_ == 1.)
	{
		p_canopy_sav_ = caSumSav / caSumTot;
		p_canopy_for_ = caSumFor / caSumTot;
	}
}

