#pragma once

#include "Globals.h"
#include "Fxn.h"

using namespace std;

class clTree
{
private:
	int Dpos_; 								// dormancy counter (days where nCGT>0)
	int Dneg_; 								// standby counter (days where nCGT=0)
	int Dormant_; 							// dormant yes =1, dormant no=0
	int competitor_;						// index of tree for lightcompetition
	int number_;							// number of tree
	int age_;								// age of tree in days
	int reprod_strategy_;					// reproduction strategy, 0=seeds, 1=root sucker
	int tree_type_;							// 0=savanna tree, 1=forest tree
	int active_days_;						// counts active days per year

	double Bl_;    							// leaf biomass, kg per plant
	double Br_;    							// root biomass, kg per plant
	double Bs_;    							// stem biomass, kg per plant

	double Bld_;								// dead leaf biomass, kg per plant
	double Bsd_;								// dead stem biomass, kg per plant
	double Brd_;								// dead root biomass, kg per plant

	double leaf_combustion_;
	double stem_combustion_;

	double gpp_;
	double Rma_;
	double Rgr_;

	double Qsum_;
	double Qi_;
	double Ci_;
	double gc_;    							// canopy stomatal conductance (mu mol/m2/s)   molar
	double gb_;								// canopy boundary layer conductanced (m/s)    non-molar
	double death_prob_;						// probability, that this tree will die
	double Droot_;   							// rooting depth m
	double nCGT_;
	double Gw_;
	double alloc_denom_;
	double Aindex_;
	double height_;							// height m
	double canopy_radius_;
	double canopy_area_;						// canopy area m^2
	double LAI_;	 							// leaf area index of this tree
	double site_param_;						// "quality" of site where plant growth
	double stem_area_;							// in m^2
	double max_root_depth_;					// maximum rooting depth

private:
	void calDroot();
	void calPlantHeight();
	void calCanopyRadius();
	void calCanopyArea();
	void calStemArea();
	void calPlantLAI();
	void calSoilMoisture(double *G_theta_weighted, double *thickness);
	void calGc(double A, double resp, double Ca, double P, double hs, double T);
	void calGb(double wind);
	void AllocateCorbon();
	void UpdateAllometry();
	void MassTurnover();
	void Respirate(double rRoot, double rStem, double rLeaf);
	double AllocLeaf();
	double AllocStem();
	double AllocRoot();
	double TopKillProb(double intensity);
	double calQsum();

public:
	clTree(double init_mass, double pCanopy, int popsize, int number, int tree_type_, double max_root_depth);
	~clTree();
	int WillIDie(int frost);
	int WillIDieAfterFire();
	void setStateAfterElephants();
	void setStateAfterFire(double intensity, double patchiness, double cc_fine, double cc_tk_helper);
	void RunPhysiology(double A0, double RmL, double mmsTOkgd, double T, double wind, double *G_theta_weighted, double G_theta_3, double Ca, double P,
			double rh, double comp_height, double C4_grass_height, double C3_grass_height, int day, double T_fac, int frost, int comp_type, double C34_ratio,
			double *thickness);
	double getEt(double T, double P, double radnet_day, double s12, double SP_HEAT, double gama12, double rho12, double VPD12);
	int getSeedProduction();
	int getRootSucker();

	int getCompetitor()
	{
		return competitor_;
	}
	int getNumber()
	{
		return number_;
	}
	int getReprodStrategy()
	{
		return reprod_strategy_;
	}
	int getDormant()
	{
		return Dormant_;
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
	double getCi()
	{
		return Ci_;
	}
	double getGw()
	{
		return Gw_;
	}
	double getGc()
	{
		return gc_;
	}
	double getGb()
	{
		return gb_;
	}
	double getDeathProb()
	{
		return death_prob_;
	}
	double getDroot()
	{
		return Droot_;
	}
	double getnCGT()
	{
		return nCGT_;
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
	double getStemArea()
	{
		return stem_area_;
	}
	double getLAI()
	{
		return LAI_;
	}
	int getTreeType()
	{
		return tree_type_;
	}
	int getActiveDays()
	{
		return active_days_;
	}

};

// -------------------------------------------------------------------------------------------------------------------
// This constructor generates a tree of random size between 0 and max_Bl, whereas Bl is the leaf biomass.
// Bs and Br are in a constant ratio to Bl. 
clTree::clTree(double init_mass, double pCanopy, int popsize, int number, int tree_type, double max_root_depth)
{
	number_ = number;
	age_ = 365 * (int) init_mass;
	Dpos_ = 0;
	Dneg_ = 0;
	Dormant_ = 1;
	competitor_ = -1;

	Bl_ = init_mass * 0.01;
	Br_ = init_mass * 0.5;
	Bs_ = init_mass * 0.49;
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
	Ci_ = 1.;
	gc_ = 0;
	gb_ = 0;
	death_prob_ = 0;
	nCGT_ = 0.;
	Gw_ = 1.;
	Aindex_ = .2;
	canopy_area_ = 0.01;
	stem_area_ = 0.00001;
	site_param_ = 0;
	tree_type_ = tree_type;
	active_days_ = 0;
	max_root_depth_ = max_root_depth;

	reprod_strategy_ = 0;
	if (MyRand() < PROB_ROOT_SUCKER[tree_type_]) reprod_strategy_ = 1;

	// compute some allometric measures, functions in Fxn
	UpdateAllometry();

	// we determine a neighbour tree for light competition.
	double has_comp = MyRand();
	if (has_comp > COMP_PAR_1 - COMP_PAR_2 * pCanopy - popsize * COMP_PAR_3) competitor_ = (int) floor(popsize * MyRand());
}

// -------------------------------------------------------------------------------------------------------------------
clTree::~clTree()
{
}

void clTree::calDroot()
{
	Droot_ = CylinderDepth(Br_, MIN_ROOT_RAD_TREE, ROOT_DENSITY_TREE, max_root_depth_);
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// Standard power function, based on fits to Niklas&Enquist data, estimates the plant height from the stem biomass.
void clTree::calPlantHeight()
{
	height_ = HEIGHT_C1_TREE[tree_type_] * pow(Bs_ + Bl_, HEIGHT_C2_TREE[tree_type_]);
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// Compute radius of the canopy
void clTree::calCanopyRadius()
{
	canopy_radius_ = height_ * GAMMA_CANOPY_TREE[tree_type_];
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// Compute the canopy area of a tree,
// assumes canopy radius is constant proportion of height
void clTree::calCanopyArea()
{
	canopy_area_ = M_PI * canopy_radius_ * canopy_radius_;

	return;
}

// -------------------------------------------------------------------------------------------------------------------
void clTree::calStemArea()
{
	stem_area_ = STEM_AREA_HELPER[tree_type_] * height_;		// from Higgins 2007 Ecology; 200 to get radius in m
	stem_area_ = stem_area_ * stem_area_ * 3.141593;		//  from Higgins 2007 Ecology;
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// compute leaf area index
void clTree::calPlantLAI()
{
	LAI_ = Bl_ * SLA_TREE[tree_type_] / MyMax(0.1, canopy_area_) + MIN_LAI;
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// compute the allometric values of the tree
void clTree::UpdateAllometry()
{
	calDroot();
	calPlantHeight();
	calCanopyRadius();
	calCanopyArea();
	calStemArea();
	calPlantLAI();

	return;
}

// ---------------------------------------------------------------------------------------------------------------
int clTree::WillIDie(int frost)
{
	int will_I = 0;

	death_prob_ = DEATH_PROB_FROST[tree_type_] * (double) frost;
	if (Dormant_ == 0 && nCGT_ < 0.) death_prob_ += DEATH_PROB_CARBON[tree_type_];
	if (competitor_ > 0) death_prob_ += DEATH_PROB_COMP[tree_type_];

	if (MyRand() < death_prob_) will_I = 1;

	if (number_ < 2 && will_I == 1)
	{
		Bld_ += 0.99 * Bl_;
		Bsd_ += 0.99 * Bs_;
		Brd_ += 0.99 * Br_;
		Bl_ *= 0.01;
		Bs_ *= 0.01;
		Br_ *= 0.01;
	}

	return will_I;
}

int clTree::WillIDieAfterFire()
{
	int will_I = 0;

	if (MyRand() > RESPROUTING_PROB[tree_type_]) will_I = 1;

#	ifdef S_NOMORT
	will_I = 0;
	#	endif

	if (number_ < 2 && will_I == 1)
	{
		Bld_ += 0.99 * Bl_;
		Bsd_ += 0.99 * Bs_;
		Brd_ += 0.99 * Br_;
		Bl_ *= 0.01;
		Bs_ *= 0.01;
		Br_ *= 0.01;
	}

	return will_I;
}

// -------------------------------------------------------------------------------------------------------------------
// allocation function of net carbon gain to the roots
double clTree::AllocRoot()
{
	return (A0_ROOT_TREE_HELPER[tree_type_] - Gw_) * alloc_denom_;
}

// -------------------------------------------------------------------------------------------------------------------
// allocation function of net carbon gain to the stem
double clTree::AllocStem()
{
	return (A0_STEM_TREE_HELPER[tree_type_] - Qi_) * alloc_denom_;
}

// -------------------------------------------------------------------------------------------------------------------
// allocation function of net carbon gain to the leaf
double clTree::AllocLeaf()
{
	return (1. - Ci_) * alloc_denom_;
}

// -------------------------------------------------------------------------------------------------------------------
// call the alloaction functions for leaf, stem and root, and add it to the tree biomass
void clTree::AllocateCorbon()
{
	if (Dormant_ == 0)
	{
		double npp = gpp_ - Rgr_;
		alloc_denom_ = 1. / (ALLOC_DENOM_HELPER[tree_type_] - Qi_ - Gw_ - Ci_);
		Bl_ += (npp * AllocLeaf());
		Bs_ += (npp * AllocStem());
		Br_ += (npp * AllocRoot());
	}

	return;
}

// -------------------------------------------------------------------------------------------------------------------
// run the tree physiology
void clTree::RunPhysiology(double A0, double RmL, double mmsTOkgd, double T, double wind, double *G_theta_weighted, double G_theta_3, double Ca, double P,
		double rh, double comp_height, double C4_grass_height, double C3_grass_height, int day, double T_fac, int frost, int comp_type, double C34_ratio,
		double *thickness)
{
	double Acs;			// canopy and canopy-stressed photosythesis
	double RmS;
	double RmR;			// maint respiration total, stem, root
	double RmLc;			// canopy leaf maint resp
	double RmLg;			// leaf maint resp in kg/day/plant

	gpp_ = 0.;
	Rma_ = 0.;
	Rgr_ = 0.;

	if (day == 0)
	{
		active_days_ = 0;
	}

	// water availability
	calSoilMoisture(G_theta_weighted, thickness);

	// light availability
	Qi_ = 1.;

	if (MyRand() > C34_ratio)
	{          // trees must compete with C3 grasses or C4 grasses
		if (C4_grass_height > height_)    // C4 grass competition
		Qi_ *= LICMP_1[GR_C4_OPN + 2][tree_type_] / C4_grass_height * height_ + LICMP_2[GR_C4_OPN + 2][tree_type_];
	}
	else
	{
		if (C3_grass_height > height_)    // C3 grass competition
		Qi_ *= LICMP_1[GR_C3_OPN + 2][tree_type_] / C3_grass_height * height_ + LICMP_2[GR_C3_OPN + 2][tree_type_];
	}

	if (comp_height > height_)           // tree-tree competition
	Qi_ *= LICMP_1[comp_type][tree_type_] / comp_height * height_ + LICMP_2[comp_type][tree_type_];

	Qsum_ = Qi_ * calQsum();
	Ci_ = MyMin(1., Bl_ / (Bl_ + Br_ + Bs_) * A0_LEAF_TREE_INV[tree_type_]);

	// conopy photosynthesis
	Acs = A0 * Qsum_ * Gw_;		// water and light stressed canopy photosynthesis (=A0*Qsum*Gw) (mmol/m^2/s)
	// carbon gain (gross primary production)
	gpp_ = Acs * mmsTOkgd * (double) (1 - Dormant_);  // water and light stressed carbon gain, kg/day/m^2, gpp=0 if plant is dormant
	gpp_ = MyMax(gpp_, gpp_ * canopy_area_);   // transform to kg/day/plant

	// components of maintenance respiration
	double bl_bs_ratio = Bl_ / Bs_;
	RmLc = R_MAINT_RESP_TR[tree_type_] * Acs; //maint resp leaf is now function of water stress and canopy photosynthesis (mmol/m^2/s)
	RmLg = RmLc * mmsTOkgd;                   //leaf maint resp in kg/day/plant
	RmLg = MyMax(RmLg, RmLg * canopy_area_);
	RmS = MaintRespFast(BETA_STEM_TREE * Bs_, BETA_N, UPSILON_STEM_TREE, T_fac);
	RmR = MaintRespFast(BETA_ROOT_TREE * (Br_ * (1. - bl_bs_ratio)), BETA_N, UPSILON_ROOT_TREE, T_fac)   // coarse roots
	+ Br_ * bl_bs_ratio * RmLg;                                   // fine roots
// 	cout << RmLg/Bl_ << endl;
	Rgr_ = gpp_ * SIGMA_GROW_RESP_TREE[tree_type_];
	nCGT_ = gpp_ - Rgr_ - RmR - RmS - RmLg;

	//canopy conductance
	calGb(wind);  // m/s
	calGc(Acs, RmLc, Ca, P, rh, T);  // (mmol/m^2/s)
	double Ti;
	if (tmp_min_month > TI_CONST) Ti = 0.;
	else
		Ti = 1.5 * (-1. + tmp_min_month / TI_CONST);

	Aindex_ = A0 * ((WATER_IND_TREE[tree_type_] * G_theta_3 + WATER_IND_TREE_1[tree_type_] * Gw_) + Ti) - RmL + site_param_;

	// if no carbon gain increment days without carbon gain
	if (Dormant_ == 0 && Aindex_ <= STRESS_INDEX_TREE[tree_type_]) Dneg_++;
	// if carbon gain then reset days without carbon gain counter
	if (Dormant_ == 0 && Aindex_ > STRESS_INDEX_TREE[tree_type_]) Dneg_ = 0;
	// if dormant and carbon gain then increment dormancy break counter
	if (Dormant_ == 1 && Aindex_ > STRESS_INDEX_TREE[tree_type_]) Dpos_++;
	// if dormant and carbon gain = 0 then reset dormancy break counter
	if (Dormant_ == 1 && Aindex_ <= STRESS_INDEX_TREE[tree_type_]) Dpos_ = 0;

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
	if (Dneg_ >= D_NEG_TREE[tree_type_])
	{
		if (number_ == 0) phen_counter++;
		Dormant_ = 1;
		Dneg_ = 0;
		Bld_ += (1. - REM_BM_TREE) * Bl_; // leaf abscission NOTE new version: permament turnover and litterfall
		Bl_ = REM_BM_TREE * Bl_;      // leaf abscission NOTE new version: permament turnover and litterfall
	}

	// if days of positive nCGT > threshold then move to metabolic state
	if (Dpos_ >= D_POS_TREE[tree_type_])
	{
		Dormant_ = 0;
		Dpos_ = 0;

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

	age_++;

	return;
}

// ---------------------------------------------------------------------------------------------------------------------------
// compute the water aviability of a tree. The value is between 0 and 1 and is influenced by the
// soil moisture and the rooting depth
// FASTER VERSION: G_theta*THICK already calculated in TreePop, this vector is equal vor all trees!
void clTree::calSoilMoisture(double *G_theta_weighted, double *thickness)
{
	int layers = 1;    // counts the number of layers, in which the plant has roots
	double depth = thickness[0];
	Gw_ = G_theta_weighted[0];

	while (Droot_ > depth)
	{
		Gw_ += G_theta_weighted[layers];
		depth += thickness[layers];
		layers++;
	}

	Gw_ /= depth;

	return;
}

// --------------------------------------------------------------------------------------------------------------
void clTree::MassTurnover()
{
//  new model for turnover: fine root to coarse root ratio equals leaf to stem ratio
//                          fine root respiration equals leaf respiration
//                          coarse root respiration equals stem respiration
	double bl_bs_ratio = Bl_ / Bs_;
	double f_root_turnover = bl_bs_ratio * Br_ * HELPER_BLD_TREE[tree_type_];
	double c_root_turnover = (1. - bl_bs_ratio) * Br_ * HELPER_BRD_TREE;

	if (Dormant_ == 0)
	{
		Bld_ += HELPER_BLD_TREE[tree_type_] * Bl_;
		Bl_ *= HELPER_BL__TREE[tree_type_];
	}

	Bsd_ += HELPER_BSD_TREE * Bs_;
	Bs_ *= HELPER_BS__TREE;

	Brd_ += (f_root_turnover + c_root_turnover);
	Br_ -= (f_root_turnover + c_root_turnover);

	return;
}

// --------------------------------------------------------------------------------------------------------------
void clTree::Respirate(double rRoot, double rStem, double rLeaf)
{
	double tmp_bm = Br_ + Bs_ + Bl_;

	Br_ = MyMax(0.005, Br_ - rRoot);
	Bs_ = MyMax(0.005, Bs_ - rStem);
	Bl_ = MyMax(0.005, Bl_ - rLeaf);

	Rma_ = tmp_bm - (Br_ + Bs_ + Bl_);

	return;
}

// --------------------------------------------------------------------------------------------------------------
void clTree::setStateAfterFire(double intensity, double patchiness, double cc_fine, double cc_tk_helper)
{
	leaf_combustion_ = 0.;
	stem_combustion_ = 0.;

	if (MyRand() < patchiness)
	{
		double omega_top_kill;
		omega_top_kill = TopKillProb(intensity);

		if (MyRand() < omega_top_kill)
		{
			double bl_rem = MyMax(0.001 * Bl_, 0.005);   // remaining live leaf biomass
			double bs_rem = MyMax(0.01 * Bs_, 0.005);   // remaining live stem biomass

			leaf_combustion_ = Bl_ * cc_fine;  // combustion
			stem_combustion_ = Bs_ * CombComplTopkillTree(cc_tk_helper, height_);  // combustion

			combustion_global += stem_combustion_ + leaf_combustion_;

			Bsd_ += Bs_ - stem_combustion_ - bs_rem;
			Bs_ = bs_rem;

			Bld_ += Bl_ - leaf_combustion_ - bl_rem;
			Bl_ = bl_rem;

			Dpos_ = 0;
			Dneg_ = 0;                           // NOTE need to simulate resprouting probability!!
			UpdateAllometry();
		}
	}

	return;
}

// --------------------------------------------------------------------------------------------------------------

void clTree::setStateAfterElephants()
{
	Bl_ = 0.05 * Bl_;
	Bs_ = 0.05 * Bs_;
	Dpos_ = 0;
	Dneg_ = 0;
	UpdateAllometry();

	return;
}

// ---------------------------------------------------------------------------------------------------------------------------
// Function to compute the topkill probability for a tree
double clTree::TopKillProb(double intensity)
{
	double prob = 0.;

	if (height_ > 0) prob = (exp(TOP_KILL_CONST[tree_type_] - TOP_KILL_H[tree_type_] * log(height_) + TOP_KILL_I[tree_type_] * sqrt(intensity)))
			/ (1. + exp(TOP_KILL_CONST[tree_type_] - TOP_KILL_H[tree_type_] * log(height_) + TOP_KILL_I[tree_type_] * sqrt(intensity)));

	return prob;
}

int clTree::getSeedProduction()
{
	int nSeeds = 0;

	if (reprod_strategy_ == 0 && age_ > 3650 && nCGT_ > 0)
	{
		nSeeds = (int) floor(nCGT_ / (SEED_WEIGHT[tree_type_]));
	}

	return nSeeds;
}

// Root suckering occurs when root biomass exceeds 20kg
int clTree::getRootSucker()
{
	int ret = 0;

	if (reprod_strategy_ == 1 && Br_ > 20.)
	{
		ret = 1;
		Br_ -= 1.;
	}

	return ret;
}

// calculate gc_ in molar units
void clTree::calGc(double A, double resp, double Ca, double P, double hs, double T)
{
	double cs;
	double gb = gMolar(T, gb_, P);    // convert gb_ from non-molar to molar

	if (gb > 0)
	{
		cs = Ca - CS_FACTOR * A * P / gb; 							//this line follows Arora
		gc_ = MyMax(0., M_TREE[tree_type_] * (A - resp) * hs / 100. * P / cs + B_TREE[tree_type_]);	//ie for negative An we produce zero gs
	}
	else
		gc_ = 0.;

	return;
}

void clTree::calGb(double wind)
{
	gb_ = GetgbCANOPY(wind, height_);
	return;
}

double clTree::getEt(double T, double P, double radnet_day, double s12, double SP_HEAT, double gama12, double rho12, double VPD12)
{
	if (Dormant_ == 1) return 0;

	double gc_nmolar = gNonMolar(T, gc_, P);  //gn
	return canopy_area_ * FOAPenMan(radnet_day, gb_, gc_nmolar, s12, SP_HEAT, gama12, rho12, VPD12);  //mm/day
}

double clTree::calQsum()
{
	return (1. - exp(-K_CAN_EXT_TREE[tree_type_] * (LAI_ + 0.4))) * K_CAN_EXT_TREE_INVERSE[tree_type_];
}

