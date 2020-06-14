#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>

using namespace std;

// for optimization, the target parameters are defined as global variables

int GLOB_YEAR;				// GLOBAL VARIABLE
double tmp_min_month = 0;	// GLOBAL VARIABLE!!!! Needed for Aindex
double phen_counter = 0;	// GLOBAL VARIABLE!!!! Counts swiches between metabolic and dormant state of tree 0
double A0_max_C3_day = 0;	// GLOBAL VARIABLE!!!! Index of day with maximum photosynthesis, needed for reproduction
double gs_C3_global = 0;
double gs_C4_global = 0;

double combustion_global = 0.;  // GLOBAL VARIABLE, biomass compustion by fire

#include "MyNrutil.h"
#include "MyMath.h"
#include "Globals.h"
#include "GlobalTypes.h"
#include "PenmanMonteith.h"
#include "Radiation.h"
#include "Leaf.h"
#include "InDataClass.h"
#include "InDataReaderClass.h"
#include "Fxn.h"
#include "Moisture.h"
#include "Fire.h"
#include "GrassPopClass.h"
#include "TreePopClass.h"
#include "YearlyDataWriterClass.h"
#include "SoilClass.h"

int main(int argc, char **argv)
{
//	abort if number of arguments is wrong
	if (argc < 9)
	{
		cerr << endl << "INI ERROR: wrong number of arguments." << endl;
		cerr << "INI ERROR: Needed arguments: ./Model     trees years lon lat fire rseed scenario output_dir" << endl;
		cerr << "INI ERROR: Abort simulation..." << endl << endl;
		return 1;
	}

	// initialize variables for simulation
	int count_years = 0;	// counter for years
	int month_in_year = 0;	// counter for month
	int day_in_year = 0;	// counter for days
	int fire_num = 0;
	int frost = 0;

	double sunhrs_day;					// sunhours per day
	double radnet_day;					// net radiation per day
	double mmsTOkgd;					// micromol to mol;mol to g; growth resp loss, carbon to dry mass; seconds in 12 hours
	double EtSite = 0;	// evapotranspiration, total
	double EtSiteGrass = 0;	// evapotranspiration, grass
	double EtSiteTrees = 0;	// evapotranspiration, trees
	double EtSiteGround = 0;	// evapotranspiration, ground

	double EtSiteRef;					// reference evapotranspiration
	double evapo_sum_year = 0;
	double rain_sum_year = 0;
	double rain_sum = 0;
	double evapo_sum = 0;
	double evapo_ref_sum = 0.;
	double evapo_grass_sum = 0;
	double evapo_soil_sum = 0;

	double total_fire_intensity = 0.;
	double total_npp = 0.;
	double total_nee = 0.;

	double T_fac = 0;  // needed for f(T) function of MaintResp

	arryYear Rain;					// array for rain distributon over the year
	arryYear Ignitions;				// stores days where ignition is possible
	for (int m = 0; m < 365; m++)
		Ignitions[m] = 0;

	clYearlyDataWriter myYearlyData(YEARLY_DATA_LENGTH);
	double yearly_output_data[YEARLY_DATA_LENGTH];

	// initialize variables with input arguments
	int num_of_trees = atoi(argv[1]);
	int years_to_run = atoi(argv[2]);
	int with_fire = atoi(argv[5]);
	double latitude = (double) atof(argv[4]); // y-coordinate
	double longitude = (double) atof(argv[3]); // x-coordinate
	double d_seed = (double) atof(argv[6]);
	int i_scen = (int) atoi(argv[7]);
	string s_seed = (string) argv[6];
	string s_scen = (string) argv[7];
	string output_directory = OUT_DATA_HOME + (string) argv[8] + "/";

	double d_params[INPUT_PARAMS];
	string s_params[INPUT_PARAMS];
	for (int i = 0; i < INPUT_PARAMS; i++)
		d_params[i] = 0;
	for (int i = 0; i < INPUT_PARAMS; i++)
		s_params[i] = "0";

	double DAILY_year;
	double DAILY_month;
	double DAILY_day;
	ifstream cs_co2_file;
	string cs_co2_file_name;
	string DAILY_file_name;

	if (i_scen == 61)
	{
		cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_RCP45.txt";     // RCP 4.5 with CO2, temp and rain
		DAILY_file_name = IN_DATA_HOME + "timeseries/sa_MPI_RCP45_b_" + (string) argv[3] + "_" + (string) argv[4] + ".dat";
	}

	if (i_scen == 71)
	{
		cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_RCP85.txt";     // RCP 8.5 with CO2, temp and rain
		DAILY_file_name = IN_DATA_HOME + "timeseries_85/sa_MPI_RCP85_b_" + (string) argv[3] + "_" + (string) argv[4] + ".dat";
	}

	ifstream DAILY_clim_file(DAILY_file_name.c_str());
	if (!DAILY_clim_file)
	{
		cout << endl << "INI ERROR: Cannot open " << DAILY_file_name << ". Exit simulation." << endl << endl;
		return 1;
	}

	cs_co2_file.open(cs_co2_file_name.c_str());
	if (!cs_co2_file)
	{
		cout << "INI ERROR: Cannot open CO2 file " << cs_co2_file_name << ", abort simulation." << endl << endl;
		return 1;
	}

	double ca_par_preassure = CA_PAR_PREASSURE;

	// read in databases/shortlist and initialize data
	clInDataReader MyReader;
	clInData IData = MyReader.getInData(longitude, latitude);

	if (IData.soil_N_[0] < 0. || IData.reh_[0] < 0. || IData.theta_wp_[0] < 0.)
	{
		cout << "INI ERROR: Cannot simulate the ocean, abort simulation: " << longitude << " " << latitude << endl << endl;
		return 1;
	}

	double soil_carbon = IData.soil_C_[0];
	clSoil mySoil(soil_carbon / 5000.); // initialize humus pool with soil C, from g/m^2 to kg/m^2

	double cs_preassure;

	cs_co2_file >> cs_preassure;
	ca_par_preassure = cs_preassure / 10.;

	IData.calcAtmospheric(latitude);
	IData.calcLeafPhoto(ca_par_preassure, soil_carbon);

	// calculate mean A0 and respiration, needed for output
	double A0C3_mean = 0;
	double A0C4_mean = 0;
	double Rl_mean = 0;
	double tmp_mean = 0;
	double sun_mean = 0;
	double tmp_min = 1000;
	double A0_max_C3 = 0;
	double hum_mean = 0;

	for (int ii = 0; ii < 12; ii++)
	{
		if (IData.A012C3_[ii] > A0_max_C3)
		{
			A0_max_C3 = IData.A012C3_[ii];
			A0_max_C3_day = ii;
		}
		A0C3_mean += IData.A012C3_[ii];
		A0C4_mean += IData.A012C4_[ii];
		Rl_mean += IData.RmL12C3_[ii];
		tmp_mean += IData.tmp_day_[ii];
		tmp_min = MyMin(IData.tmp_min_[ii], tmp_min);
		sun_mean += IData.sun_[ii];
		hum_mean += IData.reh_[ii];
// 		cout << IData.A012C3_[ii] << "  " << IData.A012C4_[ii] << endl;
	}

	A0_max_C3_day *= 31;
	A0C3_mean /= 12.;
	A0C4_mean /= 12.;
	Rl_mean /= 12.;
	tmp_mean /= 12.;
	sun_mean /= 12.;
	hum_mean /= 12.;
	gs_C3_global /= 12.;
	gs_C4_global /= 12.;

// 	cout << setw(14) << A0C3_mean << setw(14) << gs_C3_global << setw(14) << A0C4_mean << setw(14) << gs_C4_global << endl;

	GetNetRadiation(latitude, IData.sun_[month_in_year], IData.tmp_max_[month_in_year], IData.tmp_min_[month_in_year], IData.eA12_[month_in_year], day_in_year,
			&radnet_day, &sunhrs_day);

	mmsTOkgd = MMSTOKGD_HELPER * sunhrs_day;

	double diff_fc_wp[IData.soil_layers_];		// 1/(theta_fc-theta_wp), needed for GetSoilTheta

	for (size_t i = 0; i < IData.soil_layers_; i++)
	{
		if (IData.theta_fc_[i] > 0 && IData.theta_wp_[i] > 0) diff_fc_wp[i] = 1. / (IData.theta_fc_[i] - IData.theta_wp_[i]);
		else
			diff_fc_wp[i] = 0;
	}

	uint64_t nseed = (((uint64_t) (1122554162 + (int) d_seed)) << 16) + 0x330E;
	//nseed = (nseed <<16) +;
	gen.seed(nseed);

	clTreePop MyTreePop(IData.depth_[IData.soil_layers_ - 1]);
	clGrassPop MyGrassPop;

	{
		const int tree_type = MyRand() <= PROB_FOREST_TREE ? TR_FOR : TR_SAV;
		MyTreePop.addTree(100., tree_type);
	}

	if (num_of_trees > 1)
	{
		for (int count_trees = 2; count_trees <= num_of_trees; count_trees++)
		{
			const int tree_type = MyRand() <= PROB_FOREST_TREE ? TR_FOR : TR_SAV;

			MyTreePop.addTree(150. * MyRand(), tree_type);
		}
	}

	MyGrassPop.addGrass(INIT_MASS_GRASS, GR_C4_SAV);  // C4 savanna tree sub-canopy
	MyGrassPop.addGrass(INIT_MASS_GRASS, GR_C4_OPN);  // C4 between canopy
	MyGrassPop.addGrass(INIT_MASS_GRASS, GR_C4_FOR);  // C4 forest tree sub-canopy
	MyGrassPop.addGrass(INIT_MASS_GRASS, GR_C3_SAV);  // C3 savanna tree sub-canopy
	MyGrassPop.addGrass(INIT_MASS_GRASS, GR_C3_OPN);  // C3 between canopy
	MyGrassPop.addGrass(INIT_MASS_GRASS, GR_C3_FOR);  // C3 forest tree sub-canopy

// --------------------------------------------------------------------------------------------------------------
// --- START YEAR LOOP ------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------
	// loop over all years
	for (count_years = 1; count_years <= years_to_run; count_years++)
	{
		GLOB_YEAR = count_years;

#ifdef DECADE_EQUILLIBRIUM
		int YEAR_EQUI = 0;
		assert(DECADE_EQUILLIBRIUM%10==1);  //Parameter DECADE_EQUILLIBRIUM needs to be e.g. 1971, 1981, [...], 2091
		const int DECADE_EQUI = (DECADE_EQUILLIBRIUM-1961)/10;
		assert(DECADE_EQUI>=1);
		assert(DECADE_EQUI<=13);
        if( DECADE_EQUI<13 ) YEAR_EQUI   = 10*(DECADE_EQUI-1) + MyRound(MyRand()*10.);
        else                 YEAR_EQUI   = 10*(DECADE_EQUI-1) + MyRound(MyRand()* 9.);

        cs_co2_file.seekg(0, ios::beg);
        for (int jj=0; jj<YEAR_EQUI; jj++ ) cs_co2_file >> cs_preassure;
		for (int jj = 0; jj < (365 * (YEAR_EQUI - 1)); jj++)
			DAILY_clim_file >> DAILY_year >> DAILY_month >> DAILY_day >> IData.wnd_[month_in_year] >> Rain[day_in_year] >> IData.reh_[month_in_year]
					>> IData.tmp_min_[month_in_year] >> IData.tmp_max_[month_in_year];
#endif

		cs_co2_file >> cs_preassure;
		ca_par_preassure = cs_preassure / 10.;
		IData.calcLeafPhoto(ca_par_preassure, soil_carbon);

		// generate daily rain sequence and store it in Rain
		RainFallYear(DAYS_IN_YEAR, IData.rdo_, IData.pwet_, IData.ralpha_, IData.rbeta_, Rain);

		if (with_fire == 1 && count_years > 70) GetIgnitions(Ignitions, MyTreePop.getpCanopy());  // NOTE Ignitions start after 70 year spinup

		MyTreePop.setBornToZero();
		MyTreePop.setDeadToZero();

		evapo_sum_year = 0;
		rain_sum_year = 0;
		combustion_global = 0;

		// --------------------------------------------------------------------------------------------------------------
		// --- START DAY LOOP -------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------
		// loop over all days in a year
		for (day_in_year = 0; day_in_year < DAYS_IN_YEAR; day_in_year++)
		{

			month_in_year = (int) floor(((double) (day_in_year + 1)) / 30.42);

			T_fac = MaintRespTmpFac(IData.tmp_min_[month_in_year]);   // tmp function needed for respiration

			tmp_min_month = IData.tmp_min_[month_in_year];

			DAILY_clim_file >> DAILY_year >> DAILY_month >> DAILY_day >> IData.wnd_[month_in_year] >> Rain[day_in_year] >> IData.reh_[month_in_year]
					>> IData.tmp_min_[month_in_year] >> IData.tmp_max_[month_in_year];

			if (Rain[day_in_year] < 0) Rain[day_in_year] = 0;

			IData.tmp_[month_in_year] = IData.tmp_min_[month_in_year] + (IData.tmp_max_[month_in_year] - IData.tmp_min_[month_in_year]) / 2.;
			IData.tmp_day_[month_in_year] = IData.tmp_[month_in_year] + (IData.tmp_max_[month_in_year] - IData.tmp_[month_in_year]) / 2.;

			IData.calcAtmospheric(latitude);
			IData.calcLeafPhoto(ca_par_preassure, soil_carbon);
#ifndef DECADE_EQUILLIBRIUM
			//Spinup period: replay the first decade over the first 210 years
			if ((count_years % 10 == 0 && count_years <= 210 && day_in_year == 364) || (count_years == 210 && day_in_year == 364))
			{
				DAILY_clim_file.seekg(0, ios::beg);
			}
#endif

			GetNetRadiation(latitude, IData.sun_[month_in_year], IData.tmp_max_[month_in_year], IData.tmp_min_[month_in_year], IData.eA12_[month_in_year],
					day_in_year, &radnet_day, &sunhrs_day);

			mmsTOkgd = MMSTOKGD_HELPER * sunhrs_day;

			if (MyRand() < IData.frost_[month_in_year]) frost = 1;
			else
				frost = 0;

			MyTreePop.RunDeathProcess(frost);

			// compute reference evapotranspiration for soil decomposition (Yasso) ----------------------------------------
			EtSiteRef = FOAPenManRef(IData.tmp_[month_in_year], IData.tmp_[(month_in_year - 1) % 12], IData.s12_[month_in_year], radnet_day,
					IData.gama12_[month_in_year], WindAtHeight(MyTreePop.getMeanHeight(), IData.wnd_[month_in_year]), IData.eS12_[month_in_year],
					IData.eA12_[month_in_year]);

			// Evapotranspiration of soil ---------------------------------------------------------------------------------
			if (MyGrassPop.getAllDormant() <= 3)   // no respiration if 1/2 of grasses are active
			EtSiteGround = 0;
			else
			{
				double gb_soil = GetgbCANOPY(IData.wnd_[month_in_year], 2.);
				EtSiteGround = 0.0864 / 2.45 * exp(-4.28 + 11.97 * MyMin(0.35, IData.theta_[0])) * IData.VPD12_[month_in_year] / gb_soil;
			}

			// Evapotranspiration of grasses --------------------------------------------------------------------------------
			EtSiteGrass = MyGrassPop.getEt(IData.tmp_day_[month_in_year], IData.atm_press_, radnet_day, IData.s12_[month_in_year], SP_HEAT,
					IData.gama12_[month_in_year], IData.rho12_[month_in_year], IData.VPD12_[month_in_year], MyTreePop.getpCanopy());  //mm/day

			// Evapotranspiration of trees --------------------------------------------------------------------------------
			EtSiteTrees = MyTreePop.getEt(IData.tmp_day_[month_in_year], IData.atm_press_, radnet_day, IData.s12_[month_in_year], SP_HEAT,
					IData.gama12_[month_in_year], IData.rho12_[month_in_year], IData.VPD12_[month_in_year]);  //mm/day

			//-----------------------------------------------------------------------------------------
			// Evapotranspiration total
			EtSite = EtSiteGround + EtSiteTrees + EtSiteGrass;

			rain_sum_year += Rain[day_in_year];
			evapo_sum_year += EtSite;
			rain_sum += Rain[day_in_year];
			evapo_sum += EtSite;
			evapo_grass_sum += EtSiteGrass;
			evapo_soil_sum += EtSiteGround;
			evapo_ref_sum += EtSiteRef;

			if (Rain[day_in_year] > EtSite)
			{
				BucketIn(Rain[day_in_year] - EtSite, IData.theta_, IData.theta_fc_, IData.thickness_, IData.soil_layers_);
			}
			if (Rain[day_in_year] <= EtSite)
			{
				BucketOut(EtSite - Rain[day_in_year], IData.theta_, IData.theta_wp_, IData.thickness_, IData.soil_layers_);
			}

			std::vector<double> g_theta(IData.soil_layers_);
			GetSoilTheta(diff_fc_wp, IData.theta_, IData.theta_wp_, g_theta, IData.soil_layers_);

			// Run tree and grass physiology.
			MyTreePop.RunPhysiology(IData.A012C3_[month_in_year], IData.RmL12C3_[month_in_year], mmsTOkgd, IData.tmp_day_[month_in_year],
					IData.wnd_[month_in_year], g_theta.data(), ca_par_preassure, IData.atm_press_, IData.reh_[month_in_year], MyGrassPop.getHeight(GR_C4_OPN),
					MyGrassPop.getHeight(GR_C3_OPN), day_in_year, month_in_year, IData.theta_[0], IData.theta_fc_[1], IData.theta_wp_[1], T_fac, frost,
					MyGrassPop.getC34Ratio(), IData.thickness_, IData.soil_layers_);

			MyGrassPop.RunPhysiology(IData.A012C4_[month_in_year], IData.A012C3_[month_in_year], IData.RmL12C4_[month_in_year], IData.RmL12C3_[month_in_year],
					mmsTOkgd, IData.tmp_day_[month_in_year], IData.wnd_[month_in_year], g_theta.data(), ca_par_preassure, IData.atm_press_,
					IData.reh_[month_in_year], MyTreePop.getMeanSavHeight(), MyTreePop.getMeanForHeight(), T_fac, frost, MyTreePop.getpCanopySav(),
					MyTreePop.getpCanopyFor(), IData.thickness_, day_in_year);

			if (Ignitions[day_in_year] > 0)
			{
				double fire_intensity;

				double dead_fuel = MyGrassPop.getDryBiomassForFire();   // in kg/m^2

				double live_fuel = MyGrassPop.getWetBiomassForFire();

				double live_fuel_moisture = IData.reh_[month_in_year] / 100. + IData.pwet_[month_in_year];

				double dead_fuel_moisture = MyGrassPop.getFuelMoisture() * live_fuel_moisture;

				fire_intensity = LightFire(dead_fuel, live_fuel, dead_fuel_moisture, live_fuel_moisture, IData.wnd_[month_in_year], day_in_year);

				if (fire_intensity > 0)
				{
					double patchiness = Patchiness(fire_intensity);
					double scorch = ScorchHeight(fire_intensity);
					double cc_fine = CombComplFine(scorch);
					double cc_coarse = CombComplCoarse(scorch);
					double cc_heavy = CombComplHeavy(scorch);
					double cc_tk_helper = CombComplTopkillHelper(scorch);

					MyTreePop.setStateAfterFire(fire_intensity, patchiness, cc_fine, cc_coarse, cc_heavy, cc_tk_helper);
					MyGrassPop.setStateAfterFire(patchiness, cc_fine);

					fire_num++;
					total_fire_intensity += fire_intensity;

					total_nee -= (MyGrassPop.getLeafLiveCombustion() + MyGrassPop.getLeafDeadStCombustion() + MyGrassPop.getLeafDeadLyCombustion()
							+ MyTreePop.getStemLiveCombustion() + MyTreePop.getLeafLiveCombustion() + MyTreePop.getStemDeadStCombustion()
							+ MyTreePop.getStemDeadLyCombustion() + MyTreePop.getLeafDeadStCombustion() + MyTreePop.getLeafDeadLyCombustion()) * 0.44; // in kg C/m^2

				}  // if ( fire_intensity>0 )

			}  // if ( Ignitions[day_in_year]>0 )

			// Update soil pools, biomasses given in kg/m^2
			mySoil.UpdateCarbonPools(MyTreePop.getSoilNWL() + /* non woody litter    - tree fine roots and dead leaf*/
			MyGrassPop.getSoilNWL(), /* non woody litter    - grass roots and dead leaf*/
			MyTreePop.getSoilFWL(), /* fine woody litter   - tree coarse roots and fine stem */
			MyTreePop.getSoilCWL(), /* coarse woody litter - tree stem coarse and heavy */
			tmp_mean, (1. - MyTreePop.getDormant(0)) * (Rain[day_in_year] - EtSiteRef));

			total_npp += ((MyGrassPop.getGPP() - MyGrassPop.getRma() - MyGrassPop.getRgr() + MyTreePop.getGPP() - MyTreePop.getRma() - MyTreePop.getRgr())
					* 0.44); // in kg C/m^2
			total_nee += (((MyGrassPop.getGPP() - MyGrassPop.getRma() - MyGrassPop.getRgr() + MyTreePop.getGPP() - MyTreePop.getRma() - MyTreePop.getRgr())
					* 0.44) - mySoil.GetCarbonRelease()); // in kg C/m^2;

// ==========================================================================================================
//
//      IT FOLLOWS ONLY DATA OUTPUT AND THE ENDS OF THE DAY LOOP AND THE YEAR LOOP
//
// ==========================================================================================================

		} // End of day loop

// 		output at end of simulation
		{
			/**
			 * Definition of vector for annual data output, these data are written
			 * into YearlyData.***.dat in the output directory. If the length of
			 * this vector is changed, the constant YEARLY_DATA_LENGTH in the
			 * Globals.h must be changed.
			 */

			/*  1 */yearly_output_data[0] = longitude;
			/*  2 */yearly_output_data[1] = latitude;
			/*  3 */yearly_output_data[2] = count_years;
			/*  4 */yearly_output_data[3] = i_scen;
			/*  5 */yearly_output_data[4] = with_fire;
			/*  6 */yearly_output_data[5] = d_seed;
			/*  7 */yearly_output_data[6] = rain_sum_year; //rain_sum,
			/*  8 */yearly_output_data[7] = evapo_sum;
			/*  9 */yearly_output_data[8] = evapo_grass_sum;
			/* 10 */yearly_output_data[9] = evapo_soil_sum;
			/* 11 */yearly_output_data[10] = MyGrassPop.getBlMaxYearC4() * 10.; // t/ha
			/* 12 */yearly_output_data[11] = MyGrassPop.getBrMaxYearC4() * 10.; // t/ha
			/* 13 */yearly_output_data[12] = MyGrassPop.getBlMaxYearC3() * 10.; // t/ha
			/* 14 */yearly_output_data[13] = MyGrassPop.getBrMaxYearC3() * 10.; // t/ha
			/* 15 */yearly_output_data[14] = MyGrassPop.getLeafBmDeadSt() * 10. + MyGrassPop.getLeafBmDeadLy() * 10.; // t/ha
			/* 16 */
			yearly_output_data[15] = MyTreePop.getpCanopySavYearMean();
			/* 17 */yearly_output_data[16] = MyTreePop.getpCanopyForYearMean();
			/* 18 */yearly_output_data[17] = MyTreePop.getBlYearMax() * 10.; // t/ha
			/* 19 */yearly_output_data[18] = MyTreePop.getBsYearMax() * 10.; // t/ha
			/* 20 */yearly_output_data[19] = MyTreePop.getBrYearMax() * 10.; // t/ha
			/* 21 */yearly_output_data[20] = MyGrassPop.getC34RatioYearMean();
			/* 22 */yearly_output_data[21] = MyTreePop.getMeanHeightYearMean();
			/* 23 */yearly_output_data[22] = MyTreePop.getMaxHeightYearMean();
			/* 24 */yearly_output_data[23] = MyTreePop.getPopSize();
			/* 25 */yearly_output_data[24] = MyTreePop.getSavTreeNum();
			/* 26 */yearly_output_data[25] = MyTreePop.getForTreeNum();
			/* 27 */yearly_output_data[26] = MyTreePop.getMaxBasalAreaYearMean();
			/* 28 */yearly_output_data[27] = total_npp / (double) count_years; // kg C/m^2
			/* 29 */yearly_output_data[28] = total_nee / (double) count_years;
			/* 30 */yearly_output_data[29] = mySoil.GetCarbonStored();
			/* 31 */yearly_output_data[30] = phen_counter;
			/* 32 */yearly_output_data[31] = fire_num;
			/* 33 */yearly_output_data[32] = total_fire_intensity / (double) fire_num;
			/* 34 */yearly_output_data[33] = IData.soil_N_[0];
			/* 35 */yearly_output_data[34] = IData.soil_C_[0];
			/* 36 */yearly_output_data[35] = tmp_mean;
			/* 37 */yearly_output_data[36] = ca_par_preassure * 10.;
			/* 38 */yearly_output_data[37] = d_params[0];
			/* 39 */yearly_output_data[38] = d_params[1];
			/* 40 */yearly_output_data[39] = d_params[2];
			/* 41 */yearly_output_data[40] = d_params[3];
			/* 42 */yearly_output_data[41] = d_params[4];
			/* 43 */yearly_output_data[42] = d_params[5];
			/* 44 */yearly_output_data[43] = d_params[6];
			/* 45 */yearly_output_data[44] = d_params[7];
			/* 46 */yearly_output_data[45] = d_params[8];
			/* 47 */yearly_output_data[46] = combustion_global; //d_params[9];
			/* 48 */yearly_output_data[47] = MyTreePop.getActiveDays();
			/* 49 */yearly_output_data[48] = MyGrassPop.getActiveDays();
			/* 50 */yearly_output_data[49] = evapo_ref_sum;
			/* 51 */yearly_output_data[50] = A0C3_mean;
			/* 52 */yearly_output_data[51] = A0C4_mean;
			/* 53 */yearly_output_data[52] = gs_C3_global;
			/* 54 */yearly_output_data[53] = gs_C4_global;

			myYearlyData.setYearlyData(yearly_output_data);
			myYearlyData.printYearlyDataToFile(output_directory);
		}

	} // end of year loop

	return 0;
}

