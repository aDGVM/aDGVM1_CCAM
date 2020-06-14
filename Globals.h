#pragma once

// Model directories

// Home directory of model, contains the source code and the executable
const string MODEL_HOME = "./";

// Home directory of input data, contains
//    filenames.txt        Contain names of global databases,
//                         these databases must be in the same directory
//    shortlists.txt       Contains the names of shortlists that have already been read
//                         for study sites
//    shortlists           Directory that contains shortlists that have already been read
//    ClimateChange        Directory that contains files with climate change scenarios
//                         (CO2, precipitation, temperature)
//    fire_sites           Directory with fire scenarios for specific study sites, these
//                         scenarios are used when compiler flag S_FIRE_FROM_FILE is set
//    prec_sites           Directory with precipitation data for specific study sites, these
//                         data are used when compiler flag S_RAIN_FROM_FILE is set
const string IN_DATA_HOME = "./data/";

// Home directory of output data, this directory contains sub-directories (experiment_xyz).
// Output data are written in OUT_DATA_HOME/(experiment_xyz),
// the directory (experiment_xyz) is specified in the model command line.
const string OUT_DATA_HOME = "./OutputData/";

// Model constants

// Start constants for the leaf photosynthesis model

const double ABS_PHOTONS_C3 = 0.86;			// absorbtance to incident flux of photons, Collatz C3
const double ALPHA_C3 = 0.08;			// intrinsic quantum efficiency for C02 uptake, Collatz C3
const double ALPHAR_F_C4 = 0.067;		// product of alpha_r (intrinsic quantum yield of C3
const double ABS_PHOTONS_C4 = 0.80;			// absorbtance to incident flux of photons, Collatz C4
// photosynthesis) and f (fraction of absorbed photons used
// by C3 reactions), FJiC4, Collatz C4 units mol/mol
const double KAPPA_C4 = 0.7 * 1e6;		// initial slope of photosynthetic CO2 resp for C4, 0.7 mumol/m2/s
// FC4Jc, Collatz C4
const double CA_PAR_PREASSURE = 38.70;		// Atmospheric partial pressure C02

const double OI_PAR_PREASSURE = 21.0;			// O2 Partial pressure KPa, Collatz

const double MTR = 9.;

const double B_TREE[2] = { 0.01, 0.01 };   // ball berry constant conductance units (take care with reporting units here)
const double M_TREE[2] = { MTR, MTR };     // ball berry constant dimensionless

const double B_GRASS[6] = { 0.04, 0.04, 0.04, 0.04, 0.04, 0.04 }; // ball berry constant conductance units
//  const double M_GRASS[6]                 = { 3.75, 3.75, 3.75, 4.0, 4.0, 4.0 }; // ball berry constant const double M_GRASS[6]                 = { 4.0, 4.0, 4.0, 4.0, 4.0, 4.0 }; // ball berry constant dimensionless
const double M_GRASS[6] = { 4., 4., 4., MTR, MTR, MTR }; // ball berry constant dimensionless

const double KC_K25_C3 = 30.0; 	   	// c3 value of Kc at 25degC, Pa, Collatz C3
const double KC_Q10_C3 = 2.1;   		// c3 rate of change of Kc with temp, Collatz C3
const double KO_K25_C3 = 30.0;    		// c3 value of Ko at 25degC, KPa, Collatz C3
const double KO_Q10_C3 = 1.2;   		// c3 rate of change of Ko with temp, Collatz C3
const double TAU_K25_C3 = 2600.0;  		// c3 value of tau at 25degC,  Collatz C3
const double TAU_Q10_C3 = 0.57;  		// c3 rate of change of tau with temp, Collatz C3

const double KC_K25_C4 = 140.0;	   	// c4 value of Kc at 25degC, Pa, Collatz C4
const double KC_Q10_C4 = 2.1;   		// c4 rate of change of Kc with temp, Chen
const double KO_K25_C4 = 34.0; 	   	// c4 value of Ko at 25degC, KPa, Collatz C4
const double KO_Q10_C4 = 1.2;   		// c4 rate of change of Ko with temp, Chen
const double TAU_K25_C4 = 2600.0;  		// c4 value of tau at 25degC,  Collatz C4
const double TAU_Q10_C4 = 0.67;  		// c4 rate of change of tau with temp, Chen

// const double V_MAX_CONST_C3				= 0.8;
// const double V_MAX_CONST_C4				= 0.4;

// end constants for the leaf photosynthesis model

// start constants needed for leaf boundary layer conductance

const double CLD_TREE = 0.02;			// charac leaf dimension tree (assumes 2 cm leaf width)
const double CLD_GRASS = 0.005;  		// charac leaf dimesnion grass (assumes 5 mm leaf width)

// end constants needed for leaf boundary layer conductance

// start constants needed for (canopy) boundary layer conductance
const double Z_WIND = 10.;			// height at which reference windspeed is measured, New et al. 2000
const double ZD_CONST = 0.86;			// constant for estimating displacement height, Jones 1992
const double Z0_CONST = 0.06;			// constant for estimating roughness length, Jones 1992
const double VEGETATION_HEIGHT = 1.5;			// H_bar average vegetation height in m (needed to calculate wind)
const double ROUGHNESS_LENGTH = Z0_CONST * VEGETATION_HEIGHT;
const double DISPLACEMENT_HEIGHT = ZD_CONST * VEGETATION_HEIGHT;
const double REF_HEIGHT_Z = 10.;			// Height, at which wind is meassured
const double WIND_AT_HEIGHT_HELPER = 1. / log((REF_HEIGHT_Z - DISPLACEMENT_HEIGHT) / ROUGHNESS_LENGTH);

const double KARMAN_CONST = 0.41;			// von Karman constant, for canopy boundary layer conductance
const double KARMAN_CONST_QUAD = pow(KARMAN_CONST, 2.);	// quadrat of von Karman constant

// Start CANOPY constants
const double K_CAN_EXT_TREE[2] = { 0.5, 0.4 };
const double K_CAN_EXT_TREE_INVERSE[2] = { 1. / K_CAN_EXT_TREE[0], 1. / K_CAN_EXT_TREE[1] };

const double K_CAN_EXT_GRASS[6] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 }; // canopy extinction coefficient, Woodward; Arora et al. 2003 uses
const double K_CAN_EXT_GRASS_INVERSE[6] = { 1. / K_CAN_EXT_GRASS[0], 1. / K_CAN_EXT_GRASS[1], 1. / K_CAN_EXT_GRASS[2], 1. / K_CAN_EXT_GRASS[3], 1.
		/ K_CAN_EXT_GRASS[4], 1. / K_CAN_EXT_GRASS[5] };

// Start respiration constants
const double R_MAINT_RESP_TR[2] = { 0.015, 0.015 }; // tree RmL leaf  maint respiration as proportion of Vm, value from Collatz  original
const double R_MAINT_RESP_GR[6] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025 }; // grass RmL leaf  maint respiration as proportion of Vm, value from

// Start allometric constants, Table 5
const double HEIGHT_C1_TREE[2] = { 1.295720, 1.295720 };    // parameters for height(StemBiomass), Higgins2007
const double HEIGHT_C2_TREE[2] = { 0.392157, 0.392157 };    // parameters for height(StemBiomass), Higgins2007

const double HEIGHT_C1_GRASS[6] = { 3.5, 3.5, 3.5, 3.5, 3.5, 3.5 };  // parameters for height(StemBiomass), Arora2005
const double HEIGHT_C2_GRASS[6] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };  // parameters for height(StemBiomass), Arora2005

const double STEM_AREA_C1[2] = { 2.797, 2.797 };        // parameters for basal area, Higgins (2007) Ecology
const double STEM_AREA_C2[2] = { 200., 200. };        // parameters for basal area, Higgins (2007) Ecology
const double STEM_AREA_HELPER[2] = { STEM_AREA_C1[0] / STEM_AREA_C2[0], STEM_AREA_C1[1] / STEM_AREA_C2[1] };

const double GAMMA_CANOPY_TREE[2] = { 0.37, 0.42 };        // ratio of canopy radius to height for trees
const double GAMMA_CANOPY_GRASS[6] = { 0.4, 0.4, 0.4, 0.4, 0.4, 0.4 };  // ratio of canopy radius to height for grasses (estimated)

const double SLA_TREE[2] = { 10., 12. };         // specific leaf area tree, m2/kg from ???
const double SLA_GRASS[6] = { 10.9, 10.9, 10.9, 10.9, 10.9, 10.9 }; // specific leaf area grass, m2/kg from Scholes&Walker

const double MAX_ROOT_DEP_TREE = 2.;			// maximum rooting depth (m) trees -see Schenk and Jackson papers for estimates
const double ROOT_DENSITY_TREE = 0.1e03;		// density of root biomass (kg/m3) assumes root biomass is a light wood
const double MIN_ROOT_RAD_TREE = 0.015;		// minimum root radius, m

const double MAX_ROOT_DEP_GRASS = 0.3;			// maximum rooting depth (m) grasses
const double ROOT_DENSITY_GRASS = 0.1e03;		// density of root biomass (kg/m3) assumes root biomass is a light wood
const double MIN_ROOT_RAD_GRASS = 0.005;		// minimum root radius
// end allometric constants

// Start constants describing respiration, Table 7
const double BETA_N = 0.218;  		// respiration rate for roots and stems kg C per Kg N, Arora et al. 2003
// (after Keyser et al. 2000)

const double BETA_STEM_TREE = 0.025;		// fraction of respirating tissue
const double BETA_ROOT_TREE = 0.025;		// fraction of respirating tissue
const double BETA_ROOT_GRASS = 0.15;			// fraction of respirating tissue

const double UPSILON_STEM_TREE = 150.;			// C:N for tree stems note this needs to be higher than root value
const double UPSILON_ROOT_TREE = 60.;			// C:N for tree roots
const double UPSILON_STEM_GRASS[6] = { 120., 120., 120., 120., 120., 120. }; // C:N ratio, C4 have higher C:N than C3
const double UPSILON_ROOT_GRASS[6] = { 120., 120., 120., 120., 120., 120. }; // C:N ratio

const double SIGMA_GROW_RESP_TREE[2] = { 0.35, 0.35 };			// growth respiration constant (Arora =0.35)			// NEW
// const double SIGMA_GROW_RESP_GRASS[6]   = { 0.15, 0.15, 0.15, 0.35, 0.35, 0.35 }; // growth respiration constant (Arora =0.35)            // NEW
const double SIGMA_GROW_RESP_GRASS[6] = { 0.35, 0.35, 0.35, 0.35, 0.35, 0.35 }; // growth respiration constant (Arora =0.35)  // NOTE changed that for za!!
// End constants describing respiration

// Start allocation constants, Table 8
const double A0_ROOT_GRASS[6] = { 0.4, 0.4, 0.4, 0.4, 0.4, 0.4 };        //optimal resource allocation proportion to root
const double A0_STEM_GRASS[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };         //optimal resource allocation proportion to stem
const double A0_LEAF_GRASS[6] = { 0.6, 0.6, 0.6, 0.6, 0.6, 0.6 };        //optimal resource allocation proportion to leaf

const double A0_ROOT_TREE[2] = { 0.35, 0.15 };        // optimal resource allocation proportion to root
const double A0_STEM_TREE[2] = { 0.35, 0.50 };        // optimal resource allocation proportion to stem
const double A0_LEAF_TREE[2] = { 0.30, 0.35 };        // optimal resource allocation proportion to leaf

const double A0_LEAF_TREE_INV[2] = { 1. / A0_LEAF_TREE[0], 1. / A0_LEAF_TREE[1] };
const double A0_ROOT_TREE_HELPER[2] = { 1. + A0_ROOT_TREE[0], 1. + A0_ROOT_TREE[1] };
const double A0_STEM_TREE_HELPER[2] = { 1. + A0_STEM_TREE[0], 1. + A0_STEM_TREE[1] };

const double ALLOC_DENOM_HELPER[2] = { 3. + A0_ROOT_TREE[0] + A0_STEM_TREE[0], 3. + A0_ROOT_TREE[1] + A0_STEM_TREE[1] };
// End allocation constants

// Start phenology constants, Table 9
const double STRESS_INDEX_TREE[2] = { 0., 0. };    // 3. stress index for trees
const double STRESS_INDEX_GRASS[6] = { 0., 0., 0., 0., 0., 0. };        // 5. stress index for grass

const double REM_BM_GRASS = 0.1;			// proportion of grass leaf biomass kept after leaf abscission
const double REM_BM_TREE = 0.01;			// proportion of tree leaf biomass kept after leaf abscission

const int D_POS_TREE[2] = { 10, 10 };    // 5 days of positive carbon gain while dormant before growth restarts (days)
const int D_NEG_TREE[2] = { 7, 7 };    // 7 days of negative carbon gain before leaf fall (days)

const int D_POS_GRASS[6] = { 7, 7, 7, 7, 7, 7 };  // days of positive carbon gain while dormant before growth restarts (days)
const int D_NEG_GRASS[6] = { 5, 5, 5, 5, 5, 5 };  // days of negative carbon gain before leaf fall (days)

const double TI_CONST = 15.;
// End phenology constants

// Start turnover parameters, Table 10
// trees
const double EXTINC_EXP_TREE = -0.909;    // exponent for leaf longevity
const double EXTINC_FAC_TREE = 29664.;    // factor for leaf longevity
const double OMEGA_LEAF_TREE[2] = { EXTINC_FAC_TREE * pow(10. * SLA_TREE[0], EXTINC_EXP_TREE), EXTINC_FAC_TREE * pow(10. * SLA_TREE[1], EXTINC_EXP_TREE) }; // factor 10 to transform units
const double OMEGA_ROOT_TREE = 1e6;       // root longevity in days
const double OMEGA_STEM_TREE = 1e6;       // stem longevity in days
// const double OMEGA_ROOT_TREE			= 9125.;     // (root longevity in days (Gill and Jackson 2000)
// const double OMEGA_STEM_TREE			= 9125.;     // stem longevity in days (no source)
const double OMEGA_LEAF_TR_TREE = 30.;       // transition from hanging to lying dead leaf
const double OMEGA_STEM_TR_TREE = 300.;      // transition from standing to lying dead stem biomass

// grass
const double EXTINC_EXP_GRASS = -0.909;    // exponent for leaf longevity
const double EXTINC_FAC_GRASS = 29664.;    // factor for leaf longevity
const double OMEGA_LEAF_GRASS[6] = { EXTINC_FAC_GRASS * pow(10. * SLA_GRASS[0], EXTINC_EXP_GRASS), EXTINC_FAC_GRASS * pow(10. * SLA_GRASS[1], EXTINC_EXP_GRASS),
		EXTINC_FAC_GRASS * pow(10. * SLA_GRASS[2], EXTINC_EXP_GRASS), EXTINC_FAC_GRASS * pow(10. * SLA_GRASS[3], EXTINC_EXP_GRASS), EXTINC_FAC_GRASS
				* pow(10. * SLA_GRASS[4], EXTINC_EXP_GRASS), EXTINC_FAC_GRASS * pow(10. * SLA_GRASS[5], EXTINC_EXP_GRASS) };

const double OMEGA_ROOT_GRASS = 384.;      // (root longevity in days (Gill and Jackson 2000)
const double OMEGA_STEM_GRASS = 384.;      // stem longevity in days (no source)
const double OMEGA_LEAF_TR_GRASS = 150.;     // transition rate from standing to lying dead leaf biomass

// helper used for calculations

const double HELPER_BLD_TREE[2] = { 1. / OMEGA_LEAF_TREE[0], 1. / OMEGA_LEAF_TREE[1] };		// for turnover
const double HELPER_BL__TREE[2] = { 1. - 1. / OMEGA_LEAF_TREE[0], 1. - 1. / OMEGA_LEAF_TREE[1] };		// for turnover

const double HELPER_BSD_TREE = 1. / OMEGA_STEM_TREE;		// for turnover
const double HELPER_BS__TREE = 1. - 1. / OMEGA_STEM_TREE;		// for turnover

const double HELPER_BRD_TREE = 1. / OMEGA_ROOT_TREE;		// for turnover
const double HELPER_BR__TREE = 1. - 1. / OMEGA_ROOT_TREE;		// for turnover

const double HELPER_BLTL_TREE = 1. / OMEGA_LEAF_TR_TREE;		// for turnover
const double HELPER_BLT__TREE = 1. - 1. / OMEGA_LEAF_TR_TREE;		// for turnover

const double HELPER_BSTL_TREE = 1. / OMEGA_STEM_TR_TREE;		// for turnover
const double HELPER_BST__TREE = 1. - 1. / OMEGA_STEM_TR_TREE;		// for turnover

const double HELPER_BLS_GRASS[6] = { 1. / OMEGA_LEAF_GRASS[0], 1. / OMEGA_LEAF_GRASS[1], 1. / OMEGA_LEAF_GRASS[2], 1. / OMEGA_LEAF_GRASS[3], 1.
		/ OMEGA_LEAF_GRASS[4], 1. / OMEGA_LEAF_GRASS[5] };	// for turnover
const double HELPER_BL__GRASS[6] = { 1. - 1. / OMEGA_LEAF_GRASS[0], 1. - 1. / OMEGA_LEAF_GRASS[1], 1. - 1. / OMEGA_LEAF_GRASS[2], 1. - 1. / OMEGA_LEAF_GRASS[3],
		1. - 1. / OMEGA_LEAF_GRASS[4], 1. - 1. / OMEGA_LEAF_GRASS[5] };	// for turnover

const double HELPER_BSS_GRASS = 1. / OMEGA_STEM_GRASS;		// for turnover
const double HELPER_BS__GRASS = 1. - 1. / OMEGA_STEM_GRASS;		// for turnover

const double HELPER_BRD_GRASS = 1. / OMEGA_ROOT_GRASS;		// for turnover
const double HELPER_BR__GRASS = 1. - 1. / OMEGA_ROOT_GRASS;		// for turnover

const double HELPER_BLTL_GRASS = 1. / OMEGA_LEAF_TR_GRASS;	// for turnover
const double HELPER_BLT__GRASS = 1. - 1. / OMEGA_LEAF_TR_GRASS;	// for turnover
// End turnover parameters

// Start constant for transitions from dead (lying) biomass pools into Yasso soil model pools
const double OMEGA_LD_NWL_TREE = 200.;		// dead leaves to non woody litter
const double OMEGA_RF_NWL_TREE = 5.;		// fine roots to non woody litter
const double OMEGA_RC_FWL_TREE = 5.;		// coarse roots to fine woody litter
const double OMEGA_SF_FWL_TREE = 200.;		// fine stem biomass to fine woody litter
const double OMEGA_SC_CWL_TREE = 500.;		// coarse/heavy stem biomas to coarse woody litter

const double OMEGA_LD_NWL_GRASS = 200.;		// dead leaves to non woody litter
const double OMEGA_RF_NWL_GRASS = 5.;		// fine roots to non woody litter
// 
const double HELPER_LD_NWL_TREE = 1. / OMEGA_LD_NWL_TREE;
const double HELPER_RF_NWL_TREE = 1. / OMEGA_RF_NWL_TREE;
const double HELPER_RC_FWL_TREE = 1. / OMEGA_RC_FWL_TREE;
const double HELPER_SF_FWL_TREE = 1. / OMEGA_SF_FWL_TREE;
const double HELPER_SC_CWL_TREE = 1. / OMEGA_SC_CWL_TREE;

const double HELPER_LD_NWL_GRASS = 1. / OMEGA_LD_NWL_GRASS;
const double HELPER_RF_NWL_GRASS = 1. / OMEGA_RF_NWL_GRASS;
// End constant for transitions from dead biomass pools to Yasso soil model pools

// Start parameters for water index in Aindex
const double WATER_IND_TREE[2] = { 0.66, 0.66 };
const double WATER_IND_TREE_1[2] = { 1. - WATER_IND_TREE[0], 1. - WATER_IND_TREE[1] };
const double WATER_IND_GRASS[6] = { 0.66, 0.66, 0.66, 0.66, 0.66, 0.66 };
const double WATER_IND_GRASS_1[6] = { 1. - WATER_IND_GRASS[0], 1. - WATER_IND_GRASS[1], 1. - WATER_IND_GRASS[2], 1. - WATER_IND_GRASS[3], 1.
		- WATER_IND_GRASS[4], 1. - WATER_IND_GRASS[5] };
// End parameters for water index in Aindex

// Start constants for fire model, Table 14
// constants for fire intensity (Higgins 2008)
const double FIRE_H = 16890.;
const double FIRE_C = 301.;
const double FIRE_AW = 119.7;
const double FIRE_QM = 2.6e6;
const double FIRE_QV = 160749.;

// constants for topkill probability (Higgins et al. 2000)
// +- arbitrary values for forest trees that ensure that P(topkill)=1
const double TOP_KILL_CONST[2] = { 4.3, 5.3 };
const double TOP_KILL_H[2] = { 5.003, 4.003 };
const double TOP_KILL_I[2] = { 0.004408, 0.005408 };

// constants for fire ignition
// const double IGNITION_MIN_INT			= 300.;				// minimum fire intensity for ignition
const double IGNITION_MIN_INT = 150.;				// minimum fire intensity for ignition
const double IGNITION_PROB = 0.01;				// probability of fire ignition
const double IGNITION_PAR_2 = 0.1;				// parameter describing fire ignition sequence, treecover
const double IGNITION_PAR_1 = -0.083333;		// parameter describing fire ignition sequence, treecover

const double SEED_GERM_PROB[2] = { 0.25, 0.25 };
const double PROB_ROOT_SUCKER[2] = { 0., 0. };
// const double PROB_ROOT_SUCKER[2]		= { 0.25, 0.05};

const double RESPROUTING_PROB[2] = { 0.995, 0.8 };

#if defined S_UNFLAMABLE_C4
const double DESIC_COEFF[6]				= { 0.995, 0.995, 0.995, 0.995, 0.995, 0.995 };  // C4 dries out faster (Bond 2008)
#else
const double DESIC_COEFF[6] = { 0.95, 0.95, 0.95, 0.995, 0.995, 0.995 };  // C4 dries out faster (Bond 2008)
#endif
// End constants for fire model

// Start constants for light competition model, Table
const double LIGHT_COMP_TREE_TREE_1 = 0.35;		// tree shades tree
const double LIGHT_COMP_TREE_TREE_2 = 1. - LIGHT_COMP_TREE_TREE_1;
const double LIGHT_COMP_GRASS_TREE_1 = 0.25;		// grass shades tree
const double LIGHT_COMP_GRASS_TREE_2 = 1. - LIGHT_COMP_GRASS_TREE_1;
const double LIGHT_COMP_GRASS_GRASS_1 = 0.5;				// grass shades grass
const double LIGHT_COMP_GRASS_GRASS_2 = 1. - LIGHT_COMP_GRASS_GRASS_1;
const double LIGHT_COMP_TREE_GRASS_1 = 0.4;		// tree shades grass
const double LIGHT_COMP_TREE_GRASS_2 = 1. - LIGHT_COMP_TREE_GRASS_1;
// end constants for light competition model

// Start constants for light competition model, Table

// parameter L_XXX_YYY describes how competitor XXX shades target YYY
const double L_SAV_SAV = {0.5};
const double L_SAV_FOR = {0.15};
const double L_SAV_C4G = {0.5};
const double L_SAV_C3G = {0.15};

const double L_FOR_SAV = {0.5};
const double L_FOR_FOR = {0.15};
const double L_FOR_C4G = {0.5};
const double L_FOR_C3G = {0.15};

const double L_C4G_SAV = {0.5};
const double L_C4G_FOR = {0.15};
const double L_C4G_C4G = {0.5};
const double L_C4G_C3G = {0.15};

const double L_C3G_SAV = {0.5};
const double L_C3G_FOR = {0.15};
const double L_C3G_C4G = {0.5};
const double L_C3G_C3G = {0.15};

const double LICMP_1[8][8] =
		{ L_SAV_SAV, L_SAV_FOR, L_SAV_C4G, L_SAV_C4G, L_SAV_C4G, L_SAV_C3G, L_SAV_C3G, L_SAV_C3G, L_FOR_SAV, L_FOR_FOR, L_FOR_C4G, L_FOR_C4G, L_FOR_C4G,
				L_FOR_C3G, L_FOR_C3G, L_FOR_C3G, L_C4G_SAV, L_C4G_FOR, L_C4G_C4G, L_C4G_C4G, L_C4G_C4G, L_C4G_C3G, L_C4G_C3G, L_C4G_C3G, L_C4G_SAV, L_C4G_FOR,
				L_C4G_C4G, L_C4G_C4G, L_C4G_C4G, L_C4G_C3G, L_C4G_C3G, L_C4G_C3G, L_C4G_SAV, L_C4G_FOR, L_C4G_C4G, L_C4G_C4G, L_C4G_C4G, L_C4G_C3G, L_C4G_C3G,
				L_C4G_C3G, L_C3G_SAV, L_C3G_FOR, L_C3G_C4G, L_C3G_C4G, L_C3G_C4G, L_C3G_C3G, L_C3G_C3G, L_C3G_C3G, L_C3G_SAV, L_C3G_FOR, L_C3G_C4G, L_C3G_C4G,
				L_C3G_C4G, L_C3G_C3G, L_C3G_C3G, L_C3G_C3G, L_C3G_SAV, L_C3G_FOR, L_C3G_C4G, L_C3G_C4G, L_C3G_C4G, L_C3G_C3G, L_C3G_C3G, L_C3G_C3G };

const double LICMP_2[8][8] = { 1. - L_SAV_SAV, 1. - L_SAV_FOR, 1. - L_SAV_C4G, 1. - L_SAV_C4G, 1. - L_SAV_C4G, 1. - L_SAV_C3G, 1. - L_SAV_C3G, 1. - L_SAV_C3G,
		1. - L_FOR_SAV, 1. - L_FOR_FOR, 1. - L_FOR_C4G, 1. - L_FOR_C4G, 1. - L_FOR_C4G, 1. - L_FOR_C3G, 1. - L_FOR_C3G, 1. - L_FOR_C3G, 1. - L_C4G_SAV, 1.
				- L_C4G_FOR, 1. - L_C4G_C4G, 1. - L_C4G_C4G, 1. - L_C4G_C4G, 1. - L_C4G_C3G, 1. - L_C4G_C3G, 1. - L_C4G_C3G, 1. - L_C4G_SAV, 1. - L_C4G_FOR, 1.
				- L_C4G_C4G, 1. - L_C4G_C4G, 1. - L_C4G_C4G, 1. - L_C4G_C3G, 1. - L_C4G_C3G, 1. - L_C4G_C3G, 1. - L_C4G_SAV, 1. - L_C4G_FOR, 1. - L_C4G_C4G, 1.
				- L_C4G_C4G, 1. - L_C4G_C4G, 1. - L_C4G_C3G, 1. - L_C4G_C3G, 1. - L_C4G_C3G, 1. - L_C3G_SAV, 1. - L_C3G_FOR, 1. - L_C3G_C4G, 1. - L_C3G_C4G, 1.
				- L_C3G_C4G, 1. - L_C3G_C3G, 1. - L_C3G_C3G, 1. - L_C3G_C3G, 1. - L_C3G_SAV, 1. - L_C3G_FOR, 1. - L_C3G_C4G, 1. - L_C3G_C4G, 1. - L_C3G_C4G, 1.
				- L_C3G_C3G, 1. - L_C3G_C3G, 1. - L_C3G_C3G, 1. - L_C3G_SAV, 1. - L_C3G_FOR, 1. - L_C3G_C4G, 1. - L_C3G_C4G, 1. - L_C3G_C4G, 1. - L_C3G_C3G, 1.
				- L_C3G_C3G, 1. - L_C3G_C3G };
// end constants for light competition model

// Start constants for reproduction and mortality, Table 13
const double SEED_WEIGHT[2] = { 0.001, 0.001 };    // weight of a single seed in kg (Hovstadt 1999)
const double SEED_FRAC_DAY[2] = { 0.1, 0.1 };
const double SEED_DECAY_RATE[2] = { 0.9967, 0.9967 };       // Higgins 2000

const double DEATH_PROB_FROST[2] = { 0.001, 0.001 };
const double DEATH_PROB_CARBON[2] = { 0.001, 0.001 };
const double DEATH_PROB_COMP[2] = { 0.001, 0.0005 };
// End constants for reproduction and mortality, Table 13

// Start constants for soil, Table 15
const size_t    N_SOIL_LAYERS        = 12;    		// Number of soil layers
// const double DEPTH[SOIL_LAYERS] = {      0.05,     0.1,    0.2,    0.3,    0.4,    0.6,    0.8,    1.0,     1.25,     1.5,     1.75,     2.0 };
// const double THICK[SOIL_LAYERS] = { 0.05,     0.05,    0.1,    0.1,    0.1,    0.2,    0.2,    0.2,    0.25,     0.25,    0.25,     0.25,    };
// end constants for soil

// parameters for competitor probability
const double COMP_PAR_1 = 1.1;
const double COMP_PAR_2 = 1.5;
const double COMP_PAR_3 = 1e-4;
// end parameters for competitor probability

// constants to avoid zero in denominator
const double MIN_LAI = 0.000001;	// if the LAI of a plant is 0, we have a problem with our formulas
const double MIN_HEIGHT = 0.000001;

// other constants
const int DAYS_IN_YEAR = 365;			// Days in a year
const int MONTH_IN_YEAR = 12;
const double GRID_SIZE = 10000.;		// 1ha

//Penman Monteith constants, source: mostly FAO
const double SGC = 0.287;      	// specific gas constant, KJ/kg/K, FOArho}
const double SP_HEAT = 1.013e-3;   	// MJ/kg/degC = J/g/degC from FOA specific heat of moist air}
const double LAMBDA = 2.45;       	// MJ/kg from FOA latent heat of air at 20degC}
const double GSC = 0.082;       	// solar constant - extrateresstrial solar radiation MJ/m2/minute
// const double	GSC						= 0.082*24.;       	// solar constant - extrateresstrial solar radiation MJ/m2/day
const double ANGSTRONG_A = 0.25;        	// constant defining prop of extrat. radiation reaching
// earth on overcast days
const double ANGSTRONG_B = 0.50;        	// constant defining the additional prop of extrat. radiation
// reaching earth on clear days
const double ALBEDO = 0.23;        	// canopy reflection coefficient estimate for a grass
// reference crop (dimensionless)
const double SBC = 4.903e-09;   	// stephan bolzman constant MJ.k^-4.day^-1

// globals from fxn.h
// const double	E25H_CONST_1			= 0.807;
// const double	E25H_CONST_2			= 2.137;
const double CS_FACTOR = 1.4;
const double MAIN_RESP_TRANS = 3.22;
const double MAIN_RESP_FAC = 0.046;
const double MAIN_RESP_TEMP_TRANS = 20.;
const double MAIN_RESP_TEMP_FAC = 10.;
const double MAIN_RESP_CARB_PROPORT = 1.;

// initial values for trees and grass
const double INIT_MASS_GRASS = 0.001;				// 10 g per m^2
const double INIT_MASS_TREE = 0.001;				// 10 g per plant

// constants to compute moisture of wet biomass for fire, Dennison et al.
const double MOIST_CONST_A = 0.24;
const double MOIST_CONST_B = 1.09;
const double MOIST_CONST_C = 14.2;
const double MOIST_CONST_D = 0.78;
// 
// converts micromol/m^2/s into kg/day/m^2
const double MMSTOKGD_HELPER = 1e-6 * 1e-3 * 44. * (12. / 44.) * (1.544) * 3600.;

// constants needed for data input
const double INPUT_GRID_SIZE = 0.1666667 / 2.;		// half of grid size of input data ( 10min )
const int DATA_NUM = 104;

const int NUM_SIZE_CLASSES = 70;

const int YEARLY_DATA_LENGTH = 54;

const int INPUT_PARAMS = 14;

double const cs_y_vec[42] = { -38.238, -36.372, -34.507, -32.642, -30.777, -28.911, -27.046, -25.181, -23.316, -21.450, -19.585, -17.720, -15.855, -13.989,
		-12.124, -10.259, -8.394, -6.528, -4.663, -2.798, -0.933, 0.933, 2.798, 4.663, 6.528, 8.394, 10.259, 12.124, 13.989, 15.855, 17.720, 19.585, 21.450,
		23.316, 25.181, 27.046, 28.911, 30.777, 32.642, 34.507, 36.372, 38.238 };

double const SFRAC_FINE = 0.34;		// fraction of stem biomass that is fine wood
double const SFRAC_COARSE = 0.25;		// fraction of stem biomass that is coarse wood
double const SFRAC_HEAVY = 0.41;		// fraction of stem biomass that is heavy wood

//const int		TIMESTEP_DAILY_DATA			= 92;		// write daily data every NN days
const int TIMESTEP_DAILY_DATA = 1;		// write daily data every NN days

int const TR_SAV = 0;
int const TR_FOR = 1;

int const GR_C4_SAV = 0;
int const GR_C4_OPN = 1;
int const GR_C4_FOR = 2;
int const GR_C3_SAV = 3;
int const GR_C3_OPN = 4;
int const GR_C3_FOR = 5;

const double PROB_FOREST_TREE = 0.5;          // probability for forest trees (initialization)

