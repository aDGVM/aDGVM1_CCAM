#pragma once

// Global parameters for the Yasso soil model

// the concentrations of different compounds must sum to one
double const YASSO_C_NWL_EXT = 0.27;    // concentration of compound ext in litter type nwl
double const YASSO_C_NWL_CEL = 0.51;    // concentration of compound cel in litter type nwl
double const YASSO_C_NWL_LIG = 0.22;    // concentration of compound lig in litter type nwl

double const YASSO_C_FWL_EXT = 0.03;    // concentration of compound ext in litter type fwl
double const YASSO_C_FWL_CEL = 0.65;    // concentration of compound cel in litter type fwl
double const YASSO_C_FWL_LIG = 0.32;    // concentration of compound lig in litter type fwl

double const YASSO_C_CWL_EXT = 0.03;    // concentration of compound ext in litter type cwl
double const YASSO_C_CWL_CEL = 0.69;    // concentration of compound cel in litter type cwl
double const YASSO_C_CWL_LIG = 0.28;    // concentration of compound lig in litter type cwl

// docomposition rates, scaled from rate per year to rate per day
double const YASSO_A_CWL = 0.077 / 365.;   // exposure of coarse woody litter to microbial decomposition
double const YASSO_A_FWL = 0.54 / 365.;    // exposure of fine woody litter to microbial decomposition
double const YASSO_K_EXT = 0.82 / 365.;    // decomposition rate of ext
double const YASSO_K_CEL = 0.30 / 365.;    // decomposition rate of cel
double const YASSO_K_LIG = 0.22 / 365.;    // decomposition rate of lig
double const YASSO_K_HUM1 = 0.012 / 365.;   // decomposition rate of hum1 (faster humus)
double const YASSO_K_HUM2 = 0.0012 / 365.;  // decomposition rate of hum2 (slower humus)

double const YASSO_P_EXT = 0.2;     // proportion of mass decomposed in compartment ext
double const YASSO_P_CEL = 0.2;     // proportion of mass decomposed in compartment cel
double const YASSO_P_LIG = 0.2;     // proportion of mass decomposed in compartment lig
double const YASSO_P_HUM1 = 0.2;     // proportion of mass decomposed in compartment hum1

double const YASSO_BETA = 0.105;   // parameter for effects of temperature
double const YASSO_GAMMA = 0.00274; // parameter for effects of drought

double const YASSO_S_HUM = 0.9;     // sensitivity of humus (<1)

double const YASSO_T0 = 20.;     // temperature in standard conditions
double const YASSO_D0 = 20.;     // summer drought in standard conditions

