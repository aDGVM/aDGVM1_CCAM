#pragma once

#include "Globals.h"
#include "GlobalTypes.h"
#include "InDataClass.h"

using namespace std;

/**
 * This class contains routines to read input from shortlists.
 */
class clInDataReader
{
public:
	clInDataReader() :
			latitude_(0.), longitude_(0.)
	{
	}
	;
	~clInDataReader()
	{
	}
	;
	clInData getInData(double longitude, double latitude);

private:
	double latitude_;
	double longitude_;
	clInData MyInData_;

	void ReadShortList();
};

clInData clInDataReader::getInData(double longitude, double latitude)
{
	latitude_ = latitude;
	longitude_ = longitude;

	ReadShortList();

	for (size_t i = 0; i < 12; i++)
	{
		MyInData_.rdo_[i] = MyInData_.pwet_[i] * 30.42;
		MyInData_.tmp_day_[i] = MyInData_.tmp_[i] + (MyInData_.tmp_max_[i] - MyInData_.tmp_[i]) / 2.;
	}

	for (size_t i = 0; i < MyInData_.soil_layers_; i++)
		MyInData_.theta_[i] = MyInData_.theta_wp_[i] + 0.25 * (MyInData_.theta_fc_[i] - MyInData_.theta_wp_[i]);

	MyInData_.atm_press_ = 101.325 * pow((293.0 - 0.0065 * MyInData_.elev_) / 293.0, 5.26) * 1000.;
	MyInData_.latitude_ = latitude_;
	MyInData_.longitude_ = longitude_;

	return MyInData_;
}

void clInDataReader::ReadShortList()
{
	string file_name = IN_DATA_HOME + "shortlists/shortlist_" + doubleToString(longitude_) + "_" + doubleToString(latitude_) + "_0.dat";

	ifstream Quelle(file_name.c_str());

	if (!Quelle)
	{
		cerr << "INI ERROR Can't open file: " << file_name << endl;
		cerr << "INI ERROR Maybe the filename is listed in shortlists.slt but the file does not exist." << endl;
		exit(1);
	}

	for (size_t i = 0; i < 12; i++)
		Quelle >> MyInData_.tmp_min_[i];
	for (size_t i = 0; i < 12; i++)
		Quelle >> MyInData_.tmp_max_[i];
	for (size_t i = 0; i < 12; i++)
		Quelle >> MyInData_.tmp_[i];
	for (size_t i = 0; i < 12; i++)
		Quelle >> MyInData_.ralpha_[i];
	for (size_t i = 0; i < 12; i++)
		Quelle >> MyInData_.rbeta_[i];
	for (size_t i = 0; i < 12; i++)
		Quelle >> MyInData_.reh_[i];
	for (size_t i = 0; i < 12; i++)
		Quelle >> MyInData_.pwet_[i];
	for (size_t i = 0; i < 12; i++)
		Quelle >> MyInData_.sun_[i];
	for (size_t i = 0; i < 12; i++)
		Quelle >> MyInData_.wnd_[i];
	for (size_t i = 0; i < 12; i++)
		Quelle >> MyInData_.frost_[i];
	for (size_t i = 0; i < 12; i++)
		MyInData_.frost_[i] = MyInData_.frost_[i] / 31.;

	Quelle >> MyInData_.elev_;
	size_t soil_layers;
	Quelle >> soil_layers;
	assert(soil_layers == MyInData_.soil_layers_);

	for (size_t i = 0; i < MyInData_.soil_layers_; i++)
		Quelle >> MyInData_.theta_wp_[i];
	for (size_t i = 0; i < MyInData_.soil_layers_; i++)
		Quelle >> MyInData_.theta_fc_[i];

	for (size_t i = 0; i < MyInData_.soil_layers_; i++)
		Quelle >> MyInData_.bulk_dens_[i];
	for (size_t i = 0; i < MyInData_.soil_layers_; i++)
		Quelle >> MyInData_.soil_N_[i];
	for (size_t i = 0; i < MyInData_.soil_layers_; i++)
		Quelle >> MyInData_.soil_C_[i];
	for (size_t i = 0; i < MyInData_.soil_layers_; i++)
		Quelle >> MyInData_.thickness_[i];

	for (size_t i = 0; i < MyInData_.soil_layers_; i++)
	{
		if (MyInData_.soil_N_[i] <= 0.000001) MyInData_.soil_N_[i] = 0.01;
	}

	for (size_t i = 0; i < MyInData_.soil_layers_; i++)
	{
		if (MyInData_.soil_C_[i] <= 0.000001) MyInData_.soil_C_[i] = 0.01;
	}

	for (size_t i = 0; i < MyInData_.soil_layers_; i++)
	{
		if (MyInData_.soil_C_[i] > 30000.) MyInData_.soil_C_[i] = 30000.;
	}

	// calculate depth of different soil layers
	MyInData_.depth_[0] = MyInData_.thickness_[0];
	for (size_t i = 1; i < MyInData_.soil_layers_; i++)
		MyInData_.depth_[i] = MyInData_.depth_[i - 1] + MyInData_.thickness_[i];

	if (MyInData_.soil_N_[0] <= -1) //Invalid Gridcell
	{

		for (size_t i = 0; i < 12; i++)
			MyInData_.tmp_min_[i] = -10;
		for (size_t i = 0; i < 12; i++)
			MyInData_.tmp_max_[i] = -10;
		for (size_t i = 0; i < 12; i++)
			MyInData_.tmp_[i] = -10;
		for (size_t i = 0; i < 12; i++)
			MyInData_.ralpha_[i] = -10;
		for (size_t i = 0; i < 12; i++)
			MyInData_.rbeta_[i] = -10;
		for (size_t i = 0; i < 12; i++)
			MyInData_.reh_[i] = -10;
		for (size_t i = 0; i < 12; i++)
			MyInData_.pwet_[i] = -10;
		for (size_t i = 0; i < 12; i++)
			MyInData_.sun_[i] = -10;
		for (size_t i = 0; i < 12; i++)
			MyInData_.wnd_[i] = -10;
		for (size_t i = 0; i < 12; i++)
			MyInData_.frost_[i] = -10;

		MyInData_.atm_press_ = -10;
		MyInData_.elev_ = -10;

		for (size_t i = 0; i < MyInData_.soil_layers_; i++)
			MyInData_.theta_wp_[i] = -10;
		for (size_t i = 0; i < MyInData_.soil_layers_; i++)
			MyInData_.theta_fc_[i] = -10;
		for (size_t i = 0; i < MyInData_.soil_layers_; i++)
			MyInData_.bulk_dens_[i] = -10;
		for (size_t i = 0; i < MyInData_.soil_layers_; i++)
			MyInData_.soil_N_[i] = -10;
		for (size_t i = 0; i < MyInData_.soil_layers_; i++)
			MyInData_.soil_C_[i] = -10;
	}

	Quelle.close();

	return;
}

