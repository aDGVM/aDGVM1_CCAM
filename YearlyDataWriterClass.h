#pragma once

// This class provides routines to write yearly output data into the
// files YearlyData.***.dat in the output directory. Output data must
// be specified in the file yearly_output_variables.h

class clYearlyDataWriter
{
private:
	int data_length_;
	int count_years_;
	int i_scen_;
	int with_fire_;
	int d_seed_;
	double *data_set_;

public:
	clYearlyDataWriter(int data_length);
	~clYearlyDataWriter();
	void setYearlyData(double *data_set);
	void printYearlyDataToConsole();
	void printYearlyDataToFile(string output_directory);
};

// -----------------------------------------------------------------
// -- constructor allocates vector data_set and                   --
// -- initializes all values with zero                            --
// -----------------------------------------------------------------
clYearlyDataWriter::clYearlyDataWriter(int data_length)
{
	data_length_ = data_length;
	data_set_ = (double*) malloc(data_length_ * sizeof(double));

	for (int i = 0; i < data_length_; i++)
		data_set_[i] = -9999;

	count_years_ = -9999;
	i_scen_ = -9999;
	with_fire_ = -9999;
	d_seed_ = -9999;
}

// -----------------------------------------------------------------
// -- destructor, deletes vector data_set_                        --
// -----------------------------------------------------------------
clYearlyDataWriter::~clYearlyDataWriter()
{
	free(data_set_);
}

// -----------------------------------------------------------------
// -- Sets all variables                                          --
// -----------------------------------------------------------------
void clYearlyDataWriter::setYearlyData(double *data_set)
{
	for (int i = 0; i < data_length_; i++)
		data_set_[i] = data_set[i];

	count_years_ = (int) data_set_[2];
	i_scen_ = (int) data_set_[3];
	with_fire_ = (int) data_set_[4];
	d_seed_ = (int) data_set_[5];

	return;
}

// -----------------------------------------------------------------
// -- prints dataset to console                                   --
// -----------------------------------------------------------------
void clYearlyDataWriter::printYearlyDataToConsole()
{
	for (int i = 0; i < data_length_; i++)
		cout << setw(14) << data_set_[i];

	cout << endl;

	return;
}

void clYearlyDataWriter::printYearlyDataToFile(string output_directory)
{
	char file_name[250];
	sprintf(file_name, "%sYearlyData.%i.%i.%i.dat", output_directory.c_str(), count_years_, i_scen_, with_fire_);

	ofstream YearlyData(file_name, ios::app);
	if (!file_name)
	{
		cerr << "Can't open file: " << file_name << endl;
		exit(1);
	}

	for (int i = 0; i < data_length_; i++)
		YearlyData << setw(14) << data_set_[i];

	YearlyData << endl;

	YearlyData.close();

	return;
}

