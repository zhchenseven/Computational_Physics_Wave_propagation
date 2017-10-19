# include <iostream>
# include <armadillo>

# include "Data_process.h"

using namespace std;
using namespace arma;

int rLd(int argv, char ** argc)
{
	Data_process Project;
	int N = 32;
	double x_low = 0;
	double x_sup = 1;
	int M = 30;
	double a = 1;
	double nu = 1;
	double NT = 1.5;
	Project.set_para(N, x_low,x_sup, M, a, nu,  NT);
	string ini_type = "sin";
	string bnd_type = "PD";
	Project.set_IBC(ini_type, bnd_type);
	string q_type = "UW1";
	Project.set_scheme(q_type);
	string flux_type = "F1";
	Project.set_flux(flux_type);

	string vis_name = "LF";
	int index = 0;
	string data_type = "dat";
	string size_sufx = "size";
	Project.set_process_ctrl(vis_name, index, data_type, size_sufx);

	Project.exe_cetr();
	Project.process();
	Project.write_log();

	system("pause");
	return 0;
}