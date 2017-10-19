# include <iostream>
# include <armadillo>

# include "Data_process.h"

using namespace std;
using namespace arma;

int rcP2Ud(int argv, char ** argc)
{
	Data_process Project;
	int N = 32;
	double x_low = 0;
	double x_sup = 1;
	int M = 30;
	double a = 1;
	double nu = 0.25;
	double NT = 1.5;
	Project.set_para(N, x_low,x_sup, M, a, nu,  NT);
	string ini_type = "cos";
	string bnd_type = "PD";
	Project.set_IBC(ini_type, bnd_type);
	string q_type = "UW1";
	Project.set_scheme(q_type);
	string flux_type = "F2";
	Project.set_flux(flux_type);

	string vis_name = "cos_PD_F2_UW1";
	int index = 0;
	string data_type = "dat";
	string size_sufx = "size";
	Project.set_process_ctrl(vis_name, index, data_type, size_sufx);


	double xi = 0.25;
	double d = 0.1;

	Project.set_para_extra(xi, d);

	Project.exe_cetr();
	Project.process();
	Project.write_log();

	system("pause");
	return 0;
}