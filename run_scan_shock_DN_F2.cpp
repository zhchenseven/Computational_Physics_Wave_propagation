# include <iostream>
# include <armadillo>

# include "Data_process.h"

using namespace std;
using namespace arma;

int rssD2(int argv, char ** argc)
//int main(int argv, char ** argc)
{
	Data_process Project;
	rowvec nu_arr(4);
	nu_arr << 1 << 0.75 << 0.5 << 0.25 << endr;
	int N = 32;
	double x_low = 0;
	double x_sup = 1;
	int M = 30;
	double a = 1;
	double nu = 1;
	double NT = 0.5;
	Project.set_para(N, x_low,x_sup, M, a, nu,  NT);
	string ini_type = "shock";
	string bnd_type = "DN";
	Project.set_IBC(ini_type, bnd_type);
	string q_type = "UW1";
	Project.set_scheme(q_type);
	string flux_type = "F2";
	Project.set_flux(flux_type);

	string vis_name = "shock_DN_F2";
	int index = 0;
	string data_type = "dat";
	string main_fd_name = "shock_DN_F2_scan_T_0p5";
	string size_sufx = "size";
	Project.set_process_ctrl(vis_name, index, data_type, size_sufx);
	Project.set_main_fd(main_fd_name);


	int j;
	for (j = 0; j < 4; j++)
	{
		Project.creat_main_fd_ctrl = 1;
		vis_name = ini_type + "_" + bnd_type + "_nu_" + Project.dn2s(nu_arr(j));
		Project.set_vis_name(vis_name);
		cout << "vis_name is " << Project.vis_name << endl;
		Project.nu = nu_arr(j);
		if (j != 0)
			Project.creat_main_fd_ctrl = 0;
		Project.exe_cetr();
		Project.process();
		Project.write_log();
		}


	system("pause");
	return 0;
}