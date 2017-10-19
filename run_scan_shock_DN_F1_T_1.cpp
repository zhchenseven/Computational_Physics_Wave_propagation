# include <iostream>
# include <armadillo>

# include "Data_process.h"

using namespace std;
using namespace arma;

int rssD1T1(int argv, char ** argc)
//int main(int argv, char ** argc)
{
	Data_process Project;
	int N = 32;
	rowvec nu_arr(4);
	string q_arr[4] = { "LF","UW1","LW","MD" };
	nu_arr << 1 << 0.75 << 0.5 << 0.25 << endr;
	double x_low = 0;
	double x_sup = 1;
	int M = 30;
	double a = 1;
	double nu = 1;
	double NT = 1;
	Project.set_para(N, x_low,x_sup, M, a, nu,  NT);
	string ini_type = "shock";
	string bnd_type = "DN";
	Project.set_IBC(ini_type, bnd_type);
	string q_type = "UW1";
	Project.set_scheme(q_type);
	string flux_type = "F1";
	Project.set_flux(flux_type);

	string vis_name = "a";
	int index = 0;
	string main_fd_name = "shock_DN_F1_scan_T_1";
	string data_type = "dat";
	string size_sufx = "size";
	Project.set_process_ctrl(vis_name, index, data_type, size_sufx);
	Project.set_main_fd(main_fd_name);

	int i, j;
	for (i = 0; i < 4; i++)
	{
		Project.q_type = q_arr[i];
		for (j = 0; j < 4; j++)
		{
			Project.creat_main_fd_ctrl = 1;
			vis_name = ini_type + "_" + bnd_type + "_" + flux_type + "_" + q_arr[i] + "_nu_" + Project.dn2s(nu_arr(j));
			Project.set_vis_name(vis_name);
			cout << "vis_name is " << Project.vis_name << endl;
			Project.nu = nu_arr(j);
			if (i != 0 || j != 0)
				Project.creat_main_fd_ctrl = 0;
			Project.exe_cetr();
			Project.process();
			Project.write_log();
		}
	}



	system("pause");
	return 0;
}