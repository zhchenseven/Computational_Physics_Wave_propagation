#ifndef _PARA
# define _PARA

# include <iostream>
# include <armadillo>

using namespace std;
using namespace arma;

class Para
{
public:
	//N is the number of gridding points except two end points
	int N;
	//x_low is the lower limit of x
	double x_low;
	// x_sup is the upper limit of x
	double x_sup;
	//M controls the number of tracking time, which includes the initial time and the final time, we demand M>2
	int M;
	// a is the velocity
	double a;
	// T specifies the final time point
	double T;
	// nu is the CFL number 
	double nu;

	// uo is the wavefunction at initial time
	rowvec uo;

	// un 
	rowvec un;

	// dt is the time spacing
	double dx;
	// dx is the spacing
	double dt;
	// x is the range of interest
	rowvec x;

	// u_track is the matrix to record data
	mat u_track;

	// ind_track specifies the index to  record data
	uvec ind_track;

	// N_itr is the number of iterations to undergo
	int N_itr;

	//M_fac is the actual number of time track
	int M_fac;

	// NT specifies number of periods we run
	double NT;

	// u_exact is the exact solution of the wave equation corresponding to the recording time
	mat u_exact;

	// ui records the initial wave packet for the later exact solutions
	rowvec ui;

	//xi indicates the center of the third initial condition
	double xi;

	// d indicates the spatial period of the third condition
	double d;

	// uf_exact indicates the exact solution of the final time point
	rowvec uf_exact;

	// t_out is the output time points
	rowvec t_out;

	// err is the error of L1 error
	double err;


	Para(int NN = 32, double x_loww = 0, double x_supp = 0, int MM = 50, double aa = 1, double TT = 1, double nuu = 1, rowvec uoo=rowvec(), rowvec unn=rowvec(),
	double dxx=0,double dtt=0,rowvec xx=rowvec() ,mat u_trackk=mat(), uvec ind_trackk=uvec(), int N_itrr =0, int M_facc=0, double NTT=0, mat u_exactt=mat(),rowvec uii=rowvec(),
	double xii=0, double dd=0, rowvec uf_exactt=rowvec(),rowvec t_outt=rowvec(), double errr=0);

	void set_para(int NN, double x_loww, double x_supp, int MM, double aa,  double nuu, double NTT);

	void set_para_extra(double xii, double dd);

	
};

# endif
