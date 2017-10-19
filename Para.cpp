# include <iostream>
# include <armadillo>
# include "Para.h"


using namespace std;
using namespace arma;

Para::Para(int NN , double x_loww , double x_supp , int MM , double aa , double TT , double nuu , rowvec uoo , rowvec unn ,
	double dxx, double dtt, rowvec xx , mat u_trackk, uvec ind_trackk , int N_itrr, int M_facc, double NTT, mat u_exactt,rowvec uii,
	double xii , double dd, rowvec uf_exactt, rowvec t_outt, double errr )
{
	N = NN;
	x_low = x_loww;
	x_sup = x_supp;
	M = MM;
	a = aa;
	T = TT;
	nu = nuu;
	uo = uoo;
	un = unn;
	dx = dxx;
	dt = dtt;
	x = xx;
	u_track = u_trackk;
	ind_trackk = ind_track;
	N_itr = N_itrr;
	M_fac = M_facc;
	NT = NTT;
	u_exact = u_exactt;
	ui = uii;
	xi = xii;
	d = dd;
	uf_exact = uf_exactt;
	t_out = t_outt;
	err = errr;
	cout << "A class of Para is constructed." << endl << endl;
}


void Para::set_para(int NN, double x_loww, double x_supp, int MM, double aa, double nuu, double NTT )
{
	N = NN;
	x_low = x_loww;
	x_sup = x_supp;
	M = MM;
	a = aa;
	NT = NTT;
	T = NTT*(x_sup-x_low)/a;
	nu = nuu;
	x = linspace<rowvec>(x_low, x_sup, N + 2);
	dx = x(1) - x(0);
	dt = nu*dx / a;
	uo = rowvec(N + 2, fill::zeros);
	un = rowvec(N + 2, fill::zeros);
	N_itr = int(T / dt);
}


void Para::set_para_extra(double xii, double dd)
{
	xi = xii;
	d = dd;
}