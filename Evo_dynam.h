# ifndef _EVO_DYNAM
# define _EVO_DYNAM

# include <iostream>
# include <armadillo>
# include "IBC.h"

using namespace std;
using namespace arma;


class Evo_dynam :public IBC
{
public:
	// qnu is the constant of q(nu)*dx/(2*dt)
	double qnu;
	//txr is dt/dx
	double txr;
	//one_nu2 is (1-nu)/2
	double one_nu2;
	Evo_dynam(double qnuu=0, double txrr=0, double one_nu22=0);
	void evo_exe_cetr(void);
	// this is the scheme of the first flux function and the periodic condition
	void evo_F1_perd(void);
	void evo_F1_DN(void);
	void evo_F2_DN(void);
	double  q_LF(double nu);
	double q_UW1(double nu);
	double q_LW(double nu);
	double q_MD(double nu);

	typedef double(Evo_dynam::*p_q)(double);
	p_q retn_q_func(void);

	void get_track_point(void);

	double F1(rowvec & uu, int j1, int j2);

	double f(double xo);

	rowvec get_exact_sol(int j);

	double B(double x, double y);

	double max(double x, double y);
	double min(double x, double y);

	double sign(double x);
	double F2(rowvec & uu, int j);

	void evo_F2_perd(void);

	double F2(rowvec & uu, int j1,int j2, int j3);

	double eva_ini(double xo,double to=0);
};

# endif
