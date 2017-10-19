# ifndef _IBC
# define _IBC

# include <iostream>
# include <armadillo>
# include "Para.h"

using namespace std;
using namespace arma;

class IBC :public Para
{
public:
	// ini_type is the type of the initial conditions
	string ini_type;
	// bnd_type is the type of boundary conditions
	string bnd_type;
	// q_type is the choices of different schemes
	string q_type;
	// flux_type determines the type of flux
	string flux_type;

	IBC(string ini_typee = "", string bnd_typee = "", string q_typee = "", string flux_typee = "");

	void set_IBC(string ini_typee, string bnd_typee);
	void set_scheme(string q_typee);
	void set_flux(string flux_type);
	void IBC_cetr(void);
};



class IB_sin
{
public:
	IB_sin();
	void eva_IB(Para * pP);
};


class IB_shock
{
public:
	IB_shock(void);
	void eva_IB(Para *pP);
};

class IB_cos
{
public:
	IB_cos(void);
	void eva_IB(Para* pP);
	int get_closest_index( rowvec & x, double xm);
};
# endif
