# include <iostream>
# include <armadillo>

# include "IBC.h"

using namespace std;
using namespace arma;

IBC::IBC(string ini_typee , string bnd_typee , string q_typee , string flux_typee )
{
	ini_type = ini_typee;
	bnd_type = bnd_typee;
	q_type = q_typee;
	flux_type = flux_typee;
	Para();
	cout << "A class of IBC is constructed." << endl << endl;
}



void IBC::set_IBC(string ini_typee, string bnd_typee)
{
	ini_type = ini_typee;
	bnd_type = bnd_typee;
}

void IBC::set_scheme(string q_typee)
{
	q_type = q_typee;
}

void IBC::set_flux(string flux_typee)
{
	flux_type = flux_typee;
}

void IBC::IBC_cetr(void)
{
	dt = nu*dx / a;
	N_itr = int(T / dt);
	if (ini_type == "sin")
	{
		Para * pP = (Para *)this;
		IB_sin I;
		I.eva_IB(pP);
	}
	else if (ini_type == "shock")
	{

			Para * pP = (Para *)this;
			IB_shock I;
			I.eva_IB(pP);

	}
	else if (ini_type == "cos")
	{
		Para * pP = (Para *)this;
		IB_cos I;
		I.eva_IB(pP);
	}
	else
	{
		cout << "Your specified initial condition type " << ini_type << " is not available yet." << endl;
		exit(1);
	}
}


IB_sin::IB_sin()
{
	cout << "A class of sinunoidal initial condition is constructed." << endl << endl;
}

void IB_sin::eva_IB(Para * pP)
{
	double pi = datum::pi;
	pP->uo= rowvec(pP->x.n_elem, fill::zeros);
	pP->uo = sin(2 * pi*pP->x);
	pP->ui = pP->uo;
	pP->un = rowvec(pP->x.n_elem, fill::zeros);
}

IB_shock::IB_shock()
{
	cout << "A class of shock discontinuity initial condition is constructed." << endl << endl;
}

void IB_shock::eva_IB(Para *pP)
{
	int n_bd = int((pP->N-2)/ 4);
	pP->uo = rowvec(pP->x.n_elem, fill::zeros);
	pP->uo(span(0, n_bd)) = rowvec(n_bd+1,fill::ones);
	pP->ui = pP->uo;
	pP->un = rowvec(pP->x.n_elem, fill::zeros);

}

IB_cos::IB_cos(void)
{
	cout << "A class of cos initial condition is constructed." << endl << endl;
}

void IB_cos::eva_IB(Para* pP)
{
	double x_min, x_max,pi=datum::pi;
	pP->uo = rowvec(pP->x.n_elem, fill::zeros);
	x_min = pP->xi - pP->d;
	x_max = pP->xi + pP->d;
	int ind_min, ind_max;
	ind_min = get_closest_index(pP->x, x_min);
	ind_max = get_closest_index(pP->x, x_max);
	rowvec id(pP->N + 2, fill::zeros);
	id(span(ind_min, ind_max)) = rowvec(ind_max-ind_min+1,fill::ones);
	pP->ui = arma::abs(cos(pi*(pP->x - pP->xi) / (2 * pP->d)));
	pP->ui = pP->ui%id;
	pP->uo = pP->ui;
	pP->un = rowvec(pP->x.n_elem, fill::zeros);

}

int IB_cos::get_closest_index( rowvec & x, double x_ref)
{
	rowvec x_abs = arma::abs(x - x_ref);
	int ind_c;
	ind_c = x_abs.index_min();
	return ind_c;
}