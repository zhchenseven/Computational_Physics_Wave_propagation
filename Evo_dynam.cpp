# include <iostream>
# include <armadillo>

# include "Evo_dynam.h"

using namespace std;
using namespace arma;

enum q_type_code
{
	LF,
	UW1,
	LW,
	MD

};


q_type_code hash_q_type(string const & s)
{
	if (s == "LF")
		return LF;
	if (s == "UW1")
		return UW1;
	if (s == "LW")
		return LW;
	if (s == "MD")
		return MD;
	else
	{
		cout << "Your specified scheme type " << s << " is not available yet." << endl;
		exit(1);
	}
}



Evo_dynam::Evo_dynam(double qnuu , double txrr,double one_nu22)
{
	qnu = qnuu;
	txr = txrr;
	one_nu2 = one_nu22;
	IBC();
	cout << "A class of Evo_dynam is constructed." << endl << endl;
}


void Evo_dynam::evo_exe_cetr(void)
{
	p_q pq = retn_q_func();
	cout << "nu=" << nu << ",dx=" << dx << ",dt=" << dt << endl;
	qnu = (this->*pq)(nu)*dx/(2*dt);
	get_track_point();
	if (bnd_type == "PD" && flux_type == "F1")
		evo_F1_perd();
	else if (bnd_type == "DN" && flux_type == "F1")
		evo_F1_DN();
	else if (bnd_type == "DN" && flux_type == "F2")
		evo_F2_DN();
	else if (bnd_type == "PD" && flux_type == "F2")
		evo_F2_perd();
	else
	{
		cout << "Your specified flux type " << flux_type << " and boundary condition "<<bnd_type<<  " is not available yet." << endl;
		exit(1);
	}
}



typedef double(Evo_dynam::*p_q)(double);
p_q Evo_dynam::retn_q_func(void)
{
	switch (hash_q_type(q_type))
	{
	case LF:
	{
		p_q p = &Evo_dynam::q_LF;
		return p;
	}

	case UW1:
	{
		p_q p = &Evo_dynam::q_UW1;
		return p;
	}

	case LW:
	{
		p_q p = &Evo_dynam::q_LW;
		return p;
	}

	case MD:
	{
		p_q p = &Evo_dynam::q_MD;
		return p;
	}

	default:
	{
		cout << "Your specified scheme type " << q_type << " is not available yet." << endl;
		exit(1);
		break;
		break;
	}
	}
}




double  Evo_dynam::q_LF(double nu)
{
	return 1.0;
}

double Evo_dynam::q_UW1(double nu)
{
	return abs(nu);
}


double Evo_dynam::q_LW(double nu)
{
	return pow(nu, 2);
}

double Evo_dynam::q_MD(double nu)
{
	return 1.0 / 3 + 2 * pow(nu, 2) / 3;
}


void Evo_dynam::evo_F1_perd(void)
{
	int n,j;
	int n_flag = 1;
	double Fj1,Fj2;
	txr = dt / dx;
	double t_elasped;
	wall_clock timer;
	t_out(0) = 0;
	u_track.row(0) = uo;
	u_exact.row(0) = uo;
	cout << "M= " << M << ",M_fac= " << M_fac<<endl;
	timer.tic();
	for (n = 0; n < N_itr; n++)
	{
		cout << "N_itr=" << N_itr << endl << "n=" << n << endl;
		for (j = 0; j < N + 1; j++)
		{
			if (j == 0)
			{
				Fj1 = F1(uo, N, 0);
				Fj2 = F1(uo, 0, 1);
			}
			else
			{
				Fj1 = F1(uo, j-1, j);
				Fj2 = F1(uo, j, j+1);
			}
			un(j) = uo(j) - txr*(Fj2 - Fj1);
		}
		un(N + 1) = un(0);
		if (n == ind_track(n_flag))
		{
			t_out(n_flag) = (n + 1)*dt;
			u_track.row(n_flag) = un;
			u_exact.row(n_flag) = get_exact_sol(n + 1);
			n_flag++;
		}
		uo = un;
	}
	//uf_exact = u_exact.row(n_flag - 1);
	uf_exact = u_exact.row(u_exact.n_rows - 1);
	t_elasped = timer.toc();
	cout << "The F1 periodic boundary condition equation takes " << t_elasped << " s in C++." << endl;
}


void Evo_dynam::get_track_point(void)
{
	if (N_itr <= M - 1)
	{
		M_fac = N_itr+1;
		ind_track = regspace<uvec>(0, M_fac-1);
	}
	else
	{
		M_fac = M;
		int q = (N_itr - 1) / (M - 2);
		int r = N_itr % (M - 1);
		int K = N_itr - 1 + (2 - M)*q;
		int L = M - 2 - K;
		cout << "q= " << q << ",r= " << r << ",K=" << K << ",L=" << L << endl;
		ind_track = uvec(M_fac, fill::zeros);
		int i;
		for (i = 1; i <= M_fac - 1; i++)
		{
			if (i == 1)
				ind_track(i) = 1;
			else if (i <= 1 + K)
			{
				ind_track(i) = ind_track(i - 1) + q + 1;
			}
			else
				ind_track(i) = ind_track(i - 1) + q;

		}
	}
	cout << "M_fac=" << M_fac << endl<<"N_itr="<<N_itr << "ind_track" << ind_track << endl;
	ind_track -= 1;
	ind_track(0) = 0;
	t_out = rowvec(M_fac, fill::zeros);
	u_track = mat(M_fac, N + 2, fill::zeros);
	u_exact = mat(M_fac, N + 2, fill::zeros);
}


inline double Evo_dynam::f(double xo)
{
	return a*xo;
}

inline double Evo_dynam::F1(rowvec & uu, int j1, int j2)
{
	// j1 j2 typically maps to j,j+1
	double F;
	F = (f(uu(j1)) + f(uu(j2))) / 2 - qnu*(uu(j2) - uu(j1));
	return F;
}


//rowvec Evo_dynam::get_exact_sol(int j)
//{
//	int N_all = N + 2,i,pi;
//	rowvec un_exact = rowvec(N_all, fill::zeros);
//	if (bnd_type == "periodic")
//	{
//		for (i = 0; i < N_all; i++)
//		{
//			pi = i - j;
//			while (pi < 0)
//				pi += N_all;
//			un_exact(i) = ui(pi);
//		}
//	}
//	else
//	{
//		un_exact = rowvec(N_all, fill::ones);
//		int n_bd= int((N - 2) / 4)+j;
//		for (i = 0; i < N_all; i++)
//		{
//			//pi = i - j;
//			//if(pi < 0)
//			//	pi += N_all;
//			//if (pi > n_bd)
//			//	un_exact(i) = 0;
//			if (i > n_bd)
//				un_exact(i) = 0;
//		}
//	}
//	return un_exact;
//}


rowvec Evo_dynam::get_exact_sol(int j)
{
	int N_all = N + 2, i, pi;
	rowvec un_exact = rowvec(N_all, fill::zeros);
	for (i = 0; i < N_all; i++)
	{
		un_exact(i) = eva_ini(x(i) - a*j*dt);
	}
	return un_exact;
}

void Evo_dynam::evo_F1_DN(void)
{
	int n, j;
	int n_flag = 1;
	double Fj1, Fj2;
	txr = dt / dx;
	double t_elasped;
	wall_clock timer;
	t_out(0) = 0;
	u_track.row(0) = uo;
	u_exact.row(0) = uo;
	cout << "M= " << M << ",M_fac= " << M_fac<<", N_itr= "<<N_itr << endl<<"ind_track"<<ind_track<< endl;
	timer.tic();
	for (n = 0; n < N_itr; n++)
	{
		for (j = 0; j < N + 1; j++)
		{
			if (j == 0)
			{
				un(j) = 1;
				//Fj1 = F1(uo, N, 0);
				//Fj2 = F1(uo, 0, 1);
			}
			else
			{
				Fj1 = F1(uo, j - 1, j);
				Fj2 = F1(uo, j, j + 1);
				un(j) = uo(j) - txr*(Fj2 - Fj1);
			}
			
		}
		un(N + 1) = un(N);
		//cout << "n= " << n << ",n_flag= "<<n_flag<<",size ind_track= "<<ind_track.n_elem<< endl;
		if (n == ind_track(n_flag))
		{
			t_out(n_flag) = (n + 1)*dt;
			u_track.row(n_flag) = un;
			u_exact.row(n_flag) = get_exact_sol(n + 1);
			n_flag++;
		}
		uo = un;
	}
	uf_exact = u_exact.row(n_flag-1);
	//uf_exact = get_exact_sol(n);
	//u_exact.row(u_exact.n_rows-1) = uf_exact;
	t_elasped = timer.toc();
	cout << "The F1 dirichlet and neumann condition takes " << t_elasped << " s in C++." << endl;
}

double Evo_dynam::B(double x, double y)
{
	double f;
	if (x*y <= 0)
		f = 0;
	else 
	{
		double rxy = x / y;
		if (rxy <= 2 && rxy >= 0.5)
			f = sign(x)*max(abs(x), abs(y));
		else
			f = 2*sign(x)*min(abs(x), abs(y));
	}
	return f;
}



inline double Evo_dynam::max(double x, double y)
{
	return (x > y) ? x : y;
}

inline double Evo_dynam::min(double x, double y)
{
	return (x < y) ? x : y;
}

inline double Evo_dynam::sign(double x)
{
	if (x > 0)
		return 1.0;
	else if (x == 0)
		return .0;
	else
		return -1.0;
}

void Evo_dynam::evo_F2_DN(void)
{
	int n, j;
	int n_flag = 1;
	double Fj1, Fj2;
	txr = dt / dx;
	double t_elasped;
	wall_clock timer;
	one_nu2 = (1 - nu) / 2;
	u_track.row(0) = uo;
	u_exact.row(0) = uo;
	t_out(0) = 0;
	timer.tic();
	for (n = 0; n < N_itr; n++)
	{
		for (j = 0; j < N + 1; j++)
		{
			if (j == 0 || j==1)
			{
				un(j) = 1;
			}
			else
			{
				Fj1 = F2(uo, j - 1);
				Fj2 = F2(uo, j);
				un(j) = uo(j) - txr*(Fj2 - Fj1);
			}
			
		}
		un(N + 1) = un(N);
		if (n == ind_track(n_flag))
		{
			t_out(n_flag) = (n + 1)*dt;
			u_track.row(n_flag) = un;
			u_exact.row(n_flag) = get_exact_sol(n + 1);
			n_flag++;
		}
		uo = un;
	}
	uf_exact = u_exact.row(n_flag - 1);
	//uf_exact = get_exact_sol(n);
	//u_exact.row(u_exact.n_rows - 1) = uf_exact;
	t_elasped = timer.toc();
	cout << "The F2 periodic boundary condition equation takes " << t_elasped << " s in C++." << endl;
}

double Evo_dynam::F2(rowvec & uu, int j)
{
	return f(uu(j) + one_nu2*B(uu(j + 1) - uu(j), uu(j) - uu(j-1)));
}


void Evo_dynam::evo_F2_perd(void)
{
	int n, j;
	int n_flag = 1;
	double Fj1, Fj2;
	txr = dt / dx;
	double t_elasped;
	wall_clock timer;
	one_nu2 = (1 - nu) / 2;
	t_out(0) = 0;
	u_track.row(0) = uo;
	u_exact.row(0) = uo;
	cout << "M= " << M << ",M_fac= " << M_fac << endl;
	timer.tic();
	for (n = 0; n < N_itr; n++)
	{
		for (j = 0; j < N + 1; j++)
		{
			if (j == 0)
			{
				Fj2 = F2(uo, 1, 0,N);
				Fj1 = F2(uo, 0, N,N-1);
			}
			else if (j==1)
			{
				Fj2 = F2(uo, 2, 1,0);
				Fj1 = F2(uo, 1, 0,N);
			}
			else
			{
				Fj1 = F2(uo, j - 1);
				Fj2 = F2(uo, j);
			}
			un(j) = uo(j) - txr*(Fj2 - Fj1);
		}
		un(N + 1) = un(0);
		if (n == ind_track(n_flag))
		{
			t_out(n_flag) = (n + 1)*dt;
			u_track.row(n_flag) = un;
			u_exact.row(n_flag) = get_exact_sol(n + 1);
			n_flag++;
		}
		uo = un;
	}
	uf_exact = u_exact.row(n_flag - 1);
	t_elasped = timer.toc();
	cout << "The F2 periodic boundary condition equation takes " << t_elasped << " s in C++." << endl;
}


double Evo_dynam::F2(rowvec & uu, int j1, int j2, int j3)
{
	// j1 j2 j3 typically corresponds to j+1,j,j-1
	return f(uu(j2) + one_nu2*B(uu(j1) - uu(j2), uu(j2) - uu(j3)));
}

double Evo_dynam::eva_ini(double xo, double to )
{
	double x_range = x_sup - x_low,pi=datum::pi;
	while ((xo<x_low || xo>x_sup) && ini_type !="shock")
	{
		if (xo < x_low)
			xo += x_range;
		else
			xo -= x_range;
	}
	if (ini_type == "sin")
		return sin(2 * pi*xo);
	else if (ini_type == "shock")
	{
		double x_thre = x_low + x_range / 4 ;
		if (xo <= x_thre)
			return 1;
		else
			return 0;
	}
	else if (ini_type == "cos")
	{
		if (xo <= xi + d && xo >= xi - d)
			return pow(cos(pi*(xo - xi) / (2 * d)), 2);
		else
			return 0;
	}
}