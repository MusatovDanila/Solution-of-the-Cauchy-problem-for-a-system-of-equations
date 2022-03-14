#include<iostream>
#include<vector>	
#include<cmath>
#include<Eigen/Dense>
#include<Eigen/LU> 

using std::cout;
using std::vector;
template <typename T>
using v_func_v_x = vector < double(*)(T)>;
template <typename T>
using m_func_x = vector < vector <T>>;
using v_i = vector<int>;
using v_d = vector<double>;
using EVX = Eigen::VectorXf;
using EMX = Eigen::MatrixXf; 
using V_EMX = vector<Eigen::MatrixXf>;

template <typename T>
double max_delt(const T& x, const T& x_1) {
	double max = -1;
	for (unsigned int i = 0; i < x.size(); ++i) {
		max = (fabs(x[i] - x_1[i]) > max) ? fabs(x[i] - x_1[i]) : max;
	}
	return max;
}

EVX conv_v_e(const v_d& x) {
	EVX e_v_X(x.size());
	for (unsigned int i = 0; i < x.size(); ++i) {
		e_v_X(i) = x[i];
	}
	return e_v_X;
}

v_d conv_e_v(const EVX& x) {
	v_d x_v(x.size());
	for (int i = 0; i < x.size(); ++i) {
		x_v[i]=x[i];
	}
	return x_v;
}

template<typename T>
EMX W_completion(const m_func_x< double(*)(EVX x)>& d_fi, const T& e_v_x) {
	EMX W(e_v_x.size(), e_v_x.size());
	for (int i = 0; i < e_v_x.size(); ++i) {
		for (int j = 0; j < e_v_x.size(); ++j) {
			W(i, j) = d_fi[i][j](e_v_x);
		}
	}
	return W;
}


EVX Fx(const EVX& evx, const v_func_v_x<EVX>& f) {
	EVX FX(f.size());
	for (unsigned int i = 0; i < f.size(); ++i) {
		FX[i] = f[i](evx);
	}
	return FX;
}

v_d converter_e_v(const EVX& x) {
	v_d b_x(x.size());
	for (int i = 0; i < x.size(); i++) b_x[i]=x[i];
	return b_x;
}

template < typename F, typename X>
double max_f_x(const F& f_x, const X& x) {
	double max = -1;
	for (unsigned int i = 0; i < x.size(); ++i) {
		max = (fabs(f_x[i](x)) > max) ? fabs(f_x[i](x)) : max;
	}
	return max;
}

template <typename F, typename X>
void print_disc(const F& f, const X& x) {
	cout << "\n\nx[0] = "<< x[0] << "      x[1] = " << x[1] <<"        f_0(x) = " <<f[0](x) << "        f_1(x) = " << f[1](x) <<"        discrepancy = " << max_f_x<F, X>(f, x) << "\n\n";
	return;
}

template <typename F, typename X>
void print_disc(const F& f, const X& x, EMX W) {
	cout << "\n\nx[0] = " << x[0] << "      x[1] = " << x[1] << "        f_0(x) = " << f[0](x) << "        f_1(x) = " << f[1](x) << "        discrepancy = " << max_f_x<F, X>(f, x) <<std::endl <<
		"Jacobi Matrix: " << std::endl << W << "\n\n";
	return;
}

v_d Newtonw_s(const m_func_x< double(*)(EVX x)>& d_fi, 
			  const v_func_v_x<EVX>& f, 
			  const v_d& x, const double& eps) {

	EVX e_v_X(x.size());
	EVX e_v_X_1(x.size());
	e_v_X = conv_v_e(x);
	EMX W(x.size(), x.size());
	W = W_completion(d_fi, e_v_X);
	e_v_X_1 = e_v_X - W.inverse() * Fx(e_v_X, f);
	print_disc<v_func_v_x<EVX>, EVX>(f, e_v_X);
	print_disc<v_func_v_x<EVX>, EVX>(f, e_v_X_1,W);

	while (max_delt<EVX>(e_v_X, e_v_X_1) > eps) {
		cout << (max_delt<EVX>(e_v_X, e_v_X_1)) << std::endl;
		e_v_X = e_v_X_1;
		cout << "+1" ;
		for (unsigned int i = 0; i < f.size(); i++) {
			W = W_completion(d_fi, e_v_X);
			e_v_X_1 = e_v_X - W.inverse() * Fx(e_v_X, f);
		}
		print_disc<v_func_v_x<EVX>, EVX>(f, e_v_X_1, W);

	}
	cout << (max_delt<EVX>(e_v_X, e_v_X_1)) << "\n\n";

	return converter_e_v(e_v_X_1);
}