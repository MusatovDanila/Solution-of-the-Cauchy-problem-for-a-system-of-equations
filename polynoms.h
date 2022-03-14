//This code is the intellectual property of Danila Musatov .
//Contacts danilarumus2000@gmail.com
//This program contains a number of numerical methods for finding the root on the interval.
// Copyright reserved 2021 

#pragma once
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
using v_d = vector<double>;
using EVX = Eigen::VectorXf;
using EMX = Eigen::MatrixXf;


template <typename X>
double my_polynom_1(X x) {
	return x[0] * x[0] + x[1] * x[1] - 1;
}


template <typename X>
double my_polynom_2(X x) {
	return sin(x[0] - x[1]) + 0.3 * x[0] - 1;
}

template <typename X>
double my_f_1_d_1(X x) {
	return 2 * x[0];
}

template <typename X>
double my_f_1_d_2(X x) {
	return 2*x[1];
}

template <typename X>
double my_f_2_d_1(X x) {
	return cos(x[0] - x[1]) + 0.3;
}

template <typename X>
double my_f_2_d_2(X x) {
	return -cos(x[0] - x[1]);
}

template <typename X>
double my_fi_1(X x) {
	return sqrt(1 - x[1] * x[1]);
}

template <typename X>
double my_fi_2(X x) {
return x[0] - asin(1 - 0.3 * x[0]);
}
