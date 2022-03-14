//This code is the intellectual property of Danila Musatov .
//Contacts danilarumus2000@gmail.com
//This program contains a number of numerical methods for finding the root on the interval.
// Copyright reserved 2021 

#include<iostream>
#include<vector>	
#include<cmath>
#include"functions.h"
#include"polynoms.h"

int main() {
	double eps = 0.00001;
	
	v_func_v_x<v_d> my_sys_f{ &my_polynom_1<v_d> ,&my_polynom_2<v_d> };
	v_func_v_x<v_d> my_sys_fi{ & my_fi_1<v_d>, & my_fi_2<v_d>};
	m_func_x<double(*)(EVX)> d_f{ {my_f_1_d_1<EVX>, my_f_1_d_2<EVX>}, {my_f_2_d_1<EVX>, my_f_2_d_2<EVX> }};
	//std::cin >> eps;
	v_d x_0{ 1, 0 };
	v_d x_1;
	v_func_v_x<EVX> my_sys_e_f{ &my_polynom_1<EVX> ,&my_polynom_2<EVX> };
	std::cout << "\n\n Newton\n\n";
	x_1 = Newtonw_s(d_f, my_sys_e_f, x_0, eps);
	auto iter = x_1.begin();  // получаем итератор
	while (iter != x_1.end())    // пока не дойдем до конца
	{
		std::cout << *iter << std::endl;// получаем элементы через итератор
		++iter;             // перемещаемся вперед на один элемент
	}
	system("pause");
	return 0;
}