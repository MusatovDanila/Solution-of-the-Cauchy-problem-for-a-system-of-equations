import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin, cos, fabs, log
from tabulate import tabulate
from texttable import Texttable
import latextable

# Функция сохранение в таблицы tex
def savetable(array,numberOfTable):
    table = Texttable()

    table.set_cols_align(["c"] * len(array[0]))
    table.set_cols_dtype(['t'] * len(array[0]))
    table.set_deco(Texttable.HEADER | Texttable.VLINES | Texttable.HLINES)
    table.add_rows(array)

    path = "C:/Users/Danila/Documents/Study/7 semestor/Numerical methods(grid models of partial differential equations)/Courses work/Report/table_" + str(numberOfTable) + ".tex"
    my_file = open(path, 'w+')
    my_file.write(latextable.draw_latex(table))
    my_file.close()
    return

# Метод эйлера для решения системы исходя из заданных условий
def Method_Euler_sys(func1, func2, y1Start, y2Start, a, b, h):
    xs = np.arange(a, b+h, h)
    ys1, ys2 = [y1Start], [y2Start] # первые значения приближений
    for point in enumerate(xs[:-1]):
        ys1.append(ys1[-1] + h*func1(ys2[-1]))
        ys2.append(ys2[-1] + h*func2(ys1[-2]))
    return [xs, ys1, ys2]

# Правило рунге для оценки погрешности
def Runge_rule(ys1, ys2, p):
    return   max([np.abs(y2Now - y1) for y1, y2Now in zip(ys1, ys2[::2])]) / (2**p - 1)

# Вывод значений на каждой итерации и запись в таблицу tex
def print_table(xs, ys1, ys2, prec=6):
    Тable5 = []
    for i, (x,y1,y2Now) in enumerate(zip(xs, ys1, ys2)):
        Тable5.append([i, round(xs[i], prec), round(ys1[i], prec), round(ys2[i], prec)])
        print(f"Итерация {i} : x={round(x, prec)}  y1={round(y1, prec)}   y2Now={round(y2Now, prec)}")
    return Тable5
    
dy_2_func = lambda x, y, z: cos(y) # заданная функция с новой переменной
dy_1_func = lambda x, y, z: sin(z) # введение новой переменной

# определение коэффициентов K и L
K1 = lambda x, y, z, h: h*dy_1_func(x, y, z)
L1 = lambda x, y, z, h: h*dy_2_func(x, y, z)

K2 = lambda x, y, z, h: h*dy_1_func(x + h/2, y + K1(x, y, z, h)/2, z + L1(x, y, z, h)/2)
L2 = lambda x, y, z, h: h*dy_2_func(x + h/2, y + K1(x, y, z, h)/2, z + L1(x, y, z, h)/2)

K3 = lambda x, y, z, h: h*dy_1_func(x + h/2, y + K2(x, y, z, h)/2, z + L2(x, y, z, h)/2)
L3 = lambda x, y, z, h: h*dy_2_func(x + h/2, y + K2(x, y, z, h)/2, z + L2(x, y, z, h)/2)

K4 = lambda x, y, z, h: h*dy_1_func(x + h, y + K3(x, y, z, h), z + L3(x, y, z, h))
L4 = lambda x, y, z, h: h*dy_2_func(x + h, y + K3(x, y, z, h), z + L3(x, y, z, h))


def Runge_Kutta_sys(a=1, b=3, h=0.01, y_1 = 0.980444, y_2 =  0.196798):
    x = np.arange(a, b + h, h)
    y1Now = y_1 # начальное значение y1
    y2Now = y_2 # начальное значение y2
    y1 = []
    y2 = []
    y1.append(y1Now)
    y2.append(y2Now)
    mx = 0
    Table3 = []
    j = 0

    # Рассчет согласно формулам
    for index, i in enumerate(x):
        #print("!")
        k1 = K1(i, y1Now, y2Now, h)
        k2 = K2(i, y1Now, y2Now, h)
        k3 = K3(i, y1Now, y2Now, h)
        k4 = K4(i, y1Now, y2Now, h)
        l1 = L1(i, y1Now, y2Now, h)
        l2 = L2(i, y1Now, y2Now, h)
        l3 = L3(i, y1Now, y2Now, h)
        l4 = L4(i, y1Now, y2Now, h)
        #print(k1, k2, k3, k4)
        #print(l1, l2, l3, l4)
        t = y1Now + (1/6) * (k1 + 2*k2 + 2*k3 + k4)
        y2Now = y2Now + (1/6) * (l1 + 2*l2 + 2*l3 + l4)
        y1Now = t
        Table3.append([round(y1[-1],5), round(k1, 5), round(k2, 5), round(k3, 5), round(k4, 5), round(y2[-1], 5), round(l1, 5), round(l2, 5), round(l3, 5), round(l4, 5)])
        if(i != x[len(x)-1]): 
            y1.append(y1Now)
            y2.append(y2Now)
            
    # Сохранение результатов в таблицу
    savetable(Table3,31)
    
    # Вывод программы и график
    
    plt.rcParams["figure.figsize"] = [8.0, 8.0]
    plt.rcParams["figure.autolayout"] = True    
    plt.grid(True)
    #print(x)
    #print(y)
    #print(z)
    plt.plot(x, y2, label="Решение задачи Коши, h={}".format(h))
    plt.plot(x, y1, label="Решение задачи Коши, h={}".format(h))
    plt.legend()
    plt.savefig('plot.png')
    plt.show()
    return x, y1, y2

def Euler_start(a=1, b=3, h=0.01, y_1 = 0.980444, y_2 =  0.196798):    
    y1_start_e = y_1
    y2_start_e = y_2
    func1 = lambda y2 : np.sin(y2) # первая функция
    func2 = lambda y1 : np.cos(y1) # вторая функция
    h1 = 0.01
    h2 = h1/2
    xs_h1, ys1_h1, ys2_h1 = Method_Euler_sys(func1, func2, y1_start_e, y2_start_e, a, b, h1)
    xs_h2, ys1_h2, ys2_h2 = Method_Euler_sys(func1, func2, y1_start_e, y2_start_e, a, b, h2)
    
    plt.plot(xs_h1, ys1_h1, label=f"y1(x), h={h1}")
    plt.plot(xs_h1, ys2_h1, label=f"y2(x), h={h1}")
    plt.grid(True)
    plt.legend()
    plt.savefig('plot.png')
    plt.show()

    plt.plot(xs_h2, ys1_h2, label=f"y1(x), h={h2}")
    plt.plot(xs_h2, ys2_h2, label=f"y2(x), h={h2}")
    plt.grid(True)
    plt.legend()
    plt.show()


    print("Погрешность для y1 : ", Runge_rule(ys1_h1, ys1_h2, 1))
    print("Погрешность для y2 : ", Runge_rule(ys2_h1, ys2_h2, 1))

    Table_5 = print_table(xs_h1, ys1_h1, ys2_h1)
    savetable(Table_5,51)
    print_table(xs_h2, ys1_h2, ys2_h2)
    return

def Runge_Kutta_start(a=1, b=3, h=0.01, y_1 = 0.980444, y_2 =  0.196798):    
    #решение СДУ с первыми краевыми условиями
    x_h, y_1_rq_h, y_2_rq_h = Runge_Kutta_sys(a, b, h, y_1, y_2)
    x_h2, y_1_rq_h2, y_2_rq_h2 = Runge_Kutta_sys(a, b, h/2, y_1, y_2)
    
    print("Погрешность для y1 : ", Runge_rule(y_1_rq_h, y_1_rq_h2, 4))
    print("Погрешность для y2 : ", Runge_rule(y_2_rq_h, y_2_rq_h2, 4))
    
    return
    

def main():
    a = 1
    b = 3
    h_min = 0.1
    #первое решение НСЛАУ
    y_1_1=0.980444
    y_2_1=0.196798

    #второе решение НСЛАУ
    y_1_2=0.226199
    y_2_2=-0.974081

    #решение СДУ с первыми краевыми условиями методом Рунге-Кутты
    Runge_Kutta_start(a, b, h_min, y_1_1, y_2_1)

    #решение СДУ со вторыми краевыми условиями методом Рунге-Кутты
    Runge_Kutta_start(a, b, h_min, y_1_2, y_2_2)

    #решение СДУ с первыми краевыми условиями методом Эйлера
    Euler_start(a, b, h_min, y_1_1, y_2_1)

    #решение СДУ со вторыми краевыми условиями методом Эйлера
    Euler_start(a, b, h_min, y_1_2, y_2_2)  
    
    return 1

main()