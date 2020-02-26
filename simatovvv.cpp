#include "simatovvv.h"

/**
 * Введение в дисциплину
 */
void simatovvv::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void simatovvv::lab2()
{
{
for (int i = 0; i < N; i++)
x[i] = b[i];
long double m;
for (int k = 0; k < N - 1; k++)
{
for (int i = k + 1; i < N; i++)
{
m = A[i][k] / A[k][k];
for (int j = k; j < N; j++)
{
A[i][j] = A[i][j] - m * A[k][j];
}
x[i] = x[i] - m * x[k];
}
}
for (int i = N - 1; i >= 0; i--)
{
for (int j = i + 1; j < N; j++)
x[i] = x[i] - A[i][j] * x[j];
x[i] = x[i] / A[i][i];
}
}
}



/**
 * Метод прогонки
 */
void simatovvv::lab3()
{

}



/**
 * Метод простых итераций
 */
void simatovvv::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void simatovvv::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void simatovvv::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void simatovvv::lab7()
{

}


void simatovvv::lab8()
{

}


void simatovvv::lab9()
{

}


std::string simatovvv::get_name()
{
  return "simatovvv";
}
