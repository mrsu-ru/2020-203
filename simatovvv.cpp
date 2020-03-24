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
double *alpha = new double[N];
double *beta = new double[N];
alpha[0] = -A[0][1] / A[0][0];
beta[0] = b[0] / A[0][0];
for (int i = 1; i <= N-2; i++)
{
alpha[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * alpha[i - 1]);
beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / (A[i][i] + A[i][i - 1] * alpha[i - 1]);

}
beta[N - 1] = (b[N - 1] - A[N - 1][N - 2] * beta[N - 2]) / (A[N - 1][N - 1] + A[N - 1][N - 2] * alpha[N - 2]);
x[N - 1] = beta[N - 1];
for (int i = N-2; i >= 0; i--)
{
x[i] = alpha[i] * x[i + 1] + beta[i];
}
}



/**
 * Метод Холецкого
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
