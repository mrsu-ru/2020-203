#include "kotkovsn.h"

/**
 * Введение в дисциплину
 */
void kotkovsn::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kotkovsn::lab2()
{
    const double eps = 1e-12;
    for (int i = 0; i < N; i++)
    {
      if (fabs(A[i][i]) < eps)
        for (int j = i + 1; j < N; j++)
        {
          if (fabs(A[j][i]) > eps)
          {
            double *tmp = A[i];
            A[i] = A[j];
            A[j] = tmp;

            double temp = b[i];
            b[i] = b[j];
            b[j] = temp;
            break;
          }
        }

      if (fabs(A[i][i] - 1) > eps)
      {
        for (int j = i + 1; j < N; j++)
          A[i][j] /= A[i][i];
        b[i] /= A[i][i];
        A[i][i] = 1;
      }

      for (int j = 0; j < i; j++)
      {
        for (int k = i + 1; k < N; k++)
          A[j][k] -= A[i][k] * A[j][i];
        b[j] -= b[i] * A[j][i];
        A[j][i] = 0; 
      }

     for (int j = i + 1; j < N; j++)
      {
        for (int k = i + 1; k < N; k++)
          A[j][k] -= A[i][k] * A[j][i];
        b[j] -= b[i] * A[j][i];
        A[j][i] = 0; 
      }
    }

    for (int i = 0; i < N; i++)
      x[i] = b[i];
}


 
/**
 * Метод прогонки
 */
void kotkovsn::lab3()
{
  double *alpha = new double[N];
  double *beta = new double[N];

  alpha[0] = A[0][1] / A[0][0];
  beta[0] = b[0] / A[0][0];

  for (int i = 1; i < N; i++)
  {
    alpha[i] = A[i][i + 1] / (A[i][i] - A[i][i - 1] * alpha[i - 1]);
    beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / (A[i][i] - A[i][i - 1] * alpha[i - 1]);
  }

  x[N - 1] = beta[N - 1];

  for (int i = N - 2; i >= 0; i--)
    x[i] = beta[i] - alpha[i] * x[i + 1];
}



/**
 * Метод простых итераций
 */
void kotkovsn::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void kotkovsn::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void kotkovsn::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void kotkovsn::lab7()
{

}


void kotkovsn::lab8()
{

}


void kotkovsn::lab9()
{

}


std::string kotkovsn::get_name()
{
  return "Kotkov S.N.";
}
