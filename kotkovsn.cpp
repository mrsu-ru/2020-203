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

  delete []alpha;
  delete []beta;
}



/**
 * Метод Холецкого
 */
void kotkovsn::lab4()
{
  double **S = new double*[N];
  double **D = new double*[N];
  for (int i = 0; i < N; i++)
  {
    S[i] = new double[N];
    memset(S[i], 0, sizeof(double) * N);
    D[i] = new double[N];
    memset(D[i], 0, sizeof(double) * N);
  }
 
  for (int i = 0; i < N; i++)
  {
    double tmp = A[i][i];
    for (int k = 0; k <= i - 1; k++)
      tmp -= D[k][k] * S[k][i] * S[k][i];

    if (tmp > 0)   D[i][i] = 1;
    else           D[i][i] = -1;

    S[i][i] = sqrt(D[i][i] * tmp);

    for (int j = i + 1; j < N; j++)
    {
      double sum = 0;
      for (int k = 0; k <= j - 1; k++)
        sum += D[k][k] * S[k][i] * S[k][j];
      S[i][j] = (A[i][j] - sum) / (D[i][i] * S[i][i]);
    }
  }
  
  double *y = new double[N];
  for (int i = 0; i < N; i++)
  {
    b[i] /= S[i][i];
    y[i] = b[i];
    for (int j = i + 1; j < N; j++)
      b[j] -= b[i] * S[i][j];
  }

  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      S[i][j] *= D[i][i];

  for (int i = N - 1; i >= 0; i--)
  {
    y[i] /= S[i][i];
    x[i] = y[i];
    for (int j = i - 1; j >= 0; j--)
      y[j] -= y[i] * S[j][i];
  }

  delete []y;

  for (int i = 0; i < N; i++)
  {
    delete []S[i];
    delete []D[i];
  }

  delete []S;
  delete []D;
}



/**
 * Метод Якоби или Зейделя
 */
void kotkovsn::lab5()
{
  const double eps = 1e-18;
  double *z = new double[N];

  for (int k = 0; k < N; k++)
  {
    double y = 0;
    for (int i = 0; i < N; i++)
      y += A[k][i] * x[i];
    z[k] = b[k] - y;
  }

  double zNorm = 0;

  for (int k = 0; k < N; k++)
    zNorm += z[k] * z[k];
    
  while (zNorm > eps * eps)
  {
    for (int k = 0; k < N; k++)
    {
      double sum1 = 0;
      double sum2 = 0;
      for (int j = 0; j < k; j++)
        sum1 += A[k][j] * x[j];

      for (int j = k + 1; j < N; j++)
        sum2 += A[k][j] * x[j];
      
      x[k] = (b[k] - sum1 - sum2) / A[k][k];
    }

    for (int k = 0; k < N; k++)
    {
      double y = 0;
      for (int i = 0; i < N; i++)
        y += A[k][i] * x[i];
      z[k] = b[k] - y;
    }

    zNorm = 0;

    for (int k = 0; k < N; k++)
      zNorm += z[k] * z[k];
  }
  
  delete []z;
}



/**
 * Метод минимальных невязок
 */
void kotkovsn::lab6()
{
    const double eps = 1e-18;

    double *r = new double[N];    
    for (int i = 0; i < N; i++)
      r[i] = -b[i];
    for (int i = 0; i < N; i++) 
      for (int j = 0; j < N; j++)
            r[i] += A[i][j] * x[j];

    double rNorm = 0;
    for (int i = 0; i < N; i++)
      rNorm += r[i] * r[i];

    double *Ar = new double[N];
    memset(Ar, 0, sizeof(double) * N);
    for (int i = 0; i < N; i++) 
      for (int j = 0; j < N; j++)
            Ar[i] += A[i][j] * r[j];

    double ArNorm = 0;
    for (int i = 0; i < N; i++)
      ArNorm += Ar[i] * Ar[i];

    double tau = 0;
    for (int i = 0; i < N; i++)
      tau += Ar[i] * r[i];
    tau /= ArNorm;

    while (rNorm > eps * eps)
    {
      for (int i = 0; i < N; i++)
        x[i] = x[i] - tau * r[i];

      for (int i = 0; i < N; i++)
        r[i] = -b[i];
      for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++)
              r[i] += A[i][j] * x[j];      

      memset(Ar, 0, sizeof(double) * N);
      for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++)
              Ar[i] += A[i][j] * r[j];

      ArNorm = 0;
      for (int i = 0; i < N; i++)
        ArNorm += Ar[i] * Ar[i];

      tau = 0;
      for (int i = 0; i < N; i++)
        tau += Ar[i] * r[i];
      tau /= ArNorm;

      rNorm = 0;
      for (int i = 0; i < N; i++)
        rNorm += r[i] * r[i];
    }
}



/**
 * Метод сопряженных градиентов
 */
void kotkovsn::lab7()
{
  const double eps = 1e-18;

  double *xPrev = new double[N];
  memset(xPrev, 0, sizeof(double) * N);

  double *r = new double[N];
  for (int i = 0; i < N; i++)
    r[i] = -b[i];

  double *Ar = new double[N];
  memset(Ar, 0, sizeof(double) * N);
    for (int i = 0; i < N; i++) 
      for (int j = 0; j < N; j++)
            Ar[i] += A[i][j] * r[j];

  double rrScalar = 0;
  for (int i = 0; i < N; i++)
    rrScalar += r[i] * r[i];

  double rArScalar = 0;
  for (int i = 0; i < N; i++)
    rArScalar += Ar[i] * r[i];

  double alpha = 1;
  double tau = rrScalar / rArScalar;

  for (int i = 0; i < N; i++)
    x[i] = tau * b[i];
    
  for (int i = 0; i < N; i++) 
    for (int j = 0; j < N; j++)
          r[i] += A[i][j] * x[j];

  rrScalar = 0;
  for (int i = 0; i < N; i++)
    rrScalar += r[i] * r[i];

  while (rrScalar > eps * eps)
  {
    
    double rArScalarPrev = rArScalar;

    memset(Ar, 0, sizeof(double) * N);
    for (int i = 0; i < N; i++) 
      for (int j = 0; j < N; j++)
            Ar[i] += A[i][j] * r[j];

    rArScalar = 0;
    for (int i = 0; i < N; i++)
      rArScalar += Ar[i] * r[i];
    
    double tauPrev = tau;

    tau = rrScalar / rArScalar;
    alpha = 1.0 / (1 - tau / tauPrev / alpha * rrScalar / rArScalarPrev);

    for (int i = 0; i < N; i++)
    {
      double tmp = x[i];
      x[i] = alpha * x[i] + (1 - alpha) * xPrev[i] - tau * alpha * r[i];
      xPrev[i] = tmp;
    }

    for (int i = 0; i < N; i++)
      r[i] = -b[i];
    for (int i = 0; i < N; i++) 
      for (int j = 0; j < N; j++)
            r[i] += A[i][j] * x[j];

    rrScalar = 0;
    for (int i = 0; i < N; i++)
      rrScalar += r[i] * r[i];
  }

  delete []xPrev;
  delete []r;
  delete []Ar;
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
