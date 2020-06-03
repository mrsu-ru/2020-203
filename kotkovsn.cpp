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
    
    delete []r;
    delete []Ar;
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

/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void kotkovsn::lab8()
{
  const double eps = 1e-18;

  double err = 0;
  for (int i = 0; i < N; i++) 
    for (int j = i + 1; j < N; j++)
      err += A[i][j] * A[i][j] + A[j][i] * A[j][i];
  
  double **C = new double*[N];
  for (int i = 0; i < N; i++)
    C[i] = new double[N];

  while (err > eps)
  {
    double alpha = 0;

    int i = -1, j = -1;
    double max = -1e9;
    for (int ii = 0; ii < N; ii++)
    {
      for (int jj = ii + 1; jj < N; jj++)
      {
        if (fabs(A[ii][jj]) > max)
        {
           max = fabs(A[ii][jj]);
           i = ii; j = jj; 
        }
      }
    }

    if (fabs(A[i][i] - A[j][j]) < eps)   
      alpha = atan(1);
    else
      alpha = atan(2 * A[i][j] / (A[j][j] - A[i][i])) / 2;

    double s = sin(alpha), c = cos(alpha);

    for (int ii = 0; ii < N; ii++)
      for (int jj = 0; jj < N; jj++)
          C[ii][jj] = A[ii][jj];     
      
    C[i][i] = c * c * A[i][i] - 2 * s * c * A[i][j] + s * s * A[j][j];
    C[j][j] = s * s * A[i][i] + 2 * s * c * A[i][j] + c * c * A[j][j];
    C[i][j] = C[j][i] = (c * c - s * s) * A[i][j] + s * c * (A[i][i] - A[j][j]);
    for (int k = 0; k < N; k++)
    {
      if (k == i || k == j) continue;
      C[i][k] = C[k][i] = c * A[i][k] - s * A[j][k];
      C[j][k] = C[k][j] = s * A[i][k] + c * A[j][k];
    }

    err = 0;
    for (int i = 0; i < N; i++) 
      for (int j = i + 1; j < N; j++)
        err += C[i][j] * C[i][j] + C[j][i] * C[j][i];

    for (int ii = 0; ii < N; ii++)
      for (int jj = 0; jj < N; jj++)
          A[ii][jj] = C[ii][jj];
  }

  for (int i = 0; i < N; i++)
  {
      cout  << "lambda_" << i << " = " << A[i][i] << endl;  
      delete []C[i];
  }
  
  delete []C;
}

/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void kotkovsn::lab9()
{
  const double eps = 1e-3;

  double *yPrev = new double[N];
  for (int i = 0; i < N; i++)
      yPrev[i] = 1;
  
  double *y = new double[N];
  
  double delta = 1e9;
  double lambdaMax = 0;

  while (delta > eps)
  {
    memset(y, 0, sizeof(double) * N);
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        y[i] += A[i][j] * yPrev[j];

    int i = 0;
    while (i < N && fabs(y[i]) < eps && fabs(yPrev[i]) < eps) i++;

    double lambda = y[i] / yPrev[i];
    delta = fabs(lambda - lambdaMax);
    lambdaMax = lambda;
    
    for (int i = 0; i < N; i++)
      yPrev[i] = y[i];
  }

  cout << "lambdaMax = " << lambdaMax << endl;

  delete []yPrev;
  delete []y;
}


std::string kotkovsn::get_name()
{
  return "Kotkov S.N.";
}
