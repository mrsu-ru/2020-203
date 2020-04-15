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
 * Метод простых итераций
 */
void simatovvv::lab4()
{
double** L = new double*[N];
	for (int i=0; i<N; i++)
	L[i] = new double[N];
	double* y = new double[N];

   for (int i=0; i<N; i++)
   {
		for (int j=0; j<N; j++)
			{
				L[i][j]=0;
			}
   }

  for(int i=0;i<N;i++)
	{
		for(int k=0;k<=i-1;k++)

			L[i][i]+=L[i][k]*L[i][k];
			L[i][i]=sqrt(A[i][i]-L[i][i]);
			for(int j=i+1;j<N;j++)
				{
					for(int k=0;k<=i-1;k++)
						L[j][i]+=L[i][k]*L[j][k];
						L[j][i]=(A[i][j]-L[j][i])/L[i][i];
				}
  }
	double summa=0;
    for(int i=0;i<N;i++)
		{
			for(int j=0;j<i;j++){
			summa+=L[i][j]*y[j];}
			y[i]=(b[i]-summa)/L[i][i];
			summa=0;
		}
	for(int i=N-1;i>=0;i--)
		{
			for(int j=i+1;j<N;j++){
				summa+=L[j][i]*x[j];}
				x[i]=(y[i]-summa)/L[i][i];
				summa=0;
		}
delete[] y;
}



/**
 * Метод Якоби или Зейделя
 */
void simatovvv::lab5()
{
long double eps = 0.00000001;
    long double* p = new long double[N];
	long double norm;
    for (int i = 0; i < N; i++)
        x[i] = 0;
    do {
		for (int i = 0; i < N; i++)
        {
			p[i] = b[i];
			for (int j = 0; j < N; j++)
				{
				    if (i != j)
                     p[i] -= A[i][j] * x[j];
				}
			p[i] /= A[i][i];
		}
        norm = fabs(x[0] - p[0]);
		for (int h = 0; h < N; h++)
        {
			if (fabs(x[h] - p[h]) > norm)
				norm = fabs(x[h] - p[h]);
			x[h] = p[h];
		}
	} while (norm > eps);
    delete[] p;
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
