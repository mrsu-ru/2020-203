#include "manindi.h"

/**
 * Введение в дисциплину
 */
void manindi::lab1()
{
  
cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void manindi::lab2()
{
	int *sw1 = new int[N];
	for (int i = 0; i < N; i++)
	{
		sw1[i] = i;
	}
	for (int i = 0; i < N; i++)
	{
		int Maxi = i;
		for (int j = i + 1; j < N; j++)
		{
			if (abs(A[j][i]) > abs(A[Maxi][i]))
			{
				Maxi = j;
			}
		}
		if (Maxi != i)
		{
			swap(A[Maxi], A[i]);
			swap(b[Maxi], b[i]);
			swap(sw1[Maxi], sw1[i]);
		}
		b[i] /= A[i][i];
		for (int j = N - 1; j >= i; j--)
		{
			A[i][j] /= A[i][i];
		}

		for (int j = 0; j < N; j++)
		{
			if (j != i)
			{
				b[j] -= A[j][i] * b[i];
				for (int k = N - 1; k >= i; k--)
				{
					A[j][k] -= A[j][i] * A[i][k];
				}

			}

		}
	}
	for (int i = 0; i < N; i++) {
		x[i] = b[i];
	}
	for (int i = 0; i < N; i++) {
		if (sw1[i] != i) {
			swap(A[sw1[i]], A[i]);
			swap(b[sw1[i]], b[i]);
		}
	}
}



/**
 * Метод прогонки
 */
void manindi::lab3()
{
double *p = new double[N];

	double *q = new double[N];

        // Прямой ход

	p[0] = -A[0][1] / A[0][0];

	q[0] = b[0] / A[0][0];

	for (int i = 1; i < N; i++) {

		p[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * p[i - 1]);

		q[i] = (b[i] - A[i][i - 1] * q[i - 1]) / (A[i][i] + A[i][i - 1] * p[i - 1]);

	}

        // Обратный ход

	x[N - 1] = q[N - 1];

	for (int i = N - 2; i >= 0; i--) {

		x[i] = p[i] * x[i + 1] + q[i];

	}
        
}



/**
 * Метод Холецкого
 */
void manindi::lab4()
{
        double eps = 1e-20;
        double tau = 1e-5;
        double* pX = new double[N];

       while (true) {
      for (int i = 0; i < N; i++)

	pX[i] = x[i];
         for (int i = 0; i < N; i++) {
                  double sum = 0;
                    for (int j = 0; j < N; j++) 
                   sum += A[i][j] * pX[j];
                   x[i] = pX[i] - tau * (sum - b[i]);
                          }
                 double maxEr = abs(x[0] - pX[0]);
                  for (int i = 1; i < N; i++)
                    if (abs(x[i] - pX[i]) > maxEr)

				maxEr = abs(x[i] - pX[i]);

                     if (maxEr < eps)
                         break;
                                }
               
}



/**
 * Метод Якоби или Зейделя
 */
void manindi::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void manindi::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void manindi::lab7()
{

}


void manindi::lab8()
{

}


void manindi::lab9()
{

}


std::string manindi::get_name()
{
  return "Manin D.I.";
}
