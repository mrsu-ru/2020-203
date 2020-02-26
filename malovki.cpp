#include "malovki.h"

using namespace std;

/**
 * Введение в дисциплину
 */
void malovki::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void malovki::lab2()
{
	float max = 0;
	int index;
	//прямой ход
	for (int i = 0; i < N; i++)
	{
		index = i;
		max = abs(A[i][i]);
		for (int k = i; k < N; k++)
		{
			if (max < abs(A[k][i]))
			{
				max = abs(A[k][i]);
				index = k;
			}
		}
		if (i != index) {
			double* vrem;
			vrem = A[index];
			A[index] = A[i];
			A[i] = vrem;
			double temp = b[index];
			b[index] = b[i];
			b[i] = temp;
		}
		for (j = N-1; j > i; j--)
		{
			A[i][j] /= A[i][i];
		}
		b[i]/=A[i][i];
		A[i][i] = 1;
		
		for (j = i + 1; j < N; j++)
		{
			for (k = N-1; k > i; k--)
			{
				A[j][k] -= A[j][i] * A[i][k];
			}
			b[j]-=A[j][i]*b[i];
			A[j][i] = 0;
		}
	}

	//обратный ход
	x[N-1]=b[N-1];
	for (i = N - 2; i > -1; i--)
	{
		for (j = i + 1; j < N; j++)
		{
			sum += x[j] * A[i][j];
		}
		x[i] = b[i] - sum;
		sum = 0;
	}
}

/**
 * Метод прогонки
 */
void malovki::lab3()
{

}



/**
 * Метод простых итераций
 */
void malovki::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void malovki::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void malovki::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void malovki::lab7()
{

}


void malovki::lab8()
{

}


void malovki::lab9()
{

}


std::string malovki::get_name()
{
  return "Malov K. I.";
}
