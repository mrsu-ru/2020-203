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
		double* vrem;
		vrem = A[index];
		A[index] = A[i];
		A[i] = vrem;

		for (int j = N; j > i; A[i][j--] /= A[i][i]);
			A[i][i] = 1;
		for (int j = i + 1; j < N; j++)
		{
			for (int k = N; k > i; k--)
				A[j][k] -= A[i][k] * A[j][i];
				A[j][i] = 0;
		}

	}
	//обратный ход
	float ved1 = 0, ved2 = 0;
	for (int i = N - 1; i > 0; i--)
	{
		for (int k = i; k > 0; k--)
		{
			ved2 = A[k-1][i];
			for (int j = k; j < N + 1; j++)
			{
				ved1 = -A[i][j];
				A[k-1][j] += ved1 * ved2;
			}
		}
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
