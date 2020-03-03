//#include "zhalninrv.h"
#include "edelevaup.h"

/**
 * Введение в дисциплину
 */
void edelevaup::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void edelevaup::lab2()
{

	for (int i = 0; i < N; i++) {

		int max = i;
		for (int j = i + 1; j < N; j++) {
			if (abs(A[j][i]) > abs(A[max][i]))
				max = j;
		}
		double *s;
		double k;
		if (max != i) {
			s = A[i];
			A[i] = A[max];
			A[max] = s;
			k = b[i];
			b[i] = b[max];
			b[max] = k;
		}
		
		b[i] /= A[i][i];
		for (int j = N - 1; j > i; j--)
			A[i][j] /= A[i][i];
		A[i][i] = 1;

		for (int j = i + 1; j < N; j++) {

			for (int k = N; k > i; k--) {
				A[j][k] -= A[i][k] * A[j][i];
			}
			b[j] -= b[i] * A[j][i];

			A[j][i] = 0;
		}


	}

	for (int i = N - 1; i > 0; i--) {
		double s = 0, l = 0;
		for (int j = i; j < N; j++) {
			s += b[j ] * A[i - 1][j];
			A[i - 1][j] = 0;
			
		}
		b[i - 1] -= s;
		
	}

	for (int i = 0; i < N; i++) {
		x[i] = b[i];
		
		}
	

}



/**
 * Метод прогонки
 */
void edelevaup::lab3()
{

}



/**
 * Метод простых итераций
 */
void edelevaup::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void edelevaup::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void edelevaup::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void edelevaup::lab7()
{

}


void edelevaup::lab8()
{

}


void edelevaup::lab9()
{

}


std::string edelevaup::get_name()
{
  return "Edeleva U.P.";
}
