#include "zhalninrv.h"
#include "parshinad.h"

/**
 * Введение в дисциплину
 */
void parshinad::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void parshinad::lab2()
{
	double eps = 1e-12;
	for (int i = 0; i < N; i++) {
		if (fabs(A[i][i]) < eps) {
			for (int l = i + 1; l < N; l++) {
				if (fabs(A[l][i]) > eps) {
					double* tmpS = A[i];
					A[i] = A[l];
					A[l] = tmpS;
					swap(b[i], b[l]);
					break;
				}
			}
		}


		if (fabs(A[i][i] - 1.0) > eps) {
			double tmp = A[i][i];
			for (int k = i + 1; k < N; k++) {
				A[i][k] /= tmp;
			}
			A[i][i] = 1;
			b[i] /= tmp;
		}

		for (int j = i + 1; j < N; j++) {
			double tmp2 = A[j][i];
			for (int k = i; k < N + 1; k++) {
				A[j][k] -= tmp2 * A[i][k];
			}
			b[j] -= tmp2 * b[i];
		}
	}

	for (int i = 1; i < N; i++) {
		for (int j = N - i - 1; j >= 0; j--) {
			double tmp = A[j][N - i];
			b[j] -= tmp * b[N - i];
		}
	}

	for (int i = 0; i < N; i++) {
		x[i] = b[i];
	}
}



/**
 * Метод прогонки
 */
void parshinad::lab3()
{

}



/**
 * Метод простых итераций
 */
void parshinad::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void parshinad::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void parshinad::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void parshinad::lab7()
{

}


void parshinad::lab8()
{

}


void parshinad::lab9()
{

}


std::string parshinad::get_name()
{
  return "Parshin A.D.";
}
