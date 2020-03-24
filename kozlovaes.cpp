#include "kozlovaes.h"

/**
 * Введение в дисциплину
 */
void kozlovaes::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kozlovaes::lab2()
{
	for (int i = 0; i < N; i++) {
		int max = i;
		for (int j = i + 1; j < N; j++) {
			if (abs(A[j][i]) > abs(A[max][i]))
				max = j;
		}

		double *s;
		double tmp;
		if (max != i) {
			s = A[i]; A[i] = A[max]; A[max] = s;
			tmp = b[i]; b[i] = b[max]; b[max] = tmp;
		}

		for (int j = N; j > i; A[i][j--] /= A[i][i]);
		b[i] /= A[i][i];
		A[i][i] = 1;

		for (int j = i + 1; j < N; j++) {

			for (int k = N; k > i; k--) A[j][k] -= A[i][k] * A[j][i];
			b[j] -= b[i] * A[j][i];
			A[j][i] = 0;
		}
	}
	
	for (int i = N - 1; i > 0; i--) {
		double s = 0;
		for (int j = i ; j < N; j++) {
			s += b[j] * A[i - 1][j];
			A[i - 1][j] = 0;
		}
		b[i - 1] -=s;
	}
	for (int i = 0; i < N; i++) {
		x[i] = b[i];
	}
}



/**
 * Метод прогонки
 */
void kozlovaes::lab3()
{
	double *alfa = new double[N];
	double *betta = new double[N];
	int i;
	alfa[0] = A[0][1]/-A[0][0];
	betta[0] = -b[0]/-A[0][0];
 

	for(i = 1;i < N;i++){
		alfa[i] = A[i][i+1]/(-A[i][i]-alfa[i-1]*A[i][i-1]);
		betta[i] = (-b[i]+A[i][i-1]*betta[i-1])/(-A[i][i]-alfa[i-1]*A[i][i-1]);
	}
	i=N-1;
	x[N-1] = (-b[i]+A[i][i-1]*betta[i-1])/(-A[i][i]-alfa[i-1]*A[i][i-1]);

	for(int i=N-1;i>=0;i--){
		x[i] = alfa[i]*x[i+1]+betta[i];
	}
}



/**
 * Метод Холецкого
 */
void kozlovaes::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void kozlovaes::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void kozlovaes::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void kozlovaes::lab7()
{

}


void kozlovaes::lab8()
{

}


void kozlovaes::lab9()
{

}


std::string kozlovaes::get_name()
{
  return "Kozlova E.S.";
}
