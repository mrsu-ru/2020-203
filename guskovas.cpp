#include "guskovas.h"

/**
 * Введение в дисциплину
 */
void guskovas::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void guskovas::lab2()//int n, int m, double e, double** arr, double* x) {
{

	for(int i = 0; i < N; i++){
		int indexMax = i; 
		for (int j = i + 1; j < N; j++) { 
			if (abs(A[j][i]) > abs(A[indexMax][i])) indexMax = j; 
		} 
		if (indexMax != i) { 
			for (int k = 0; k < N; k++) { 
				swap(A[i][k], A[indexMax][k]); 
			} 
		} 
	}

	//ТУДА
	int n = N;

	for(int k = 0; k < n-1; k++){
		for(int i = k+1; i < n; i++){
			double c = A[i][k] / A[k][k];
			//A[i][k] = 0;
			for(int j = k+1; j < n; j++){
				A[i][j] -= A[k][j] * c; 
			}
			b[i] -= b[k] * c;
		}
	}

	for (int i = 0; i < n; i++){
			x[i] = b[i];
		}

	//Обратно

	for(int i = n-1; i >= 0; i--){
		double S = 0;
		for(int j = i+1; j < n; j++){
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
	}
	
}



/**
 * Метод прогонки
 */
void guskovas::lab3()//N, A, b, x
{
	double *ALFA = new double[N];
	double *BETA = new double[N];

	// ТУДА
	//	находим ALFA и BETA
	ALFA[0] = -A[0][1] / A[0][0];
	BETA[0] = b[0] / A[0][0];
	for (int i = 1; i < N; i++) {
		ALFA[i] = -A[i][i + 1] 
					/ 
					(
						A[i][i] + 
						A[i][i - 1] * ALFA[i - 1]
					);
		BETA[i] =   (b[i] - A[i][i - 1] * BETA[i - 1]) 
					/ 
					(A[i][i] + A[i][i - 1] * ALFA[i - 1]);
	}

	// Обратная прогонка
	//после общего вида формул подставляем иксы 
	x[N - 1] = BETA[N - 1];
	for (int i = N - 2; i >= 0; i--) {
		x[i] = ALFA[i] * x[i + 1] + BETA[i];
	}
}



/**
 * Метод Холецкого
 */
void guskovas::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void guskovas::lab5()//якоби
{
	double *f = new double[N];

	double norma = 0;//error
	do {
		for (int i = 0; i < N; i++) f[i] = x[i];

		for (int i = 0; i < N; i++) {
			double result = b[i];

			for (int j = 0; j < N; j++) result = i != j ? 
					result - (A[i][j] * f[j]) 
					:
					result
			;

			x[i] = result / A[i][i];
		}

		norma = 0;
		for (int i = 0; i < N; i++) norma = abs(f[i] - x[i]) > norma ? 
			abs(f[i] - x[i]) 
			: 
			norma
		;

	} while (norma > 1e-20);

}



/**
 * Метод минимальных невязок
 */
void guskovas::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void guskovas::lab7()
{

}


void guskovas::lab8()
{

}


void guskovas::lab9()
{

}


std::string guskovas::get_name()
{
  return "Guskov A.S.";
}
