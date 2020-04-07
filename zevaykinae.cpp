#include "zevaykinae.h"

/**
 * Введение в дисциплину
 */
void zevaykinae::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void zevaykinae::lab2()
{

	for (int i = 0; i < N; i++) 
		x[i] = b[i]; 
	long double m; 
	for (int k = 0; k < N - 1; k++) { 
		for (int i = k + 1; i < N; i++) { 
			m = A[i][k] / A[k][k]; 	
			for (int j = k; j < N; j++) { 
				A[i][j] = A[i][j] - m * A[k][j]; 
			} 
			x[i] = x[i] - m * x[k]; 
		} 
	} 
	for (int i = N - 1; i >= 0; i--) { 
		for (int j = i + 1; j < N; j++) 
		x[i] = x[i] - A[i][j] * x[j]; 
		x[i] = x[i] / A[i][i]; 
	} 
}



/**
 * Метод прогонки
 */
void zevaykinae::lab3()
{
	double *alpha = new double [N];
	double *beta = new double [N];

	alpha[0] = -A[0][1]/A[0][0];
	beta[0] = b[0]/A[0][0];

	for(int i=1; i<N; i++) //здесь определяются прогоночные коэффициенты
	{
		alpha[i] = A[i][i+1]/(-A[i][i] - A[i][i-1]*alpha[i-1]);
		beta[i] = (-b[i] + A[i][i-1]*beta[i-1])/(-A[i][i] - A[i][i-1]*alpha[i-1]);
	}

	x[N-1] = beta[N-1];
	for(int i=N-2; i>=0; i--) //решение
		x[i] = alpha[i]*x[i+1] + beta[i];

	delete [] alpha;
	delete [] beta;
}



/**
 * Метод простых итераций
 */
void zevaykinae::lab4()
{
	double** S = new double* [N];
	for (int i = 0; i < N; i++) {
		S[i] = new double[N];
		for(int j = 0; j < N; j++)
			S[i][j] = 0;
	}
	int* D = new int[N];
	for (int i = 0; i < N; i++)
		D[i] = 0;
	double temp;
	for(int i = 0; i < N; i++){
		temp = A[i][i];
		for (int j = 0; j < i; j++)
			temp -= D[j] * S[j][i] * S[j][i];
		D[i] = (temp > 0)? 1: -1;
		S[i][i] = sqrt(D[i] * temp);
		double nakSum;
		for (int j = i + 1; j < N; j++) {
			nakSum = 0;
			for (int k = 0; k < j; k++) 
				nakSum += D[k] * S[k][i] * S[k][j];
			S[i][j] = (A[i][j] - nakSum) / (D[i] * S[i][i]);
		}
	}
	double* y = new double[N];
	for (int i = 0; i < N; i++) {
		b[i] /= S[i][i];
		y[i] = b[i];
		for (int j = i + 1; j < N; j++)
			b[j] -= b[i] * S[i][j];
	}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			S[i][j] *= D[i];
	for (int i = N - 1; i >= 0; i--) {
		y[i] /= S[i][i];
		x[i] = y[i];
		for(int j = i - 1; j >= 0; j--)
			y[j] -= y[i] * S[j][i];
	}
	
	for (int i = 0; i < N; i++)
		delete[]S[i];
	delete[]S;
	delete[]D;
	delete[]y;
}



/**
 * Метод Якоби
 */
void zevaykinae::lab5()
{
	double eps = 1e-14;
	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double *p_x = new double[N];
	double norma = 0;
	do {
		for (int i = 0; i < N; i++) {
			p_x[i] = x[i];
		}
		for (int i = 0; i < N; i++) {
			double sum = 0;
			for (int j = 0; j < N; j++) {
				if (i != j) {
					sum += (A[i][j] * p_x[j]);
				}
			}

			x[i] = (b[i] - sum) / A[i][i];
		}
		norma = 0;
		for (int i = 0; i < N; i++) {
			if ((p_x[i] - x[i]) > norma) {
				norma = abs(p_x[i] - x[i]);
			}
		}
	} while (norma > eps);

	delete[] p_x;
}



/**
 * Метод минимальных невязок
 */
void zevaykinae::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void zevaykinae::lab7()
{

}


void zevaykinae::lab8()
{

}


void zevaykinae::lab9()
{

}


std::string zevaykinae::get_name()
{
  return "Zevaykin A.E.";
}
