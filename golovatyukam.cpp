#include "golovatyukam.h"

/**
 * Введение в дисциплину
 */
void golovatyukam::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void golovatyukam::lab2()
{
	int m = N + 1;
	double eps =0.000001;
	
	double** arr = new double*[N];
	for (int i = 0; i < N; i++) {
		arr[i] = new double[N + 1];
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			arr[i][j] = A[i][j];
		}
	}

	for (int i = 0; i < N; i++) {
		arr[i][N] = b[i];
	}

	//cout << "\nПрямой ход: " << endl;
	for (int i = 0; i < N; i++) {

		//////////////////////////////////////////////////////////////
		//проаверяем наличие решения
		int max = i;
		for (int j = i + 1; j < N; j++) {
			if (abs(arr[j][i]) > abs(arr[max][i]))
				max = j;
		}
		if (max != i) {
			for (int k = 0; k < N + 1; k++) {
				swap(arr[i][k], arr[max][k]);
			}
		}

		if (abs(arr[i][i]) <= eps) {
			cout << "Данная система не имеет решений " << endl;
			exit(0);
		}
		/////////////////////////////////////////////////////////////

		double ved = arr[i][i];
		for (int j = i; j < N + 1; j++) {
			arr[i][j] /= ved;
		}
		for (int j = i + 1; j < N; j++) {
			ved = arr[j][i];
			arr[j][i] = 0;
			for (int k = i + 1; k < N + 1; k++)
				arr[j][k] -= arr[i][k] * ved;
		}
	}

	//cout << "\n\nОбратный ход: " << endl;

	for (int i = N - 1; i > 0; i--) {
		double s = 0;
		for (int j = i + 1; j < m; j++) {
			double ved = arr[i - 1][j - 1];
			arr[i - 1][j - 1] = 0;
			s += arr[j - 1][m - 1] * ved;
		}
		arr[i - 1][m - 1] = arr[i - 1][m - 1] - s;
	}
	
	for (int i = 0; i < N; i++) {
		x[i] = arr[i][N];
	}

}



/**
 * Метод прогонки
 */
void golovatyukam::lab3()
{
	int n =N;
	
	double *alfa = new double[n];
	double *betta = new double[n];
	double *F = new double[n];

	alfa[0] = -A[0][1] / A[0][0];
	betta[0] = b[0] / A[0][0];
	F[0] = -b[0];

	for (int i = 1; i < n - 1; i++) {
		F[i] = -1 * b[i];

		alfa[i] = -A[i][i + 1] / (A[i][i - 1] * alfa[i - 1] + A[i][i]);
		betta[i] = (-A[i][i - 1] * betta[i - 1] - F[i]) / (A[i][i - 1] * alfa[i - 1] + A[i][i]);
	}

	betta[n - 1] = (-A[n - 1][n - 2] * betta[n - 2] + b[n - 1]) / (A[n - 1][n - 2] * alfa[n - 2] + A[n - 1][n - 1]);

	x[n - 1] = betta[n - 1];

	for (int i = n - 2; i >= 0; i--) {
		x[i] = betta[i] + alfa[i] * x[i + 1];
	}
}



/**
 * Метод простых итераций
 */
void golovatyukam::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void golovatyukam::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void golovatyukam::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void golovatyukam::lab7()
{

}


void golovatyukam::lab8()
{

}


void golovatyukam::lab9()
{

}


std::string golovatyukam::get_name()
{
  return "Golovatyukam A.M.";
}
