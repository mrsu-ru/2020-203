#include "gorbunovaa.h"


void gorbunovaa::lab1()
{
  cout << "Hello world!" << endl;
}

void gorbunovaa::lab2()
{
	const double e = 1e-10;
	int n = N;

	for (int i = 0; i < n; i++) {
		if (fabs(A[i][i]) < e) { // Проверка диагонали на элемент меньше е (т.е. нулевой)
			for (int j = i + 1; j < n; j++) {
				if (fabs(A[j][i]) > e) { // Если найден такой эл., то заменям на строку с не нулевым эл.
					double *tmp = A[i]; 
					for (int k = 0; k < n; k++) { 
						A[i][k] = A[j][k];
						A[j][k] = tmp[k];
					}

					double temp = b[i];
					b[i] = b[j];
					b[j] = temp;
					break;
				}
			}
		}

		if (fabs(A[i][i] - 1) > e) {		// Сделанно для оптимизации, 
			for (int j = i + 1; j < n; j++) // что бы элемент равный 1, лишний раз не делил все остальные эл.
				A[i][j] /= A[i][i];
			b[i] /= A[i][i];
			A[i][i] = 1;
		}

		for (int j = 0; j < i; j++) {
			for (int k = i + 1; k < n; k++)
				A[j][k] -= A[i][k];
			b[j] -= b[i] * A[j][i];
			A[j][i] = 0;
		}

		for (int j = i + 1; j < n; j++) {
			for (int k = i + 1; k < n; k++)
				A[j][k] -= A[i][k] * A[j][i];
			b[j] -= b[i] * A[j][i];
			A[j][i] = 0;
		}

		for (int j = i +1; j < n; j++) {
			for (int k = i + 1; k < n; k++)
				A[j][k] -= A[i][k] * A[j][i];
			b[j] -= b[i] * A[j][i];
			A[j][i] = 0;
		}
	}

	for (int i = 0; i < n; i++)
		x[i] = b[i];

	for (int i = 0; i < n; i++) {
		cout << "x[" << i << "] = " << x[i] << endl;
	}
}


void gorbunovaa::lab3()
{
	int n = N;
	double *a = new double[n];
	double *b = new double[n];
	double **matr = A;
	

	// cout << matr[0][1] << endl;

	// Прямой ход
	a[1] = -matr[0][1] / matr[0][0]; 
	b[1] = matr[0][n + 1] / matr[0][0];
	//F[0] = -b[0];

	for (int i = 1; i < n - 1; i++) {
		//F[i] = -b[i];
		a[i + 1] = -matr[i][i + 1] / (matr[i][i] + matr[i][i - 1] * a[i]);
		b[i + 1] = (-matr[i][i - 1] * b[i] + matr[i][n + 1]) / (matr[i][i] + matr[i][i - 1] * a[i]);
	}

	// Обрантый ход
	double *x = new double[n];

	x[n - 1] = (-matr[n - 1][n - 2] * b[n - 1] + matr[n - 1][n + 1]) / (matr[n - 1][n - 1] + matr[n - 1][n - 2] * a[n - 1]);

	for (int i = n - 2; i >= 0; i--) {
		x[i] = a[i + 1] * x[i + 1] + b[i + 1];
	}

	// Вывод
	for (int i = 0; i < n; i++) {
		cout << "x[" << i + 1 << "] = " << x[i] << endl;
	}
}

void gorbunovaa::lab4()
{
	double **L = new double*[N];
	for (int i = 0; i < N; i++) {
		L[i] = new double[N];
		for (int j = 0; j < N; j++)
			L[i][j] = 0;
	}
	double S = 0;


	for (int i = 0; i < N; i++) {

		for (int k = 0; k < i; k++) {
			S += L[i][k] * L[i][k];
		}

		L[i][i] = sqrt(A[i][i] - S);
		S = 0;

		for (int k = i + 1; k < N; k++) {
			for (int m = 0; m < i; m++) {
				S += L[k][m] * L[i][m];
			}
			L[k][i] = (A[k][i] - S) / L[i][i];
			S = 0;
		}
	}

	double **Lt = new double*[N];
	for (int i = 0; i < N; i++) {
		Lt[i] = new double[N];
		for (int j = 0; j < N; j++)
			Lt[i][j] = 0;
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			Lt[i][j] = L[j][i];
		}
	}

	cout << "A = L * Lt = " << endl;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << L[i][j] << "  ";
		}
		cout << endl;
	}
	cout << endl;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << Lt[i][j] << "  ";
		}
		cout << endl;
	}
	cout << endl;

	double *y = new double[N];

	for (int i = 0; i < N; i++)
	{
		S = 0;
		for (int k = 0; k < i; k++) {
			S += L[i][k] * y[k];
		}
		y[i] = (b[i] - S) / L[i][i];
	}


	/*for (int i = 0; i < N; i++) {
		cout << "y[" << i << "] = " << y[i] << endl;
	} */

	double *x = new double[N];
	for (int k = N - 1; k >= 0; k--)
	{
		double res = 0;
		for (int i = k + 1; i < N; i++)
			res += Lt[k][i] * x[i];
		x[k] = (y[k] - res) / Lt[k][k];
	}

	for (int i = 0; i < N; i++) {
		cout << "x[" << i << "] = " << x[i] << endl;
	}
}




void gorbunovaa::lab5()
{
	double x[N], y[N];
	for (int i = 0; i < N; i++)
	{
		x[i] = 0; // начальное приближение
	}
	double temp = 0.0;
	double eps = 1e-20;
	int k = 0;
	do
	{
		k++;
		temp = 0.0;
		for (int i = 0; i < N; i++)
			y[i] = x[i];
		for (int i = 0; i < N; i++)
		{
			double s = 0;
			for (int j = 0; j < i; j++)
				s += A[i][j] * x[j];
			for (int j = i + 1; j < N; j++)
				s += A[i][j] * y[j];
			x[i] = (b[i] - s) / A[i][i];
		}
		for (int i = 0; i < N; i++)
		{
			if (fabs(y[i] - x[i]) > temp)
				temp = fabs(y[i] - x[i]);
		}

	} while (temp >= eps);

}




void gorbunovaa::lab6()
{
	double eps = 1e-9;

	for (int i = 0; i < N; i++) {
		double maxel = A[i][i];
		int indRow = i;
		for (int j = i + 1; j < N; j++)
			if (maxel < abs(A[j][i])) {
				indRow = j;
				maxel = abs(A[j][i]);
			}


		if (indRow != i) {
			for (int j = i; j < N; j++) {
				swap(A[i][j], A[indRow][j]);
			}
			swap(b[i], b[indRow]);
		}

		double summ = 0;
		for (int j = 0; j < N; j++)
			summ += abs(A[i][j]);
		if (2 * abs(A[i][i]) < summ) {
			cout << "Error" << endl;
		}

		maxel = A[i][i];
		b[i] /= A[i][i];
		A[i][i] = 0;
		for (int j = 0; j < N; j++)
			if (j != i)
				A[i][j] /= maxel;

	}

	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}

	double xx = b[0];
	double *results = new double[N];
	do {
		for (int i = 0; i < N; i++) {
			results[i] = 0;
			for (int k = 0; k < N; k++)
				results[i] -= A[i][k] * x[k];
			results[i] += b[i];
		}
		xx = x[0];
		for (int i = 0; i < N; i++) {
			x[i] = results[i];
		}
	} while (abs(x[0] - xx) > eps);
}



void gorbunovaa::lab7()
{

}


void gorbunovaa::lab8()
{

}


void gorbunovaa::lab9()
{

}


std::string gorbunovaa::get_name()
{
  return "Gorbunov A. A.";
}