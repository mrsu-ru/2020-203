#include "gorbunovaa.h"

/**
 * Ââåäåíèå â äèñöèïëèíó
 */
void gorbunovaa::lab1()
{
  cout << "Hello world!" << endl;
}


/**
 * Ìåòîä Ãàóññà ñ âûáîðîì ãëàâíîãî ýëåìåíòà
 */
void gorbunovaa::lab2()
{
	const double e = 1e-10;

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



/**
 * Метод прогонки
 */
void gorbunovaa::lab3()
{
	float *a = new float[n];
	float *b = new float[n];

	cout << matr[0][1] << endl;

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
	float *x = new float[n];

	x[n - 1] = (-matr[n - 1][n - 2] * b[n - 1] + matr[n - 1][n + 1]) / (matr[n - 1][n - 1] + matr[n - 1][n - 2] * a[n - 1]);

	for (int i = n - 2; i >= 0; i--) {
		x[i] = a[i + 1] * x[i + 1] + b[i + 1];
	}

	// Вывод
	for (int i = 0; i < n; i++) {
		cout << "x[" << i + 1 << "] = " << x[i] << endl;
	}
}



/**
 * Ìåòîä ïðîñòûõ èòåðàöèé
 */
void gorbunovaa::lab4()
{

}



/**
 * Ìåòîä ßêîáè èëè Çåéäåëÿ
 */
void gorbunovaa::lab5()
{

}



/**
 * Ìåòîä ìèíèìàëüíûõ íåâÿçîê
 */
void gorbunovaa::lab6()
{

}



/**
 * Ìåòîä ñîïðÿæåííûõ ãðàäèåíòîâ
 */
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