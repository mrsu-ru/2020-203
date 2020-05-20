#include "gorbunovaa.h"

/**
 * Ââåäåíèå â äèñöèïëèíó
 */
void gorbunovaa::lab1()
{
  cout << "hello world!" << endl;
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
 * Ìåòîä ïðîãîíêè
 */
void gorbunovaa::lab3()
{

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