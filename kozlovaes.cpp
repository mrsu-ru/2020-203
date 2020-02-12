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
cout << "\nПрямой ход: " << endl;
	for (int i = 0; i < N; i++) {
		int max = i;
		for (int j = i + 1; j < N; j++) {
			if (abs(A[j][i]) > abs(A[max][i]))
				max = j;
		}

		if (abs(A[i][i]) <= eps) {
			cout<<"Нет решений!";
		}

		double *s;
		double tmp;
		if (max != i) {
			s = A[i]; A[i] = A[max]; A[max] = s;
			tmp = b[i]; b[i] = b[max]; b[max] = tmp;
		}
		cout << endl;
		print(N, A);
		print1(N, b);

		for (int j = N; j > i; A[i][j--] /= A[i][i]);
		b[i] /= A[i][i];
		A[i][i] = 1;
		cout << endl;
		print(N, A);
		print1(N, b);

		for (int j = i + 1; j < N; j++) {

			for (int k = N; k > i; A[j][k] -= A[i][k--] * A[j][i]);
			b[j] -= b[i] * A[j][i];
			A[j][i] = 0;
		}
		print1(N, b);

		cout << endl;
		print(N, A);
	}
	cout << "\n\nОбратный ход: " << endl;
	for (int i = N - 1; i > 0; i--) {
		double s = 0;
		for (int j = i ; j < N; j++) {
			s += b[j] * A[i - 1][j];
			A[i - 1][j] = 0;
		}
		b[i - 1] -=s;
		cout << endl;
		print1(N, b);
		print(N, A);
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

}



/**
 * Метод простых итераций
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
