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

}



/**
 * Метод простых итераций
 */
void zevaykinae::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void zevaykinae::lab5()
{

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
