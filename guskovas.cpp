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
void guskovas::lab3()
{

}



/**
 * Метод простых итераций
 */
void guskovas::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void guskovas::lab5()
{

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
