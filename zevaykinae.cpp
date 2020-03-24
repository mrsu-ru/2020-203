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
double eps = 1e-15; //чтобы выйти из цикла
double t = 1e-5; //(приближенный парметр) достаточно малое число, ПРИ ПОИМОМЩИ КОТОРОГО МЫ ПОЛУЧАЕМ КАЖДЫЙ РАЗ ВСЁ БОЛЕЕ И БОЛЕЕ ТОЧНОЕ ЧИСЛО

for (int i = 0; i < N; i++) //берем первое приближенное значение и прогоняем его
    {
		x[i] = 0;
	}

	double x1;
	double *xr = new double[N];//рассчитываем новые значения x-ов
	int step = 0;

	do {
		step++;

		for (int i = 0; i < N; i++) //рассчитываем новое точное значение в цикле
        {
			xr[i] = x[i];
			for (int k = 0; k < N; k++)
				xr[i] -= t*A[i][k] * x[k];
			xr[i] += t * b[i];

		}

		x1 = 0;
		for (int i = 0; i < N; i++) { //рассчитываем норму (с её помощью ведем учет того, что значение достаточно близкое к необходимое)
			x1 += (xr[i]-x[i])*(xr[i]-x[i]);
		}

		for (int i = 0; i < N; i++)
        {
			x[i] = xr[i]; //осуществляем перепресвоение, новое становится старым и мы возобновляем цикл
		}
	}
	while (sqrt(x1)>eps); //как меньш8е, то выходим из цимкла
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
