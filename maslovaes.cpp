#include "maslovaes.h"

/**
 * Введение в дисциплину
 */
void maslovaes::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void maslovaes::lab2()
{
for (int k = 0; k < N; k++) {
	int idmax = 0;
	for (int i=0; i<N; i++) {
		if(abs(A[i][k]) > abs(A[idmax][k])) idmax = i; 
	}

	for (int i = 0; i < N; i++)
	{
		swap(A[idmax][i], A[k][i]);
	}
	swap(b[idmax], b[k]);

	//down
	for (int i = k + 1; i < N; i++)
	{
		double tmp = A[i][k]/A[k][k];
		for (int j = k; j < N; j++) {
			A[i][j] -= A[k][j]*tmp;
		}
		b[i] -= b[k]*tmp;
	}
}

for(int i = 0; i<N; i++){
    x[i]=b[i];
}

//up
for (int i = N-1; i >= 0; i--){
	for (int j = i+1; j < N; j++) {
		x[i] -= A[i][j]*x[j];
	}
	x[i] /= A[i][i];
}
}



/**
 * Метод прогонки
 */
void maslovaes::lab3()
{
double* alpha = new double[N];
double* beta = new double[N];
alpha[0]= - A[0][1]/A[0][0];
beta[0]=b[0]/A[0][0];
for(int i=1; i<N; i++){
    alpha[i] = - A[i][i + 1]/(A[i][i] + A[i][i - 1]*alpha[i - 1]);
    beta[i] = (- A[i][i - 1]*beta[i - 1] + b[i])/(A[i][i] + A[i][i - 1]*alpha[i - 1]);
}
x[N - 1] = beta[N - 1];
for(int i=N - 2; i>=0; i--){
    x[i]=alpha[i]*x[i + 1] + beta[i];
}
}



/**
 * Метод Холецкого
 */
void maslovaes::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void maslovaes::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void maslovaes::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void maslovaes::lab7()
{

}


void maslovaes::lab8()
{

}


void maslovaes::lab9()
{

}


std::string maslovaes::get_name()
{
  return "Maslova E.S.";
}
