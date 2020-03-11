#include "zhalninrv.h"
#include "puzinva.h"

/**
 * Введение в дисциплину
 */
void puzinva::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void puzinva::lab2()
{

}



/**
 * Метод прогонки
 */
void puzinva::lab3()
{
	double k1[N - 1];
	double k2[N];
    k1[0] =- A[0][1] / A[0][0];
    k2[0] = b[0] / A[0][0];

    for (int i = 1; i < N - 1; i++) {
        k1[i] =- A[i][i + 1] / (A[i][i] + A[i][i - 1] * k1[i - 1]);
        k2[i] = (b[i] - A[i][i - 1] * k2[i - 1]) / (A[i][i] + A[i][i - 1] * k1[i - 1]);
    }
	
    k2[N - 1] = (b[N - 1] - A[N - 1][N - 2] * k2[N - 2]) / (A[N - 1][N - 1] + A[N - 1][N - 2] * k1[N - 2]);
    x[N - 1] = k2[N - 1];
	
    for (int i = N - 2; i >= 0; i--) {
        x[i] = k1[i] * x[i + 1] + k2[i];
    }
}



/**
 * Метод простых итераций
 */
void puzinva::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void puzinva::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void puzinva::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void puzinva::lab7()
{

}


void puzinva::lab8()
{

}


void puzinva::lab9()
{

}


std::string puzinva::get_name()
{
  return "Puzin V.A.";
}
