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
	double eps = 0.00001;
    int lead;
    double maxElem; 
	double keep;
    int change[N];
	
    for (int i = 0; i < N; i++) {
        lead = i;
        for (int j = i + 1; j < N; j++) {
            if (abs(A[j][i]) > abs(A[lead][i])) {
                lead = j;
            }
        }
		
        if (abs(A[lead][i]) < eps) {
            cout << "Нет решений";
            break;
        }
		
        change[i]=lead;
		
        if (lead != i) {
            swap(A[i], A[lead]);
            swap(b[i], b[lead]);
        }
		
        maxElem = A[i][i];
        A[i][i] = 1;
		
        for (int j = i + 1; j < N; j++) {
            A[i][j] /= maxElem;
        }
		
        b[i] /= maxElem;
		
        for (int j = i + 1; j < N; j++) {
            keep = A[j][i];
            A[j][i]=0;
			
            for (int k = i + 1; k < N; k++) {
                A[j][k] -= A[i][k] * keep;
            }
			
            b[j] -= b[i] * keep;
        }
    }
	
    for (int i = N - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            keep = A[j][i];
            A[j][i] = 0;
			
            for (int k =i - 1; k >= 0; k--) {
                A[j][k] -= A[i][k] * keep;
            }
			
            b[j] -= b[i] * keep;
        }
    }
	
    for (int i = 0; i < N; i++) {
        x[i] = b[i];
    }
	
    for (int i = N - 1; i >= 0; i--) {
        if (change[i] != i) {
            swap(A[change[i]], A[i]);
            swap(b[change[i]], b[i]);
        }
    }
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
 * Метод Холецкого
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
