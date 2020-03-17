//#include "zhalninrv.h"
#include "edelevaup.h"

/**
 * Введение в дисциплину
 */
void edelevaup::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void edelevaup::lab2()
{

	for (int i = 0; i < N; i++) {

		int max = i;
		for (int j = i + 1; j < N; j++) {
			if (abs(A[j][i]) > abs(A[max][i]))
				max = j;
		}
		double *s;
		double k;
		if (max != i) {
			s = A[i];
			A[i] = A[max];
			A[max] = s;
			k = b[i];
			b[i] = b[max];
			b[max] = k;
		}
		
		b[i] /= A[i][i];
		for (int j = N - 1; j > i; j--)
			A[i][j] /= A[i][i];
		A[i][i] = 1;

		for (int j = i + 1; j < N; j++) {

			for (int k = N; k > i; k--) {
				A[j][k] -= A[i][k] * A[j][i];
			}
			b[j] -= b[i] * A[j][i];

			A[j][i] = 0;
		}


	}

	for (int i = N - 1; i > 0; i--) {
		double s = 0, l = 0;
		for (int j = i; j < N; j++) {
			s += b[j ] * A[i - 1][j];
			A[i - 1][j] = 0;
			
		}
		b[i - 1] -= s;
		
	}

	for (int i = 0; i < N; i++) {
		x[i] = b[i];
		
		}
	

}



/**
 * Метод прогонки
 */
void edelevaup::lab3()
{
	double *al = new double[N];
    double *bet = new double[N];
    al[0]=A[0][1]/-A[0][0]; 
	bet[0]=b[0]/A[0][0];
	for( int i=1; i<N; i++){
		al[i]=A[i][i+1]/(-A[i][i]-A[i][i-1]*al[i-1]);
		bet[i]=(A[i][i-1]*bet[i-1]-b[i])/(-A[i][i]-A[i][i-1]*al[i-1]);
}
 
  x[N-1]=bet[N-1];
  for (int i=N-2; i>=0; i--){
	  x[i]=al[i]*x[i+1]+bet[i];
  }
 
}
/**
 * Метод простых итераций
 */
void edelevaup::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void edelevaup::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void edelevaup::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void edelevaup::lab7()
{

}


void edelevaup::lab8()
{

}


void edelevaup::lab9()
{

}


std::string edelevaup::get_name()
{
  return "Edeleva U.P.";
}
