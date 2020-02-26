#include "kozinasa.h"



/**
 * Введение в дисциплину
 */
void kozinasa::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kozinasa::lab2(){
  double y;
  for (int k=0; k<N; k++) {
	int mEl=k;
	for(int i=k+1;i<N;i++)
	  if(abs(A[i][k]) > abs(A[mEl][k])) mEl=i;
	for(int i=0;i<N;i++)
	std::swap(A[k][i],A[mEl][i]);
	std::swap(b[k],b[mEl]);

	y = A[k][k];
	for (int j=0; j<N; j++)
	  A[k][j] = A[k][j] / y;
    b[k] = b[k]/y;

    for (int i=k+1; i<N; i++){
	  y = A[i][k];
	  for (int j=0; j< N; j++){
		  A[i][j] =A[i][j]- A[k][j] * y;
	  }
    b[i] =b[i]- b[k] * y;
    }
  }

  for (int k=N-1; k>0; k--){
  for (int i=k-1; i>=0; i--){
    y = A[i][k];

    for (int j=0; j<N; j++)
      A[i][j] = A[i][j] - A[k][j] * y;
    b[i] = b[i] - b[k] * y;
    }
  }

  for(int i=0; i<N; i++)
    x[i] = b[i];
}


/**
 * Метод прогонки
 */
void kozinasa::lab3()
{
	int n = N, i;
  double *P, *Q;
  P = new double[n];
  Q = new double[n];

  P[0]=A[0][1]/-A[0][0];
  Q[0]=-b[0]/-A[0][0];
  cout << P[0] << " " << Q[0] << endl;
  for (i=1;i<n-1;i++){
  	P[i]=A[i][i+1]/(-A[i][i] - A[i][i-1]*P[i-1]);
  	Q[i]=(A[i][i-1]*Q[i-1] - b[i])/(-A[i][i] - A[i][i-1]*P[i-1]);
  
  }
   x[n-1] = (A[i][i-1]*Q[i-1] - b[i])/(-A[i][i] - A[i][i-1]*P[i-1]);
   for (int i=n-2;i>=0;i--){
   	x[i] = P[i]*x[i+1] + Q[i];
   }
}



/**
 * Метод простых итераций
 */
void kozinasa::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void kozinasa::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void kozinasa::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void kozinasa::lab7()
{

}


void kozinasa::lab8()
{

}


void kozinasa::lab9()
{

}


std::string kozinasa::get_name()
{
  return "Kozina S.A.";
}
