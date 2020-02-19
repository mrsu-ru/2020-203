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
