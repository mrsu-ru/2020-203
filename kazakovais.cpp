#include "kazakovais.h"

/**
 * Введение в дисциплину
 */
void kazakovais::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kazakovais::lab2()
{
int i,j,k,max;

  for (i=0; i<N; b[i]=-b[i], i++)

  for (i=0; i<N; i++) {
      max=i;
      double* s; double s1;
      for (j=i+1; j<N; j++) if (A[j][i]*A[j][i]>A[max][i]*A[max][i]) max=j;
      if (max!=i) {
         s=A[i]; s1=b[i];
         A[i]=A[max]; b[i]=b[max];
         A[max]=s; b[max]=s1;
      }

      for (j=N-1; j>i; A[i][j]/=A[i][i], j--);
      b[i]/=A[i][i];
      A[i][i]=1;

      for (j=N-1; j>i; j--) {
          for (k=i+1; k<N; A[j][k]-=A[i][k]*A[j][i], k++);
          b[j]-=b[i]*A[j][i];
          A[j][i]=0;
      }


  }

  double sum;
  x[N-1]=b[N-1];
  for (i=N-2;i>=0;i--){
    sum=0;
    for (j=i+1;j<N;j++) sum+=x[j]*A[i][j];
    x[i]=b[i]-sum;
  }
	
}



/**
 * Метод прогонки
 */
void kazakovais::lab3()
{

}



/**
 * Метод простых итераций
 */
void kazakovais::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void kazakovais::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void kazakovais::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void kazakovais::lab7()
{

}


void kazakovais::lab8()
{

}


void kazakovais::lab9()
{

}


std::string kazakovais::get_name()
{
  return "Kazakova I.S.";
}
