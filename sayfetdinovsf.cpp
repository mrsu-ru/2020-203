#include "sayfetdinovsf.h"

/**
 * Введение в дисциплину
 */
void sayfetdinovsf::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void sayfetdinovsf::lab2()
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
void sayfetdinovsf::lab3()
{
	double *p = new double[N];
	double *q = new double[N];

	// Прямой ход
	p[0] = -A[0][1] / A[0][0];
	q[0] = b[0] / A[0][0];
	for (int i = 1; i < N; i++) {
		p[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * p[i - 1]);
		q[i] = (b[i] - A[i][i - 1] * q[i - 1]) / (A[i][i] + A[i][i - 1] * p[i - 1]);
	}

	// Обратный ход
	x[N - 1] = q[N - 1];
	for (int i = N - 2; i >= 0; i--) {
		x[i] = p[i] * x[i + 1] + q[i];
	}

	delete[] p;
	delete[] q;
}



/**
 * Метод простых итераций
 */
void sayfetdinovsf::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void sayfetdinovsf::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void sayfetdinovsf::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void sayfetdinovsf::lab7()
{

}


void sayfetdinovsf::lab8()
{

}


void sayfetdinovsf::lab9()
{

}


std::string sayfetdinovsf::get_name()
{
  return "sayfetdinovsf";
}
