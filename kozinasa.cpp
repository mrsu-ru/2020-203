#include "kozinasa.h"


void kozinasa::lab1()
{
}
/**
 * Метод Гаусса с выбором главного элемента
 */
void kozinasa::lab2()
{
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

void kozinasa::lab3(){
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
   * Метод квадратного корня (метод Холецкого)
   */
void kozinasa::lab4(){ 
int n=N; double sum=0;
  double *y;
  double** L; 
  
  y = new double[n];

  L = new double*[n];
  
  for (int i = 0; i < n; i++){
    L[i] = new double[n];
    for (int j=0; j<n; j++){
	  L[i][j] = 0;
  	}
  	y[i] = 0;
  	x[i] = 0;
  }
  
  
  for (int i=0; i<n; i++){
  	for(int k = 0; k <= i; k++){
        sum += L[i][k] * L[i][k];
    }
    L[i][i] = sqrt(A[i][i] - sum);
    sum = 0;
    for(int j = i + 1; j < n; j++){
        for(int k = 0; k <= j - 1; k++){
            sum += L[j][k] * L[i][k];
        }
		L[j][i] = (A[j][i] - sum) / L[i][i];
        sum = 0;
    }
  }
  
  y[0] = b[0]/L[0][0];
  
  for (int i=1; i<n; i++){
  	y[i] = b[i];
  	for (int j=0; j<i; j++){
  		y[i] -= L[i][j]*y[j];
	  }
	y[i] /= L[i][i];
  }

  x[n-1] = y[n-1]/L[n-1][n-1];
  
  for (int i=n-2; i>=0; i--){
  	x[i] = y[i];
  	for (int j=n-1; j>i; j--){
  		x[i] -= L[j][i]*x[j];
	  }
	x[i] /= L[i][i];
  }
}



/**
 * Метод Холецкого
 */

void kozinasa::lab5()
{ int n=N; double e = pow(10,-35); //0.0000000000000000000000000000000000000000001;
  for (int i=0;i<n;i++){
  	x[i]=b[i];
  }
	
	double xxx = 0;
  while (abs(xxx-x[n-1])>e){
  	xxx=x[n-1];
  	for (int i=0;i<n;i++){
  	  double s1=0, s2=0;
  	  for (int j=0; j<i;j++)
  		s1 += A[i][j]*x[j];
  	  for (int j=i+1; j<n;j++)
  		  s2 += A[i][j]*x[j];
	  x[i] = (b[i] - s1 - s2)/A[i][i];
	  }
	  
  }
  
}



/**
 * Метод Якоби или Зейделя
 */
void kozinasa::lab6(){
 int n=N,i,j; double e = 1.0e-19;
 double *r, *h, *x1;
 
  x1 = new double[n];
  r = new double[n];
  h = new double[n];
  
  for (int i=0;i<n;i++){
  	x[i]=b[i];
  	x1[i]=0;
  }
  
  
  double xxx = 1;
  while (xxx>e){
	//���������� rn
	for (i=0;i<n;i++){
  	  double s = 0;
  	  for (j=0;j<n;j++){
  		 s += A[i][j]*x[j];
	   }
	  r[i]=s-b[i];
      }
    for (i=0;i<n;i++){
  	  double s = 0;
  	  for (j=0;j<n;j++){
  		 s += A[i][j]*r[j];
	   }
	  h[i]=s;
	  //cout<<h[i]<<" ";
      } //cout<<endl;
    //��������� ��� ���������� <���>
    double sum1 = 0;
  	for (i=0;i<n;i++){
  		sum1 += r[i]*h[i];
  	}
  	//����������� ��� ���������� <���>
  	double sum2 = 0;
  	for (i=0;i<n;i++){
  		sum2 += h[i]*h[i];
  	}
  	double ta = sum1/sum2;
  	for(i=0;i<n;i++){
  		x[i]=x[i]-ta*r[i];
		  //cout<<x[i]<<" ";
	  }
	//cout<<endl;
	xxx=0;
	for (i=0;i<n;i++){
		xxx+=(x1[i]-x[i])*(x1[i]-x[i]);
  		x1[i]=x[i];
  	}
  	xxx = sqrt(xxx);
  }
  
}



/**
 * Метод минимальных невязок
 */
void kozinasa::lab7()
{

}



/**
 * Метод сопряженных градиентов
 */
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
