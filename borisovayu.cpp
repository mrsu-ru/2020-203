#include "borisovayu.h"

/**
 * Введение в дисциплину
 */
void borisovayu::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void borisovayu::lab2()
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
void borisovayu::lab3(){
	double *alpha;
	alpha = new double[N];
	double *beta;
	beta = new double[N];
	int i;
	// прямая прогонка
	alpha[0] = A[0][1]/A[0][0];
	beta[0] = b[0]/A[0][0];
	for (i=1; i<N; i++){
		alpha[i] = A[i][i+1] / (A[i][i] - A[i][i-1]*alpha[i-1]);
		beta[i] = (b[i] - beta[i-1] * A[1][0]) / (A[i][i] - A[i][i-1]*alpha[i-1]) ;		
	}
	
	x[N-1] = beta[N-1];
	for (i=N-2; i>=0; i--){
		x[i] = beta[i] - alpha[i]*x[i+1];
	}
	
}



/**
 * Метод Холецкого
 */
void borisovayu::lab4()
{
     int i,j,k;
     double tmp, tmp2;
	 
     int *D;
     D = new int[N];
	 
     double **S;
     S = new double*[N];
     for (i=0; i<N; i++){
	S[i] = new double[N];
	for (j=0; j<N; j++) S[i][j] = 0;
     }
	 
     for (i=0; i<N; i++){
        tmp = A[i][i];
	for (j=0; j<i; j++) tmp-=D[j]*S[j][i]*S[j][i];
		 
	if (tmp>0) D[i] = 1;
	else       D[i] =-1;
		 
	S[i][i] = sqrt(D[i]*tmp);
		
	for (j=i+1; j<N; j++){
		tmp2 = A[i][j];
		for (k=0; k<j; k++) tmp2 -= D[k]*S[k][i]*S[k][j];
		S[i][j] = tmp2 / (D[i] * S[i][i]);
	}
      }
	 
     double *y;
     y = new double[N];
     y[0] = b[0]/S[0][0];
	 
     for (i=1; i<N; i++){
	for (j=i; j<N; j++) b[j]-=S[i-1][j]*y[i-1];
	y[i] = b[i]/S[i][i];
     }
	 
     for (i=0; i<N; i++)
       for (j=0; j<N; j++)
         S[i][j] *= D[i];

	 
     x[N-1] = y[N-1]/S[N-1][N-1];
	 
     for (i=N-2; i>=0; i--){
	for (j=i ; j>=0; j--) y[j]-=S[j][i+1]*x[i+1];
	x[i] = y[i]/S[i][i];
     }	 
	 
	 
     delete []y;
     for (int i = 0; i < N; i++)  
	delete []S[i];
     delete []S;
     delete []D;
}



/**
 * Метод Якоби или Зейделя
 */
void borisovayu::lab5()
{
    int i,j,k;
    double *z;
    z = new double[N];
    for (i=0; i<N; i++)
	    z[i]=b[i];
   
    for (i=0; i<N; i++)
	    for (j=0; j<N; j++) 
		    z[i]-=A[i][j]*x[j];
		
    double eps = 1e-17;
    double Z=0;
	
    for (i=0; i<N; i++)
	Z += z[i]*z[i];
    
    while (Z > eps * eps){
	for (i=0; i<N; i++){
		x[i] = b[i];
	        for (j=0; j<i; j++)
			x[i] -= A[i][j]*x[j];
		for (j=i+1; j<N; j++)
			x[i] -= A[i][j]*x[j];
		x[i]/=A[i][i];
		}
		
    
    for (i=0; i<N; i++)
	    z[i] = b[i];
    for (i=0; i<N; i++)
	    for (j=0; j<N; j++) 
		    z[i] -= A[i][j]*x[j];		
		
    Z=0;
    for (i=0; i<N; i++)
	Z += z[i]*z[i];    	
    }
	
    delete []z;
}



/**
 * Метод минимальных невязок
 */
void borisovayu::lab6()
{
     double eps = 1e-17;    
     double *z, *Az;
     int i,j;
	
     z = new double[N];
     Az = new double[N];
     for (i=0; i<N; i++) x[i]=0;
     for (i=0; i<N; i++) z[i]=b[i];
	
     for (i=0; i<N; i++)
	    for (j=0; j<N; j++) 
		    z[i]-=A[i][j]*x[j];

     double Z=0;
     double t=0;
     double scalar=0;
	
     for (i=0; i<N; i++)
	Z += z[i]*z[i];
    
     while (Z > eps*eps){
	for (i=0; i<N; i++){
	    Az[i] = 0;
	    for (j=0; j<N; j++)	Az[i] += A[i][j]*z[j];
	}
		
	scalar = 0;
	for (i=0; i<N; i++) scalar += Az[i]*z[i];
	t = -scalar;
	scalar = 0;
	for (i=0; i<N; i++) scalar += Az[i]*Az[i];
	t/=scalar;
	
	for (i=0; i<N; i++) x[i] = x[i] - t*z[i];
			
	for (i=0; i<N; i++)  z[i] = b[i];
        for (i=0; i<N; i++)
	        for (j=0; j<N; j++) 
		        z[i] -= A[i][j]*x[j];		
		
	Z=0;
	for (i=0; i<N; i++)  Z += z[i]*z[i];
	}
	
	delete[] Az;
	delete[] z;
}



/**
 * Метод сопряженных градиентов
 */
void borisovayu::lab7()
{

}


void borisovayu::lab8()
{

}


void borisovayu::lab9()
{

}


std::string borisovayu::get_name()
{
  return "Borisov A. Yu.";
}
