#include "borisovayu.h"
#include "cmath"

using namespace std;

double scal(double *z1, double *z2, int N)
{
	double result = 0;
	for (int i=0; i<N; i++) 
		result += z1[i]*z2[i];
	return result;
}
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
    double *z;
	double *r;
	double *Az;
	double alpha, beta, eps = 1e-17;
	int i,j, sum = 1, k=0;

    z = new double[N];
    r = new double[N];
	Az = new double[N];
	
	for (i=0; i<N; i++){
		r[i] = b[i];
        x[i] = 0;
        z[i] = r[i];		
	}
	sum = scal(r, r, N);
    do {
		for (i=0; i<N; i++){
			Az[i] = 0;
			for (j=0; j<N; j++)  Az[i] += A[i][j] * z[j];	
		}
		alpha = scal(r, r, N)/scal(Az, z, N);
		for (i=0; i<N; i++) x[i] = x[i] + alpha * z[i];	
		for (i=0; i<N; i++) r[i] = r[i] - alpha * Az[i];
		beta = sum/scal(r, r, N);
		for (i=0; i<N; i++) z[i] = r[i] + beta * z[i];
		sum = scal(r, r, N);
		k++;	
	} while (scal(r,r,N) > eps*eps);
	
}


void borisovayu::lab8()
{
    double err = 0;
	double eps = 1e-6;
	double max;
	double alpha;
	double c,s;
	int i, j, I, J;
    for (int i = 0; i < N; i++) 
        for (int j = i + 1; j < N; j++)
           err += A[i][j] * A[i][j] + A[j][i] * A[j][i];
		
    double **C;
    C = new double*[N];
    for (i=0; i<N; i++){
	    C[i] = new double[N];
    }

    while (err > eps){
		alpha = 0;
		max = abs(A[0][1]);
		I = 0;
		J = 1;
		for (i = 0; i<N; i++)
			for (j = i+1; j<N; j++)
				if (abs(A[i][j]) > max){
					max = abs(A[i][j]);
					I = i;  J = j;
			    }
				
		if (abs(A[I][I] - A[J][J]) > eps) alpha = atan(2*A[I][J] / (A[J][J] - A[I][I])) / 2;
		else alpha = atan(1);
		
        c = cos(alpha);
        s = sin(alpha);
		
		for (i=0; i<N; i++)
			for (j=0; j<N; j++)
				if (i!=I && j!=I && i!=J && j!=J) C[i][j] = A[i][j];		
		
        C[I][I] = c*c*A[I][I] - 2*s*c*A[I][J] + s*s*A[J][J];
        C[J][J] = s*s*A[I][I] + 2*s*c*A[I][J] + c*c*A[J][J];
        C[I][J] = (c*c - s*s)*A[I][J] + s*c*(A[I][I] - A[J][J]);
		C[J][I] = C[I][J];
		for (int k=0; k<N; k++)
			if (k!=I && k!=J){
				C[I][k] = c*A[I][k] - s*A[J][k];
				C[k][I] = C[I][k];
				C[J][k] = s*A[I][k] + c*A[J][k];	
				C[k][J] = C[I][k];				
			}
				
		for (i=0; i<N; i++)
			for (j=0; j<N; j++)
				A[i][j] = C[i][j];
			
	    err = 0;
        for (i = 0; i < N; i++) 
           for (j = i + 1; j < N; j++)
              err += C[i][j] * C[i][j] + C[j][i] * C[j][i];

	}

    for (i=0; i<N; i++) x[i] = A[i][i];	
	
	for (i = 0; i < N; i++) delete []C[i];
	delete []C;
}


void borisovayu::lab9()
{
	double eps = 1e-3; 
	double *yk = new double[N];
	double *yk1 = new double[N]; 
	int i,j;

	for (int i = 0; i < N; i++)	yk[i] = 1; 

	double err = 0; 
	double lambda1 = 0; 
	double lambda2 = 0;

	do{
		for (i = 0; i < N; i++) 
			for (j = 0; j < N; j++) 
				yk1[i] += A[i][j] * yk[j]; 
			
		for (i = 0; i < N; i++)
			if(fabs(yk[i]) > eps && fabs(yk1[i]) > eps){
				lambda1 = yk1[i]/yk[i];
				break;
			}
		
		err = fabs(lambda1 - lambda2); 
		lambda2 = lambda1; 

		for (i = 0; i < N; i++) {
			yk[i] = yk1[i];
            yk1[i] = 0;			
		}
	}
	while(err > eps);
	cout<<"max lambda = "<<lambda2<<endl;
}


std::string borisovayu::get_name()
{
  return "Borisov A. Yu.";
}
