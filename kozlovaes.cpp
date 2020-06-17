#include "kozlovaes.h"

/**
 * Введение в дисциплину
 */
void kozlovaes::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kozlovaes::lab2()
{
	for (int i = 0; i <N; i++) {
		int max = i;
		for (int j = i + 1; j <N; j++) {
			if (abs(A[j][i]) > abs(A[max][i]))
				max = j;
		}

		double *s;
		double tmp;
		if (max != i) {
			s = A[i]; A[i] = A[max]; A[max] = s;
			tmp = b[i]; b[i] = b[max]; b[max] = tmp;
		}

		for (int j = N; j > i; A[i][j--] /= A[i][i]);
		b[i] /= A[i][i];
		A[i][i] = 1;

		for (int j = i + 1; j <N; j++) {

			for (int k = N; k > i; k--) A[j][k] -= A[i][k] * A[j][i];
			b[j] -= b[i] * A[j][i];
			A[j][i] = 0;
		}
	}
	double s;
	for (int i = N - 1; i > 0; i--) {
		s = 0;
		for (int j = i ; j <N; j++) {
			s += b[j] * A[i - 1][j];
			A[i - 1][j] = 0;
		}
		b[i - 1] -=s;
	}
	for (int i = 0; i <N; i++) {
		x[i] = b[i];
	}
}



/**
 * Метод прогонки
 */
void kozlovaes::lab3()
{
	double *alfa = new double[N];
	double *betta = new double[N];
	int i;
	alfa[0] = A[0][1]/-A[0][0];
	betta[0] = b[0]/A[0][0];
 

	for(i = 1;i <N;i++){
		alfa[i] = A[i][i+1]/(-A[i][i]-alfa[i-1]*A[i][i-1]);
		betta[i] = (-b[i]+A[i][i-1]*betta[i-1])/(-A[i][i]-alfa[i-1]*A[i][i-1]);
	}
	i=N-1;
	x[N-1] = (-b[i]+A[i][i-1]*betta[i-1])/(-A[i][i]-alfa[i-1]*A[i][i-1]);

	for(int i=N-1;i>=0;i--){
		x[i] = alfa[i]*x[i+1]+betta[i];
	}
}



/**
 * Метод Холецкого
 */
void kozlovaes::lab4()
{
double* d=new double[N];
double** S = new double*[N];
for (int i = 0; i <N; i++) { 
S[i] = new double[N];
}

if(A[0][0]>0) d[0]=1;
else d[0]=-1;
S[0][0]=sqrt(abs(A[0][0]));

for(int i=1;i<N;i++){
S[0][i]=A[0][i]/(S[0][0]*d[0]);
}

///////Вычисление матрицы S
double sumd=0;
for(int i=1; i <N; i++){
for(int k=0;k<i;k++){
sumd += pow(S[k][i],2)*d[k];
}
if((A[i][i]-sumd)>0){
	d[i]=1;
} 
else{
	d[i]=-1;
}
S[i][i]=sqrt(d[i]*(A[i][i]-sumd));
sumd=0;
double sumS=0;
for(int j=i+1;j <N;j++){
for(int k=0;k<j;k++){
sumS +=d[k]*S[k][i]*S[k][j];
}
S[i][j]=(A[i][j]-sumS)/(d[i]*S[i][i]);
sumS=0;
}
}

///////Решение уравнения S^t*y=b

double* y=new double [N];
y[0]=b[0]/S[0][0];

double sumS=0;
for(int i=1;i<N;i++){
for(int j=0;j<i;j++){
sumS +=S[j][i]*y[j];
}
y[i]=(b[i]-sumS)/S[i][i];
sumS=0;
}

////////Решение уравнения (SD)*x=y

x[N-1]=y[N-1]/(S[N-1][N-1]*d[N-1]);

double sumSDx=0;
for(int i=N-2;i>=0;i--){
for(int k=i+1;k<N;k++){
sumSDx +=S[i][k]*x[k];
}
x[i]=(y[i]-sumSDx)/(S[i][i]*d[i]);
sumSDx=0;
}
}



/**
 * Метод Якоби или Зейделя
 */
void kozlovaes::lab5()
{
double *f = new double[N];
double eps = 1.e-30;
double norm =0;

do{
	for(int i=0;i<N;i++){
	f[i]=x[i];
	}
	
	for(int i=0;i<N;i++){
		double sum1=0,sum2=0;
		for(int j=0;j<i;j++) sum1+=A[i][j]*x[j];	
		for(int j=i+1;j<N;j++) sum2+=A[i][j]*x[j];
		x[i]=(b[i]-sum1-sum2)/A[i][i];
	}
	
norm=0;
for(int i=0;i<N;i++){
	if(abs(x[i]-f[i])>norm)	norm=abs(x[i]-f[i]);
}

}while(norm>eps);

}



/**
 * Метод минимальных невязок
 */
void kozlovaes::lab6()
{
    double eps = 1.e-20;
    double *r0 = new double [N];
    double *Fi = new double [N];
	double *Ar0 = new double [N];
    double *TempX = new double[N];
    double a, norma=0, u=0, Tu=0;
	
	do
    {
		for (int i=0; i<N; i++){
			TempX[i]=x[i];//первое приближение
		}
		
		double Temp=0;
        for(int i=0; i<N; i++)
        {
			Temp=0;
			for(int j=0; j<N; j++){
				Temp+= A[i][j]*TempX[j];
			}
			r0[i]=Temp-b[i];//Вектор невязок
			Fi[i]=2*r0[i];
        }
		
		Temp=0;
        for(int i=0; i<N; i++)
        {
			for(int j=0; j<N; j++){
				Temp+= A[i][j]*r0[j];
			}
			Ar0[i]=Temp;
		}
		
        u=0.0; 
        Tu=0.0;
        for(int i=0; i<N; i++)
        {
            u+=abs(Ar0[i]*r0[i]);
            Tu+=abs(Ar0[i]*Ar0[i]);
        } 
		
        a=u/(2*Tu); 
        for(int i=0; i<N; i++){
			x[i]=TempX[i]-a*Fi[i];
		}
		
        norma=0;
        for(int i=0; i<N; i++)
        {
            norma+=(TempX[i]-x[i]);
        }
		
    }while (norma>eps);
}



/**
 * Метод сопряженных градиентов
 */
void kozlovaes::lab7()
{
  double* r = new double[N];
  double* r0 = new double[N];
  double* Az = new double[N];
  double* z = new double[N];
  double  eps = 1e-21, norma_b,norma_r,alpha,betta, norma;


  for(int i = 0; i<N; i++)
  norma_b += b[i]*b[i]; 

 double tmp=0;
  for (int i = 0; i<N; i++){
	r[i] = b[i];
    for(int j = 0; j<N; j++){ 
	  r[i]-=A[i][j]*x[j];
    }
    z[i] = r[i];
  }

  for (int i = 0; i<N; i++){
    tmp=0;
    for(int j = 0; j<N; j++){
      tmp +=A[i][j]*z[j]; 
    }
    Az[i] = tmp;
  }
 
  double tmp1 = 0, tmp2 = 0;
    for (int i = 0; i<N; i++){
      tmp1 += abs(r[i]*r[i]);
      tmp2 += abs(Az[i]*z[i]); 
    }
  alpha = tmp1/tmp2;

  //запоминаем предыдущие r[i]
  for(int i=0; i<N; i++){
    r0[i] = r[i];
  }

  for(int i=0; i<N; i++){
    x[i] = x[i] + alpha*z[i];
    r[i] = r[i] - alpha*Az[i];
  }

   for (int i = 0; i<N; i++){
      tmp1 += abs(r[i]*r[i]);
      tmp2 += abs(r0[i]*r0[i]); 
    }
  betta = tmp1/tmp2;

  for(int i = 0; i<N; i++){
    z[i] = r[i] + betta*z[i];
  }
  norma_r = 0;
  for(int i = 0; i<N; i++){
    norma_r +=r[i]*r[i];
  }

  norma = sqrt(norma_r/norma_b);

do{
  double tmp=0;
  for (int i = 0; i<N; i++){
	tmp=0;
    for(int j = 0; j<N; j++){
      tmp +=A[i][j]*z[j]; 
    }
    Az[i] = tmp;
  }

  double tmp1 = 0, tmp2 = 0;
    for (int i = 0; i<N; i++){
      tmp1 += abs(r[i]*r[i]);
      tmp2 += abs(Az[i]*z[i]); 
    }
  alpha = tmp1/tmp2;

  for(int i=0; i<N; i++){
    r0[i] = r[i];
  }

  for(int i=0; i<N; i++){
    x[i] = x[i] + alpha*z[i];
    r[i] = r[i] - alpha*Az[i];
  }

   for (int i = 0; i<N; i++){
      tmp1 += abs(r[i]*r[i]);
      tmp2 += abs(r0[i]*r0[i]); 
    }
  betta = tmp1/tmp2;

  for(int i = 0; i<N; i++){
    z[i] = r[i] + betta*z[i];
  }
  norma_r= 0;
  for(int i = 0; i<N; i++){
    norma_r +=r[i]*r[i];
  }

  norma = sqrt(norma_r/norma_b);

}while(norma>=eps);

}


void kozlovaes::lab8()
{
 double eps = 1.0e-6;
double errC = 0, c, s,fi;
int maxi, maxj;


double** C = new double*[N];;
  for (int i=0; i<N; i++){
  	C[i] = new double[N];
  }
  
   for (int i=0; i<N; i++){
    for (int j=i+1; j<N; j++){
	  C[i][j]=0;
	}
  }
  
  for (int i=0; i<N; i++){
	  for (int j=0; j<N; j++){
		  if (i!=j){
			  errC += A[i][j] * A[i][j];

		    }
	    }
    }
	
  while (sqrt(errC) > eps){

  	fi = 0; 
	maxi = 0, maxj = 1;
  	for (int i=0; i<N; i++){
      for (int j=i+1; j<N; j++){
  	    if (abs(A[i][j]) > abs(A[maxi][maxj])){
  	    	maxi = i;
  	    	maxj = j;
		  }
	  }
	}
	if (A[maxi][maxi] == A[maxj][maxj]) fi = M_PI/4;
	else fi = 0.5 *(atan(2*A[maxi][maxj] / (A[maxi][maxi]-A[maxj][maxj])));

	c = cos(fi), s = sin(fi);
	
	C[maxi][maxi]	= c*c*A[maxi][maxi] - 2*s*c*A[maxi][maxj] + s*s*A[maxj][maxj];
	C[maxj][maxj]	= s*s*A[maxi][maxi] + 2*s*c*A[maxi][maxj] + c*c*A[maxj][maxj];
	C[maxi][maxj] = (c*c - s*s)*A[maxi][maxj] + s*c*(A[maxj][maxj] - A[maxi][maxi]);
	C[maxj][maxi] = C[maxi][maxj];



	for (int k=0; k<N; k++){
	  if (k != maxi  &&  k != maxj){
	  	C[maxi][k] = c*A[maxi][k] - s*c*A[maxj][k];
	  	C[k][maxi] = C[maxi][k];
	  	C[maxj][k] = s*A[maxi][k] + c*A[maxj][k];
	  	C[k][maxj] = C[maxj][k];
	  } 
	  for (int l=0; l<N; l++){
	    if (k != maxi  &&  k != maxj  &&  l != maxi  &&  l != maxj)
	      C[k][l] = A[k][l];

	  }
	}

	errC = 0.;
	for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
  	    if (i != j){
		  errC += C[i][j] * C[i][j];
	    }

	  }
    }
	for (int i=0; i<N; i++)
      for (int j=0; j<N; j++)
	  A[i][j] = C[i][j];
  }

  for (int i=0; i<N; i++)
	  x[i] = A[i][i];
}


void kozlovaes::lab9()
{
double eps = 1e-3;
double* y0 = new double[N];	
double* y1 = new double[N]; 
double lambda1 = 0;
double lambda2 = 0;

for (int i = 0; i < N; i++){
		y0[i] = 1;
		y1[i] = 0;
	}
	
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			y1[i] += A[i][j] * y0[j];
		}
	}
	for (int i = 0; i < N; i++){
		if ((y1[i] != 0) && (y0[i]!=0)){
			lambda1 = y1[i]/y0[i];
			break;
		}
	}
	
	while(abs(lambda2-lambda1)>eps){
		
	lambda1 = lambda2;
	lambda2 = 0;
		
	for (int i = 0; i < N; i++){
		y0[i] = y1[i];
		y1[i] = 0;
	}
	
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			y1[i] += A[i][j] * y0[j];
		}
	}
	for (int i = 0; i < N; i++){
		if ((y1[i] != 0) && (y0[i]!=0)){
			lambda2 = y1[i]/y0[i];
			break;
		}
	}	
}
cout << "The largest in modulus eigenvalue : " << lambda2;
}


std::string kozlovaes::get_name()
{
  return "Kozlova E.S.";
}
