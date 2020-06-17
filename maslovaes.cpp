#include "maslovaes.h"

/**
 * Введение в дисциплину
 */
void maslovaes::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void maslovaes::lab2()
{
for (int k = 0; k < N; k++) {
	int idmax = 0;
	for (int i=0; i<N; i++) {
		if(abs(A[i][k]) > abs(A[idmax][k])) idmax = i; 
	}

	for (int i = 0; i < N; i++)
	{
		swap(A[idmax][i], A[k][i]);
	}
	swap(b[idmax], b[k]);

	//down
	for (int i = k + 1; i < N; i++)
	{
		double tmp = A[i][k]/A[k][k];
		for (int j = k; j < N; j++) {
			A[i][j] -= A[k][j]*tmp;
		}
		b[i] -= b[k]*tmp;
	}
}

for(int i = 0; i<N; i++){
    x[i]=b[i];
}

//up
for (int i = N-1; i >= 0; i--){
	for (int j = i+1; j < N; j++) {
		x[i] -= A[i][j]*x[j];
	}
	x[i] /= A[i][i];
}
}



/**
 * Метод прогонки
 */
void maslovaes::lab3()
{
double* alpha = new double[N];
double* beta = new double[N];
alpha[0]= - A[0][1]/A[0][0];
beta[0]=b[0]/A[0][0];
for(int i=1; i<N; i++){
    alpha[i] = - A[i][i + 1]/(A[i][i] + A[i][i - 1]*alpha[i - 1]);
    beta[i] = (- A[i][i - 1]*beta[i - 1] + b[i])/(A[i][i] + A[i][i - 1]*alpha[i - 1]);
}
x[N - 1] = beta[N - 1];
for(int i=N - 2; i>=0; i--){
    x[i]=alpha[i]*x[i + 1] + beta[i];
}
}



/**
 * Метод Холецкого
 */
void maslovaes::lab4()
{
double *D = new double[N];
double **S = new double* [N];  
for (int i=0; i<N; i++){
	S[i]=new double[N];
	for(int j=0; j<N; j++)
		S[i][j]=0;
}
if (A[0][0]>0) D[0] = 1;
else D[0] = -1;

S[0][0]=sqrt(fabs(A[0][0]));
for (int j = 1; j<N ;j++)
	S[0][j]=A[0][j]/(D[0]*S[0][0]);
	
for (int i=1; i<N; i++){
	double tmp =0;
    for (int j=0; j<i; j++){
		tmp+=D[j]*S[j][i]*S[j][i];
	}
	D[i] = copysign(1, A[i][i] - tmp);
	S[i][i]=sqrt(D[i]*(A[i][i]-tmp));
	  
	for (int j=i+1; j<N; j++) {
		double sum =0;
		for (int k=0; k<j; k++){
			sum+=D[k]*S[k][i]*S[k][j];
	    }	   
	S[i][j]=(A[i][j]-sum)/(D[i]*S[i][i]);
	}
}
double* y = new double[N];
y[0]=b[0]/S[0][0];
for (int i=1; i<N; i++){
	double tmp =0;
	for (int j=0; j<i; j++){
		tmp+=y[j]*S[j][i];
	}
		y[i]=(b[i]-tmp)/S[i][i];
}		
x[N-1]=y[N-1]/(D[N-1]*S[N-1][N-1]);
	 
for (int i=N-2; i>=0; i--){
		  double tmp =0;
	for (int j=i+1; j<N; j++){
		tmp+=x[j]*D[j]*S[i][j];
	}
	x[i]=(y[i]-tmp)/(D[i]*S[i][i]);
}
}



/**
 * Метод Якоби или Зейделя
 */
void maslovaes::lab5()
{
double e = 1e-30;
double *f = new double [N];
double tmp;
  do{
    for(int i = 0; i<N; i++)
      f[i]=x[i];

    for(int i = 0; i<N; i++){
      double sum1 = 0, sum2 = 0;
      for(int j = i + 1; j < N; j++)
      sum1 += A[i][j]*x[j];

      for(int j = i-1; j>= 0; j--)
      sum2 += A[i][j]*x[j];

      x[i] = (b[i] - sum1 - sum2)/A[i][i];
    }
	tmp = 0;
	  for (int i = 0; i < N; i++)
    {
	    tmp += abs (x[i] - f[i]);
    }
  } while (tmp>e);
}



/**
 * Метод минимальных невязок
 */
void maslovaes::lab6()
{
double* F = new double[N];
double* r = new double[N];
double norma, eps = 1e-15;

do {
	for (int i = 0; i < N; i++) {
		double tmp = 0;
		for (int j = 0; j < N; j++) {
			tmp += A[i][j] * x[j];
		}
	r[i] = tmp - b[i];
		F[i] = 2 * r[i];
	}
	double* A1 = new double[N];
	for (int i = 0; i < N; i++) {
		double temp = 0;
		for (int j = 0; j < N; j++) {
			temp += A[i][j] * r[j];
		}
		A1[i] = temp;
	}
	double t1 = 0, t2 = 0;
	for (int i = 0; i < N; i++) {
		t1 += abs(A1[i] * r[i]);
		t2 += abs(A1[i] * A1[i]);
	}
	double a = t1 / (2 * t2);

	double*y = new double[N];
	for (int i = 0; i < N; i++) {
		y[i] = x[i];
	}
	for (int i = 0; i < N; i++) {
		x[i] = x[i] - a * F[i];
	}

	norma = 0;
	for (int i = 0; i < N; i++) {
		norma += (y[i] - x[i])*(y[i] - x[i]);
	}

} while (sqrt(norma) > eps);
}



/**
 * Метод сопряженных градиентов
 */
void maslovaes::lab7()
{
 int i,j; double e = 1.0e-21;
 double *Az, *r, *z, *r1;
 
  z = new double[N];
  Az = new double[N];
  r = new double[N];
  r1 = new double[N];
  
  for (int i=0;i<N;i++){
  	x[i]=b[i];
  }
  
	for (i=0;i<N;i++){
  	  double s = 0;
  	  for (j=0;j<N;j++){
  		 s += A[i][j]*x[j];
	   }
	  r[i]=-s+b[i];
	  z[i]=r[i];
      }
  
  double norma = 1;
  double k=0;
  while (norma>e){
	double alfa1=0,alfa2=0;
	k++;

  	for (i=0;i<N;i++){
  		alfa1 += r[i]*r[i];
  	}
  	
  	for (i=0;i<N;i++){
  	  double s = 0;
  	  for (j=0;j<N;j++){
  		 s += A[i][j]*z[j];
	   }
	  Az[i]=s;
      }

    for (i=0;i<N;i++){
  		alfa2 += Az[i]*z[i];
  	}
   double alfa=alfa1/alfa2;
   
  
   for (i=0;i<N;i++){
  	  x[i]+=alfa*z[i];
    }
   
   for (i=0;i<N;i++){
  	  r1[i] = r[i] - alfa*Az[i];
    }
    double betta1=0;
  	for (i=0;i<N;i++){
  		betta1 += r1[i]*r1[i];
  	}
  	double betta=betta1/alfa1;
  	double u1=0,u2=0;
  	for (i=0;i<N;i++){
  	  z[i] = r1[i] + betta*z[i];
  	  r[i] = r1[i];
  	  u1+=r1[i]*r1[i];
  	  u2+=b[i]*b[i];
    }
    norma = sqrt(u1)/sqrt(u2);
  
  }
}


void maslovaes::lab8()
{
double eps = 1.0e-6;
    double errc = 0;
    double c = 0, s = 0, alpha = 0;
    int max_i = 0, max_j = 1;
    double** C;

    C = new double*[N];
 
		for (int i = 0; i < N; i++)
	{
	     C[i] = new double[N];
        for (int j = i+1; j < N; j++)
        {
          
            if (i != j){
                errc += A[i][j]*A[i][j];
            }
            	C[i][j] = 0;
        }
    }
    while (sqrt(errc) > eps)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = i + 1; j < N; j++)
            {
                if (abs(A[i][j]) >= abs(A[max_i][max_j]))
                {

                    max_i = i;
                    max_j = j;
                }
            }
        }
        if (A[max_i][max_i] == A[max_j][max_j])
            alpha = M_PI / 4;
		else alpha = 0.5*atan((2 * A[max_i][max_j]) / (A[max_j][max_j] - A[max_i][max_i]));
        c = cos(alpha);
        s = sin(alpha);

        C[max_i][max_i]= c*c*A[max_i][max_i] - 2*s*c*A[max_i][max_j]+s*s*A[max_j][max_j];
        C[max_j][max_j]= s*s*A[max_i][max_i]+2*s*c*A[max_i][max_j]+c*c*A[max_j][max_j];
        C[max_i][max_j]= (c*c-s*s)*A[max_i][max_j]+s*c*(A[max_i][max_i]-A[max_j][max_j]);
        C[max_j][max_i]= C[max_i][max_j];

        for (int k = 0; k < N; k++)
        {
            if (k != max_i && k != max_j)
            {
                C[max_i][k] = c * A[max_i][k] - s * c * A[max_j][k];
                C[k][max_i] = C[max_i][k];
                C[max_j][k] = s * A[max_i][k] + c * A[max_j][k];
                C[k][max_j] = C[max_j][k];
            }
            for (int l = 0; l < N; l++)
            {
                if (k != max_i && k != max_j && l != max_i && l != max_j)
                    C[k][l] = A[k][l];
            }
        }
    
    errc = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            if (i != j)
            {
                errc += C[i][j] * C[i][j];
            }
        }
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i][j] = C[i][j];
        }
    }

    for (int i = 0; i < N; i++)
        x[i] = A[i][i];
}
}


void maslovaes::lab9()
{
double eps = 1.0e-3;
  double *yk, *yk_1;
  double lambda_k1 = 0, lambda_k = 1;
  int n=N;
 
  yk_1 = new double[n];
  yk = new double[n];

  for (int i=0; i<n; i++)
  	yk_1[i]=1;
  
  while(fabs(lambda_k - lambda_k1) > eps){
  	lambda_k1 = lambda_k;
  	
	for(int i=0; i<n; i++){
	  double s = 0;
      for(int j=0; j<n; j++){
        s += A[i][j]*yk_1[j];
        }
	  yk[i]=s;
    }
    
   for (int i=0; i<n; i++){
   	if (yk[0] != 0 && yk_1[0] != 0){
   		lambda_k = yk[i]/yk_1[i];
	   }
   }
   
   for (int i=0; i<n; i++)
   	yk_1[i]=yk[i];
   	
  }
  cout<<"Answer: " << lambda_k << endl;
}


std::string maslovaes::get_name()
{
  return "Maslova E.S.";
}
