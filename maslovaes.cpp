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

}


void maslovaes::lab8()
{

}


void maslovaes::lab9()
{

}


std::string maslovaes::get_name()
{
  return "Maslova E.S.";
}
