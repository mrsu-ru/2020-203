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
	for (int i = 0; i < N; i++) {
		int max = i;
		for (int j = i + 1; j < N; j++) {
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

		for (int j = i + 1; j < N; j++) {

			for (int k = N; k > i; k--) A[j][k] -= A[i][k] * A[j][i];
			b[j] -= b[i] * A[j][i];
			A[j][i] = 0;
		}
	}
	double s;
	for (int i = N - 1; i > 0; i--) {
		s = 0;
		for (int j = i ; j < N; j++) {
			s += b[j] * A[i - 1][j];
			A[i - 1][j] = 0;
		}
		b[i - 1] -=s;
	}
	for (int i = 0; i < N; i++) {
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
 

	for(i = 1;i < N;i++){
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
 * Метод квадратного корня (метод Холецкого)
 */
void kozlovaes::lab4()
{
double* d=new double[N];
double** S = new double*[N];
for (int i = 0; i < N; i++) { 
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
for(int i=1; i < N; i++){
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
for(int j=i+1;j < N;j++){
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

}



/**
 * Метод сопряженных градиентов
 */
void kozlovaes::lab7()
{

}


void kozlovaes::lab8()
{

}


void kozlovaes::lab9()
{

}


std::string kozlovaes::get_name()
{
  return "Kozlova E.S.";
}
