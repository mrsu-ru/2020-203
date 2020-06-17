#include "dvoryaninovada.h"

/**
 * Введение в дисциплину
 */
void dvoryaninovada ::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void dvoryaninovada::lab2()
{
	double maxi;
	int k, u;
	double eps=0.000001;

	k=0;
	while (k<N)
	{
		maxi=abs(A[k][k]);
		u=k;
		for (int i=k+1;i<N;i++) {
			if (abs(A[i][k])>maxi) {
				maxi=abs(A[i][k]);
				u=i;
			}
		}

		for (int j=0;j<N;j++) {
			swap(A[k][j], A[u][j]);
		}

		swap(b[k], b[u]);
		for (int i=k;i<N;i++) {
			double c=A[i][k];
			if (abs(c)<eps) continue;
			for (int j=0;j<N;j++) {
				A[i][j]=A[i][j]/c;
			}

			b[i]=b[i]/c;
			if (i==k)	continue;
			for (int j=0;j<N;j++) {
				A[i][j]=A[i][j]-A[k][j];
			}

			b[i]=b[i]-b[k];
		}

		k++;
	}

	for (k=N-1;k>=0;k--) {
		x[k]=b[k];
		for (int i=k-1; i>=0;i--) {
			double c= A[i][k];
			for (int j=0;j<N; j++) { 
				A[i][j]=A[k][j]*c+A[i][j]; 
			}
			b[i]=-b[k]*c+b[i];
		}
	}
}



/**
 * Метод прогонки
 */
void dvoryaninovada::lab3()
{
	double *P=new double[N]; 
	double *Q=new double[N];

	for (int i=0;i<N;i++) {
		P[i]=0;
		Q[i]=0;
	}

	P[0]=(-A[0][1])/A[0][0];
	Q[0]=b[0]/A[0][0];

	for (int i=1;i<N;i++) {
		P[i]=A[i][i+1]/(-A[i][i]-A[i][i-1]*P[i-1]);
		Q[i]=(-b[i]+A[i][i-1]*Q[i-1])/(-A[i][i]-A[i][i-1]*P[i-1]);
	}

	x[N-1]=Q[N-1];
	for (int i=N-2;i>=0;i--) {
		x[i]=P[i]*x[i+1]+Q[i];
	}

	delete[] P;
	delete[] Q;
}



/**
 * Метод простых итераций
 */
void dvoryaninovada::lab4()
{
	double **S=new double *[N];
	for (int i=0;i<N;i++)
		S[i]=new double[N];

	double *y=new double[N];


	for (int i=0;i<N;i++){
		x[i]=0;
		y[i]=0;
		for (int j=0;j<N;j++){
			S[i][j]=0;
		}
	}

	double t=0;
	for (int i=0;i<N;i++) {
		for (int k=0;k<=i-1;k++) t+=S[i][k]*S[i][k];

		S[i][i]=sqrt(A[i][i]-t);
		t=0;
		for (int j=i+1;j<N;j++)	{
			for (int k=0;k<=i-1;k++) t+=S[i][k]*S[j][k];

			S[j][i]=(A[i][j]-t)/S[i][i];
			t=0;
		}
	}


	for (int i=0;i<N;i++) {
		t=0;
		for (int j=0;j<i;j++) t+=S[i][j]*y[j];

		y[i]=(b[i]-t)/S[i][i];
	}


	for (int i=N-1;i>=0;i--) {
		t=0;
		for (int j=i+1;j<N;j++) t+=S[j][i]*x[j];

		x[i]=(y[i]-t)/S[i][i];
	}

	delete[] y;
	for (int i=0;i<N;i++) delete[] S[i];
	delete[] S;
}



/**
 * Метод Якоби или Зейделя
 */
void dvoryaninovada::lab5()
{
	double eps=0.000001;

	double* y=new double[N];
	double r=0; 
	double var=0;

	for (int i=0;i<N;i++){
		x[i]=0;
	}

	do{
		for (int i=0;i<N;i++)	y[i]=x[i];

		for (int i=0;i<N;i++){
			var=0; r=0;
			for (int j=0;j<i;j++) var+=A[i][j]*x[j];
			for (int j=i+1;j<N;j++)	var+=A[i][j]*x[j];

			x[i]=(b[i]-var)/A[i][i];
			for (int i=0;i<N;i++) r+=sqrt((x[i]-y[i])*(x[i]-y[i]));
		}
	} while (r>=eps);
	delete[] y;
}


void dvoryaninovada::MatrVect(double **M, double *V, double *R)
{
	for (int i = 0; i < N; i++)
	{
		R[i] = 0;
		for (int j = 0; j < N; j++)
			R[i] += M[i][j] * V[j];
	}
}

double dvoryaninovada::ScalarVect(double* v1, double* v2)
{
	double result = 0;
	for (int i = 0; i < N; i++)
		result += (v1[i] * v2[i]);
	return result;
}
/**
 * Метод минимальных невязок
 */
void dvoryaninovada::lab6()
{
	double eps=0.000001;

	int count=0; 
	double *U=new double[N];
	double *r=new double[N];
	double *TempX=new double[N];
	double *p=new double[N];
	double Tau=0.0;

	for (int i=0;i<N;i++) TempX[i]=0; 

	do	{
		MatrVect(A,TempX,U);
		for (int i=0;i<N;i++)	{
			r[i]=U[i]-b[i]; 
		}

		MatrVect(A,r,U);

		double TempTau1=ScalarVect(U,r);
		double TempTau2=ScalarVect(U,U);
		if (TempTau2==0) break;

		Tau=TempTau1/TempTau2;

		for (int i=0;i<N;i++) x[i]=TempX[i]-Tau*r[i];

		for (int i=0;i<N;i++) p[i]=x[i]-TempX[i];

		count++;
	} while ((sqrt(ScalarVect(p,p))>=eps)&&(count<500000));

	delete[] U;
	delete[] r;
	delete[] p;
	delete[] TempX;
}



/**
 * Метод сопряженных градиентов
 */
void dvoryaninovada::lab7()
{
	double Eps=0.000005;
	double Del,s,sAbs;//погрешность итерации, скалярный шаг, модуль шага


	double *K=new double[N];
	double *L=new double[N];
	double *M=new double[N];
	double *xrez=new double[N];


	 
	for (int i=0;i<N;i++) xrez[i]=0; //начальное приближение


	do {
		
		for (int i=0;i<N;i++) {
			K[i]=0;
			for (int j=0;j<N;j++) K[i]+=A[i][j]*xrez[j]; //скалярное произведении матрицы системы и вектор приближенного решения
		}

		
		for (int i=0;i<N;i++) L[i]=K[i]-b[i];//градиент

		
		for (int i=0;i<N;i++) {
			K[i]=0;
			for (int j=0;j<N;j++) K[i]+=A[i][j]*L[j]; //скалярное произведение матрицы системы и градиента
		}


		for (int i=0;i<N;i++) {
			M[i]=0;
			for (int j=0;j<N;j++) M[i]+=A[i][j]*K[j];
		}

		s=0;
		sAbs=0;

		for (int i=0;i<N;i++) {
			s+=K[i]*L[i];
			sAbs+=M[i]*K[i]; //величина смещения по направлению градиент
		}

		if (s==sAbs) s=1;
		else
			s=s/sAbs;
		
		for (int i=0;i<N;i++) x[i]=xrez[i]-s*L[i];// новое приближенное решение

		//проверка на уменьшение погрешности
		Del=abs(x[0]-xrez[0]);

		for (int i=0;i<N;i++) {
			if (abs(x[i]-xrez[i])>Del)	Del=abs(x[i]-xrez[i]);
			xrez[i]=x[i];
		}
	} while (Eps<Del);
}


void dvoryaninovada::lab8()
{
	double **B = new double*[N];
	for (int i = 0; i < N; i++)
		B[i] = new double[N];
	double eps = 1.e-10;
	int i_var, j_var;
	for (;;) {
		i_var = 0;
		j_var = 1;
		double n = 0;
		for (int i = 0; i < N - 1; i++)
			for (int j = i + 1; j < N; j++) {
				if (abs(A[i][j]) > abs(A[i_var][j_var])) {
					i_var = i;
					j_var = j;
				}
				n += A[i][j] * A[i][j];
			}
		if (sqrt(n) < eps) break;
		double fi = 0.5*atan(2 * A[i_var][j_var] / (A[i_var][i_var] - A[j_var][j_var]));
		for (int i = 0; i < N; i++) {
			B[i][i_var] = A[i][i_var] * cos(fi) + A[i][j_var] * sin(fi);
			B[i][j_var] = -A[i][i_var] * sin(fi) + A[i][j_var] * cos(fi);
		}
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				if (j != j_var && j != i_var) B[i][j] = A[i][j];
		for (int i = 0; i < N; i++) {
			A[i_var][i] = B[i_var][i] * cos(fi) + B[j_var][i] * sin(fi);
			A[j_var][i] = -B[i_var][i] * sin(fi) + B[j_var][i] * cos(fi);
		}
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				if (i != j_var && i != i_var) A[i][j] = B[i][j];
	}
	for (int i = 0; i < N; i++)
		x[i] = A[i][i];

	for (int i = 0; i < N; i++)
		delete[] B[i];
	delete[] B;
}


void dvoryaninovada::lab9()
{

	double * Y = new double[N];//первый вектор приближения
	double * M = new double[N];//второй вектор приближения
	double maxL, L, sum;
	double EPS = 1e-15;

	for (int i = 0; i < N; i++) Y[i] = 0;
	Y[0] = 1;

	do {
		sum = 0;
		for (int i = 0; i < N; i++)	sum += Y[i] * Y[i];

		L = sqrt(sum);

		//построение последовательности векторов
		for (int i = 0; i < N; i++)	{
			M[i] = 0;
			for (int j = 0; j < N; j++)	M[i] += A[i][j] * Y[j] / L;
		}
		sum = 0;

		//сравнение нормы полученного вектора с заданной погрешностью
		for (int i = 0; i < N; i++)	sum += M[i] * M[i];
		maxL = sqrt(sum);

		for (int i = 0; i < N; i++)	Y[i] = M[i];
	} while (abs(maxL - L) > EPS);

	cout << maxL << endl;
}


std::string dvoryaninovada::get_name()
{
  return "Dvoryaninova D. A.";
}
