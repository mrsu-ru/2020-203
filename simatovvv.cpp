#include "simatovvv.h"

/**
 * Введение в дисциплину
 */
void simatovvv::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void simatovvv::lab2()
{
{
for (int i = 0; i < N; i++)
x[i] = b[i];
long double m;
for (int k = 0; k < N - 1; k++)
{
for (int i = k + 1; i < N; i++)
{
m = A[i][k] / A[k][k];
for (int j = k; j < N; j++)
{
A[i][j] = A[i][j] - m * A[k][j];
}
x[i] = x[i] - m * x[k];
}
}
for (int i = N - 1; i >= 0; i--)
{
for (int j = i + 1; j < N; j++)
x[i] = x[i] - A[i][j] * x[j];
x[i] = x[i] / A[i][i];
}
}
}



/**
 * Метод прогонки
 */
void simatovvv::lab3()
{
double *alpha = new double[N];
double *beta = new double[N];
alpha[0] = -A[0][1] / A[0][0];
beta[0] = b[0] / A[0][0];
for (int i = 1; i <= N-2; i++)
{
alpha[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * alpha[i - 1]);
beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / (A[i][i] + A[i][i - 1] * alpha[i - 1]);

}
beta[N - 1] = (b[N - 1] - A[N - 1][N - 2] * beta[N - 2]) / (A[N - 1][N - 1] + A[N - 1][N - 2] * alpha[N - 2]);
x[N - 1] = beta[N - 1];
for (int i = N-2; i >= 0; i--)
{
x[i] = alpha[i] * x[i + 1] + beta[i];
}
}



/**
 * Метод Холецкого
 */
void simatovvv::lab4()
{
double** L = new double*[N];
	for (int i=0; i<N; i++)
	L[i] = new double[N];
	double* y = new double[N];

   for (int i=0; i<N; i++)
   {
		for (int j=0; j<N; j++)
			{
				L[i][j]=0;
			}
   }

  for(int i=0;i<N;i++)
	{
		for(int k=0;k<=i-1;k++)

			L[i][i]+=L[i][k]*L[i][k];
			L[i][i]=sqrt(A[i][i]-L[i][i]);
			for(int j=i+1;j<N;j++)
				{
					for(int k=0;k<=i-1;k++)
						L[j][i]+=L[i][k]*L[j][k];
						L[j][i]=(A[i][j]-L[j][i])/L[i][i];
				}
  }
	double summa=0;
    for(int i=0;i<N;i++)
		{
			for(int j=0;j<i;j++){
			summa+=L[i][j]*y[j];}
			y[i]=(b[i]-summa)/L[i][i];
			summa=0;
		}
	for(int i=N-1;i>=0;i--)
		{
			for(int j=i+1;j<N;j++){
				summa+=L[j][i]*x[j];}
				x[i]=(y[i]-summa)/L[i][i];
				summa=0;
		}
delete[] y;
}



/**
 * Метод Якоби или Зейделя
 */
void simatovvv::lab5()
{
long double eps = 0.00000001;
    long double* p = new long double[N];
	long double norm;
    for (int i = 0; i < N; i++)
        x[i] = 0;
    do {
		for (int i = 0; i < N; i++)
        {
			p[i] = b[i];
			for (int j = 0; j < N; j++)
				{
				    if (i != j)
                     p[i] -= A[i][j] * x[j];
				}
			p[i] /= A[i][i];
		}
        norm = fabs(x[0] - p[0]);
		for (int h = 0; h < N; h++)
        {
			if (fabs(x[h] - p[h]) > norm)
				norm = fabs(x[h] - p[h]);
			x[h] = p[h];
		}
	} while (norm > eps);
    delete[] p;
}



/**
 * Метод минимальных невязок
 */
void simatovvv::lab6()
{
double *prevX = new double[N];
	double *y = new double[N];
	double tau = 0.0, err = 0.0, Ay = 0.0, denom = 0.0;

	do{

		for (int i = 0; i < N; i++) {
			y[i] = b[i];
			for (int j = 0; j < N; j++) y[i] -= A[i][j] * prevX[j]; 
		}

		tau = 0.0; denom = 0.0;

		for (int i = 0; i < N; i++) {
			Ay = 0.0;

			for (int j = 0; j < N; j++) Ay += A[i][j] * y[j];

			tau += Ay * y[i]; denom += Ay * Ay;
		}
		tau /= denom; 

		for (int i = 0; i < N; i++) x[i] = prevX[i] + (tau * y[i]);
		
		err = 0.0;
		for (int i = 0; i < N; i++) if (abs(x[i] - prevX[i]) > err) err = abs(x[i] - prevX[i]);
			
		for (int i = 0; i < N; i++) prevX[i] = x[i];
	}
	while(err > 1e-20);
}



/**
 * Метод сопряженных градиентов
 */
void simatovvv::lab7()
{
        double* prevX = new double[N];
	double* prevR = new double[N];
	double* r = new double[N];
	double* z = new double[N];

	for (int i = 0; i < N; i++) {
		r[i] = b[i];
		z[i] = b[i];
	}

	double err = 0;
	double a = 0;
	double det = 0;
	double beta = 0;
	double Az = 0;

	do {

		for (int i = 0; i < N; i++) {
			prevR[i] = r[i];
			prevX[i] = x[i];
		}

		a = 0; det = 0;

		for (int i = 0; i < N; i++) {
			Az = 0;
			
			for (int j = 0; j < N; j++) Az += (A[i][j] * z[j]);
			
			a += (prevR[i] * prevR[i]); det += (Az * z[i]);
		}
		a /= det;

		for (int i = 0; i < N; i++)	x[i] = prevX[i] + (a * z[i]);
		

		err = abs(x[0] - prevX[0]);

		for (int i = 1; i < N; i++)
			if (abs(x[i] - prevX[i]) > err)
				err = abs(x[i] - prevX[i]);

		for (int i = 0; i < N; i++) {
			Az = 0;

			for (int j = 0; j < N; j++) Az += (A[i][j] * z[j]);

			r[i] = prevR[i] - (a * Az);
		}

		beta = 0; det = 0;

		for (int i = 0; i < N; i++) {
			beta += (r[i] * r[i]);
			det += (prevR[i] * prevR[i]);
		}
		beta /= det;

		for (int i = 0; i < N; i++) z[i] = r[i] + (beta * z[i]);

	}while( !(err < 1e-20) );
}


void simatovvv::lab8()
{
double eps=1.e-20;
	double B[N][N], norm;
	int imax, jmax;
	for(;;){
		imax=0; jmax=1;
		norm=0;
		for (int i=0; i<N-1; i++) {
			for (int j=i+1; j<N; j++) {
				if (abs(A[i][j])>abs(A[imax][jmax])) {
					imax=i;
					jmax=j;
				}
				norm+=A[i][j]*A[i][j];
			}
		}

		if (sqrt(norm)<eps) {
			break;
		}

        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                B[i][j]=A[i][j];
            }
        }

		double fi=0.5*atan(2*A[imax][jmax]/(A[imax][imax]-A[jmax][jmax]));
		for (int i=0; i<N; i++) {
			B[i][imax]=A[i][imax]*cos(fi)+A[i][jmax]*sin(fi);
			B[i][jmax]=-A[i][imax]*sin(fi)+A[i][jmax]*cos(fi);
		}

		for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                A[i][j]=B[i][j];
            }
        }

		for (int i=0; i<N; i++) {
			A[imax][i]=B[imax][i]*cos(fi)+B[jmax][i]*sin(fi);
			A[jmax][i]=-B[imax][i]*sin(fi)+B[jmax][i]*cos(fi);
		}
	}

	for (int i = 0; i < N; i++) {
		x[i]=A[i][i];
	}

}


void simatovvv::lab9()
{
double eps = 1e-3;
	double* yPrev = new double[N];	
	double* yNext = new double[N]; 
	double y0 = 0;
	double y1 = 0;
	double lambdaPrev = 0;
	double lambdaNext = 0;
	double delta = 0;

	for (int i = 0; i < N; i++){
		yPrev[i] = 1;
		yNext[i] = 0;
	}

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			yNext[i] += A[i][j] * yPrev[j];
		}
	}

	y0 = yPrev[0];

	for (int i = 0; i < N; i++){
		if (yNext[i] != 0){
			y1 = yNext[i];
			break;
		}
	}

	lambdaPrev = y1 / y0;
	delta = lambdaPrev;

	while (delta > eps){
		for (int i = 0; i < N; i++){
			yPrev[i] = yNext[i];
			yNext[i] = 0;
		}

		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				yNext[i] += A[i][j] * yPrev[j];
			}
		}

		for (int i = 0; i < N; i++){
			if ((yNext[i] != 0) && (yPrev[i] != 0)){
				y0 = yPrev[i];
				y1 = yNext[i];
				break;
			}
		}

		lambdaNext = y1 / y0;
		delta = fabs(lambdaNext - lambdaPrev);
		lambdaPrev = lambdaNext;
	}
	cout << "Max self value: " << lambdaNext << endl;
	delete[]yPrev;
	delete[]yNext;

}


std::string simatovvv::get_name()
{
  return "simatovvv";
}
