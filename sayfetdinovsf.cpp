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
double eps = 1e-15;
double t = 1e-5;
for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double x1;
	double *xr = new double[N];
	int step = 0;

	do {
		step++;
		for (int i = 0; i < N; i++) {
			xr[i] = x[i];
			for (int k = 0; k < N; k++)
				xr[i] -= t*A[i][k] * x[k];
			xr[i] += t * b[i];

		}
		x1 = 0.;
		for (int i = 0; i < N; i++) {
			x1 += (xr[i]-x[i])*(xr[i]-x[i]);
		}

		for (int i = 0; i < N; i++) {
			x[i] = xr[i];
		}
	} while (sqrt(x1)>eps);
}



/**
 * Метод Якоби или Зейделя
 */
void sayfetdinovsf::lab5()
{
//Метод Якоби
	double *oldx = new double[N]; 
	for (int i=0; i<N; i++) { 
		x[i]=0;} // заполняем решение нулями 
	double Err=0.0; 
	double eps=1e-20; // погрешность
	int k=0; 
	do { 
		k++; 
		Err=0.0; 
		for(int i=0; i<N; i++) 
			oldx[i]=x[i]; // предыдущее решение 
			for(int i=0; i<N; i++) 
			{ 
				double s=0;
				for(int j=0; j<i; j++) 
					s += A[i][j] * oldx[j]; //под главной диагональю
				for(int j=i+1; j<N; j++) 
					s += A[i][j] * oldx[j]; //над главной диагональю
				x[i]=(b[i] - s)/A[i][i]; // вычисляется новое решение 
			}			 
			Err=std::abs(oldx[0]-x[0]); 
			for(int i=0; i<N; i++) 
			{ 
				if(std::abs(oldx[i]-x[i]) > Err) 
				Err = std::abs(oldx[i]-x[i]);//максимальная разница между предыдущим решением и текущим. 
			} 
	} while(Err >= eps); 
delete [] oldx; 
}



/**
 * Метод минимальных невязок
 */
void sayfetdinovsf::lab6()
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
void sayfetdinovsf::lab7()
{
	double eps = 1e-20;
	double* prevX = new double[N];
	double* prevR = new double[N];
	double* r = new double[N];
	double* z = new double[N];
	for (int i = 0; i < N; i++) {
		r[i] = b[i];
		z[i] = r[i];
	}

	while (true) {

		for (int i = 0; i < N; i++) {
			prevR[i] = r[i];
			prevX[i] = x[i];
		}

		double alpha = 0, denAlpha = 0;

		for (int i = 0; i < N; i++) {
			double Az = 0;
			for (int j = 0; j < N; j++) {
				Az += A[i][j] * z[j];
			}
			alpha += prevR[i] * prevR[i];
			denAlpha += Az * z[i];
		}
		alpha /= denAlpha;

		for (int i = 0; i < N; i++) {
			x[i] = prevX[i] + alpha * z[i];
		}

		double maxErr = abs(x[0] - prevX[0]);
		for (int i = 1; i < N; i++)
			if (abs(x[i] - prevX[i]) > maxErr)
				maxErr = abs(x[i] - prevX[i]);

		if (maxErr < eps)
			break;

		for (int i = 0; i < N; i++) {
			double Az = 0;

			for (int j = 0; j < N; j++) {
				Az += A[i][j] * z[j];
			}

			r[i] = prevR[i] - alpha * Az;
		}

		double beta = 0, denBeta = 0;
		for (int i = 0; i < N; i++) {
			beta += r[i] * r[i];
			denBeta += prevR[i] * prevR[i];
		}
		beta /= denBeta;

		for (int i = 0; i < N; i++) {
			z[i] = r[i] + beta * z[i];
		}
	}

	delete[] prevX;
	delete[] r;
	delete[] prevR;
	delete[] z;
}


void sayfetdinovsf::lab8()
{
double eps = 1e-20;
	double **H = new double*[N];
	for (int i = 0; i < N; i++) 
		H[i] = new double[N];

	
	while(true){
		double n = 0;
		int i_max = 0, j_max = 1;
		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
			{
				if (abs(A[i][j]) > abs(A[i_max][j_max]))
				{
					i_max = i;
					j_max = j;
				}

				n += A[i][j] * A[i][j];
			}

		if (sqrt(n) < eps) break;

		double fi = 0.5 * atan(2 * A[i_max][j_max] / (A[i_max][i_max] - A[j_max][j_max]));
		for (int i = 0; i < N; i++)
		{
			H[i][i_max] = A[i][i_max] * cos(fi) + A[i][j_max] * sin(fi);
			H[i][j_max] = A[i][j_max] * cos(fi) - A[i][i_max] * sin(fi);
		}

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				if (j != i_max && j != j_max) H[i][j] = A[i][j];

		for (int j = 0; j < N; j++)
		{
			A[i_max][j] = H[i_max][j] * cos(fi) + H[j_max][j] * sin(fi);
			A[j_max][j] = H[j_max][j] * cos(fi) - H[i_max][j] * sin(fi);
		}

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				if (i != i_max && i != j_max) A[i][j] = H[i][j];

	}

	for (int i = 0; i < N; i++) 
		x[i] = A[i][i];

	for (int i = 0; i < N; i++) 
		delete[] H[i];
	delete[] H;
}


void sayfetdinovsf::lab9()
{
double eps = 1e-20;
	double lambda = 0;
	double plambda = 0;

	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	x[0] = 1;

	double *y = new double[N];

	while (true) {

		for (int i = 0; i < N; i++) {
			y[i] = 0;
			for (int j = 0; j < N; j++) {
				y[i] += (A[i][j] * x[j]);
			}
		}

		lambda = 0;
		for (int i = 0; i < N; i++) {
			lambda += (x[i] * y[i]);
		}

		if (abs(plambda - lambda) < eps) {
			break;
		}

		plambda = lambda;


		double norma_y = 0;
		for (int i = 0; i < N; i++) {
			norma_y += (y[i] * y[i]);
		}
		norma_y = sqrt(norma_y);

		for (int i = 0; i < N; i++) {
			x[i] = y[i] / norma_y;
		}
	}


	cout << lambda << endl;

	delete[] y;
}


std::string sayfetdinovsf::get_name()
{
  return "sayfetdinovsf";
}
