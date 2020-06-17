#include "ashryatovarr.h"

/**
 * Введение в дисциплину
 */
void ashryatovarr::lab1()
{
    cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void ashryatovarr::lab2()
{

    for (int k = 0; k < N; k++)
    {
        int idmax = 0;
        for (int i=0; i<N; i++)
        {
            if(abs(A[i][k]) > abs(A[idmax][k]))
                idmax = i;
        }

        for (int i = 0; i < N; i++)
        {
            swap(A[idmax][i], A[k][i]);
        }
        swap(b[idmax], b[k]);

        for (int i = k + 1; i < N; i++)
        {
            double tmp = A[i][k]/A[k][k];
            for (int j = k; j < N; j++)
            {
                A[i][j] -= A[k][j]*tmp;
            }
            b[i] -= b[k]*tmp;
        }
    }

    for(int i = 0; i<N; i++)
    {
        x[i]=b[i];
    }

    for (int i = N-1; i >= 0; i--)
    {
        for (int j = i+1; j < N; j++)
        {
            x[i] -= A[i][j]*x[j];
        }
        x[i] /= A[i][i];
    }

}



/**
 * Метод прогонки
 */
void ashryatovarr::lab3()
{
    double *alfa = new double[N];
    double *beta = new double[N];
    alfa[0]=A[0][1]/-A[0][0];
	beta[0]=b[0]/A[0][0];
	for( int i=1; i<N; i++){
		alfa[i]=A[i][i+1]/(-A[i][i]-A[i][i-1]*alfa[i-1]);
		beta[i]=(A[i][i-1]*beta[i-1]-b[i])/(-A[i][i]-A[i][i-1]*alfa[i-1]);
}

  x[N-1]=beta[N-1];
  for (int i=N-2; i>=0; i--){
	  x[i]=alfa[i]*x[i+1]+beta[i];
  }

}



/**
 * Метод Холецкого
 */
void ashryatovarr::lab4()
{
    double **S = new double*[N];
    for (int i=0; i<N; i++)
    {
        S[i]=new double[N];
        for(int j=0; j<N; j++)
            S[i][j]=0;
    }
    double *D = new double[N];
    if (A[0][0]>0)
        D[0]=1;
    else
        D[0]=-1;
    S[0][0]=sqrt(fabs(A[0][0]));

    for (int i=1; i<N; i++)
    {
        S[0][i]=A[0][i]/(D[0]*S[0][0]);
    }

    for (int i=1; i<N; i++)
    {
        double temp =0;
        for (int j=0; j<i; j++)
            temp+=D[j]*S[j][i]*S[j][i];
        if (A[i][i]-temp>=0)
            D[i]=1;
        else
            D[i]=-1;
        S[i][i]=sqrt(D[i]*(A[i][i]-temp));

        for (int j=i+1; j<N; j++)
        {
            double l = 0;
            for (int k=0; k<j; k++)
                l+=D[k]*S[k][i]*S[k][j];

            S[i][j]=(A[i][j]-l)/(D[i]*S[i][i]);
        }
    }
    double *y = new double[N];
    y[0]=b[0]/S[0][0];
    for (int i=1; i<N; i++)
    {
        double temp = 0;
        for (int j=0; j<i; j++)
            temp+=y[j]*S[j][i];
        y[i]=(b[i]-temp)/S[i][i];
    }
    x[N-1]=y[N-1]/(D[N-1]*S[N-1][N-1]);

    for (int i=N-2; i>=0; i--)
    {
        double temp =0;
        for (int j=i+1; j<N; j++)
            temp+=x[j]*D[j]*S[i][j];      
        x[i]=(y[i]-temp)/(D[i]*S[i][i]);
    }

}



/**
 * Метод Якоби или Зейделя
 */
void ashryatovarr::lab5()
{
    int n = N;
    double eps = 1e-69;
    double norma;
    double *y = new double[n];
    do{
        for(int i = 0; i<n; i++)
            y[i]=x[i];
        norma = 0;
        for(int i = 0; i<n; i++)
        {
            double sum1 = 0, sum2 = 0;
            for(int j = i + 1; j < n; j++)
                sum1 += A[i][j]*x[j];
            for(int j = i-1; j>= 0; j--)
                sum2 += A[i][j]*x[j];
            x[i] = (b[i] - sum1 - sum2)/A[i][i];
        }

        for (int i = 0; i < n; i++)
            norma += abs(x[i] - y[i]);
    } while (norma>eps);

}



/**
 * Метод минимальных невязок
 */
void ashryatovarr::lab6()
{
	int n = N;
	double* F = new double[n];
	double* r = new double[n];
	double* alfa = new double[n];
	double a, norma, k = 0;
	double eps = 1e-19;
	do {
		for (int i = 0; i < n; i++) {
			double tmp = 0;
			for (int j = 0; j < n; j++)
				tmp += A[i][j] * x[j];
			r[i] = tmp - b[i];
			F[i] = 2 * r[i];
		}
		double* Ar = new double[n];
		for (int i = 0; i < n; i++) {
			double tmps = 0;
			for (int j = 0; j < n; j++)
				tmps += A[i][j] * r[j];
			Ar[i] = tmps;
		}
		double ts1 = 0, ts2 = 0;
		for (int i = 0; i < n; i++) {
			ts1 += abs(Ar[i] * r[i]);
			ts2 += abs(Ar[i] * Ar[i]);
		}
		a = ts1 / (2 * ts2);

		double*y = new double[n];
		for (int i = 0; i < n; i++)
			y[i] = x[i];
		for (int i = 0; i < n; i++)
			x[i] = x[i] - a * F[i];

		norma = 0;
		for (int i = 0; i < n; i++)
			norma += (y[i] - x[i])*(y[i] - x[i]);

	}while(sqrt(norma)>eps);

}



/**
 * Метод сопряженных градиентов
 */
void ashryatovarr::lab7()
{
		
	double *xrez = new double[N];

    for (int i = 0; i<N; i++)
        xrez[i] = 0;

    double Del, s, sAbs;
    double eps = 1.e-10;
    double *K = new double[N];
    double *L = new double[N];
    double *M = new double[N];

    do
    {
        for (int i = 0; i < N; i++)
        {
            K[i] = 0;
            for (int j = 0; j < N; j++)
                K[i] += A[i][j] * xrez[j];
        }

        for (int i = 0; i < N; i++)
            L[i] = K[i] - b[i];


        for (int i = 0; i < N; i++)
        {
            K[i] = 0;
            for (int j = 0; j < N; j++)
                K[i] += A[i][j] * L[j];
        }


        for (int i = 0; i < N; i++)
        {
            M[i] = 0;
            for (int j = 0; j < N; j++)
                M[i] += A[i][j] * K[j];
        }

        s = 0;
        sAbs = 0;

        for (int i = 0; i < N; i++)
        {
            s += K[i] * L[i];
            sAbs += M[i] * K[i];
        }
        if (s == sAbs)
            s = 1;
        else
            s = s / sAbs;

        for (int i = 0; i < N; i++)
            x[i] = xrez[i] - s*L[i];

        Del = abs(x[0] - xrez[0]);

        for (int i = 0; i < N; i++)
        {
            if (abs(x[i] - xrez[i])>Del)
                Del = abs(x[i] - xrez[i]);
            xrez[i] = x[i];
        }
    }
    while (eps < Del);

}


void ashryatovarr::lab8()
{
	double eps = 1e-20;
	double** B = new double* [N];
	for (int i = 0; i < N; i++) {
		B[i] = new double[N];
	}

	while (true) {
		double norm = 0;
		int imax = 0;
		int jmax = 1;
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				if (abs(A[i][j]) > abs(A[imax][jmax])) {
					imax = i;
					jmax = j;
				}
				norm += A[i][j] * A[i][j];
			}
		}

		if (sqrt(norm) < eps) {
			break;
		}

		double fi = 0.5 * atan(2 * A[imax][jmax] / (A[imax][imax] - A[jmax][jmax]));

		for (int i = 0; i < N; i++) {
			B[i][imax] = A[i][imax] * cos(fi) + A[i][jmax] * sin(fi);
			B[i][jmax] = A[i][jmax] * cos(fi) - A[i][imax] * sin(fi);
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (j != imax && j != jmax) {
					B[i][j] = A[i][j];
				}
			}
		}

		for (int j = 0; j < N; j++) {
			A[imax][j] = B[imax][j] * cos(fi) + B[jmax][j] * sin(fi);
			A[jmax][j] = B[jmax][j] * cos(fi) - B[imax][j] * sin(fi);
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (i != imax && i != jmax) {
					A[i][j] = B[i][j];
				}
			}
		}
	}

	for (int i = 0; i < N; i++) {
		x[i] = A[i][i];
	}
}


void ashryatovarr::lab9()
{
	
  int n = N;
  double *y = new double[n];
  double *y_next = new double[n];
  double eps =1e-2;
  double lyambda = 1;
  double lyambda_next = 0;

  for(int i = 0; i<n; i++)
    y[i] = b[i];

  do{
  for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
        y_next[i] += A[i][j]*y[j];
      }
  }
  lyambda = lyambda_next;

  for(int i=0; i<n; i++){
    if(y[i]!= 0 && y_next[i] != 0){
      lyambda_next = y_next[i]/y[i];
      break;
    }
  }
  for (int i=0; i<n; i++)
   	y[i]=y_next[i];
}while(fabs(lyambda_next - lyambda)>eps);

cout<<"Result: "<<lyambda_next << endl;

}


std::string ashryatovarr::get_name()
{
    return "Ashryatova R.R.";
}
