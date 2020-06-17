#include "manindi.h"

/**
 * Введение в дисциплину
 */
void manindi::lab1()
{
  
cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void manindi::lab2()
{
	int *sw1 = new int[N];
	for (int i = 0; i < N; i++)
	{
		sw1[i] = i;
	}
	for (int i = 0; i < N; i++)
	{
		int Maxi = i;
		for (int j = i + 1; j < N; j++)
		{
			if (abs(A[j][i]) > abs(A[Maxi][i]))
			{
				Maxi = j;
			}
		}
		if (Maxi != i)
		{
			swap(A[Maxi], A[i]);
			swap(b[Maxi], b[i]);
			swap(sw1[Maxi], sw1[i]);
		}
		b[i] /= A[i][i];
		for (int j = N - 1; j >= i; j--)
		{
			A[i][j] /= A[i][i];
		}

		for (int j = 0; j < N; j++)
		{
			if (j != i)
			{
				b[j] -= A[j][i] * b[i];
				for (int k = N - 1; k >= i; k--)
				{
					A[j][k] -= A[j][i] * A[i][k];
				}

			}

		}
	}
	for (int i = 0; i < N; i++) {
		x[i] = b[i];
	}
	for (int i = 0; i < N; i++) {
		if (sw1[i] != i) {
			swap(A[sw1[i]], A[i]);
			swap(b[sw1[i]], b[i]);
		}
	}
}



/**
 * Метод прогонки
 */
void manindi::lab3()
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
        
}



/**
 * Метод Холецкого
 */
void manindi::lab4()
{
  double *d;
  double *sum;
  double vs = 0;
  double *m;
  int i,j,k;
  d = new double[N];
  sum = new double[N];
  m = new double[N];
  double **s = new double*[N];

    for (i=0; i<N;i++){
       s[i] = new double[N];
    }

    for (i = 0; i < N; i++){
       sum[i] = A[i][i];
    }
      if (sum[0]>0){
        d[0]=1;
    }else {
        d[0]=-1;
    }
  s[0][0]=sqrt(fabs(sum[0]));

    for (j=1;j<N;j++){
       s[0][j]=A[0][j]/(d[0]*s[0][0]);
    }

    for (i = 1; i < N; i++){

    for (k = 0; k < i; k++){
       sum[i] -= d[k] * pow(s[k][i],2);
    }
       if (sum[i] > 0){
         d[i] = 1;
     }else{
         d[i] = -1;
    }
  s[i][i] = sqrt(fabs(sum[i]));

    for (j = i+1; j < N; j++){
    
    for (k = 0; k < i; k++){
        vs += d[k] * s[k][i] * s[k][j];
    }

  s[i][j] = (A[i][j]-vs)/(d[i]*s[i][i]);

    vs = 0;
    }}
  m[0] = b[0] / s[0][0];
    
     for (i = 1; i < N; i++){
     
      for (k = 0; k < i; k++){
      vs += s[k][i] * m[k];
   }
    m[i] = (b[i] - vs) / s[i][i];

     vs = 0;
   }
  x[N - 1] = m[N - 1] / (d[N - 1] * s[N - 1][N - 1]);
    for (i = N - 2; i >= 0; i--){
    for (k = i + 1; k < N; k++){
      vs += s[i][k] * x[k];
   }
   x[i] = (m[i] - d[i] * vs) / (d[i] * s[i][i]);
   vs = 0;
    }

               
}



/**
 * Метод Якоби или Зейделя
 */
void manindi::lab5()
{
   double epsilon = 1e-20;

  for (int i = 0; i < N; i++) {
     x[i] = 0;
   }
   double *prev = new double[N];
   double nor = 0;
      do {
      for (int i = 0; i < N; i++) {
      prev[i] = x[i];
   }
    for (int i = 0; i < N; i++) {
     double result = b[i];

      for (int j = 0; j < N; j++) {
      if (i != j) {
     result -= (A[i][j] * prev[j]);
          } 
      }
     x[i] = result / A[i][i];
   }
     nor = 0;
     for (int i = 0; i < N; i++) {
      if (abs(prev[i] - x[i]) > nor) {
      nor = abs(prev[i] - x[i]);
       }
     }
    } while (nor > epsilon);
      delete[] prev;
}



/**
 * Метод минимальных невязок
 */
void manindi::lab6()
{
   double eps = 1e-20;
   double *px = new double[N];
      
     for (int i = 0; i < N; i++) {
           px[i] = 0;
	}

   double *r = new double[N];
           while (true) {

      for (int i = 0; i < N; i++) {
        r[i] = b[i];
      for (int j = 0; j < N; j++) {
        r[i] -= (A[i][j] * px[j]);
         }
        }
        double tau = 0;
        double tmp = 0;

      for (int i = 0; i < N; i++) {
        double Ar = 0;
      for (int j = 0; j < N; j++) {
         Ar += (A[i][j] * r[j]);
         }
         tau += (Ar * r[i]);
         tmp += (Ar * Ar);
         }
         tau /= tmp;

       for (int i = 0; i < N; i++) {
        x[i] = px[i] + tau * r[i];

	}
       double error = 0;

        for (int i = 0; i < N; i++) {
        if (abs(x[i] - px[i]) > error) {
            error = abs(x[i] - px[i]);
           }
          }
        if (error < eps) {

         break;
       }

       for (int i = 0; i < N; i++) {
       px[i] = x[i];
        }
      }
}



/**
 * Метод сопряженных градиентов
 */
void manindi::lab7()
{
const double eps = 1e-18;
double *xPrev = new double[N];

	memset(xPrev, 0, sizeof(double) * N);

double *r = new double[N];

	for (int i = 0; i < N; i++)
		r[i] = -b[i];
double *Ar = new double[N];
	
	memset(Ar, 0, sizeof(double) * N);

	for (int i = 0; i < N; i++) 

		for (int j = 0; j < N; j++)
			Ar[i] += A[i][j] * r[j];

double rrScalar = 0;

	for (int i = 0; i < N; i++)
		rrScalar += r[i] * r[i];

double rArScalar = 0;

	for (int i = 0; i < N; i++)
		rArScalar += Ar[i] * r[i];

double alpha = 1;
double tau = rrScalar / rArScalar;

	for (int i = 0; i < N; i++)
		x[i] = tau * b[i];

	for (int i = 0; i < N; i++) 

		for (int j = 0; j < N; j++)
			r[i] += A[i][j] * x[j];

rrScalar = 0;

	for (int i = 0; i < N; i++)
		rrScalar += r[i] * r[i];

	while (rrScalar > eps * eps)
	{
double rArScalarPrev = rArScalar;

memset(Ar, 0, sizeof(double) * N);

	for (int i = 0; i < N; i++) 

		for (int j = 0; j < N; j++)
			Ar[i] += A[i][j] * r[j];
rArScalar = 0;

	for (int i = 0; i < N; i++)
		rArScalar += Ar[i] * r[i];

double tauPrev = tau;
tau = rrScalar / rArScalar;
alpha = 1.0 / (1 - tau / tauPrev / alpha * rrScalar / rArScalarPrev);

	for (int i = 0; i < N; i++)
	{
double tmp = x[i];
x[i] = alpha * x[i] + (1 - alpha) * xPrev[i] - tau * alpha * r[i];
xPrev[i] = tmp;

	}

	for (int i = 0; i < N; i++)
		r[i] = -b[i];

	for (int i = 0; i < N; i++) 
	
		for (int j = 0; j < N; j++)
			r[i] += A[i][j] * x[j];

		rrScalar = 0;
	
	for (int i = 0; i < N; i++)
		rrScalar += r[i] * r[i];

	}

delete []xPrev;

delete []r;

delete []Ar;
}


void manindi::lab8()
{
double eps = 1e-20;
double** H = new double* [N];

	for (int i = 0; i < N; i++) {
	H[i] = new double[N];
	}

	while (true) {
	double n = 0;
	int i_max = 0;
	int j_max = 1;
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
			if (abs(A[i][j]) > abs(A[i_max][j_max])) {
			i_max = i;
			j_max = j;
			}
		n += A[i][j] * A[i][j];
		}
	}

	if (sqrt(n) < eps) {

	break;

	}

double fi = 0.5 * atan(2 * A[i_max][j_max] / (A[i_max][i_max] - A[j_max][j_max]));

	for (int i = 0; i < N; i++) {
		H[i][i_max] = A[i][i_max] * cos(fi) + A[i][j_max] * sin(fi);
		H[i][j_max] = A[i][j_max] * cos(fi) - A[i][i_max] * sin(fi);
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
		if (j != i_max && j != j_max) {
		H[i][j] = A[i][j];
				}
		}
	}

	for (int j = 0; j < N; j++) {
		A[i_max][j] = H[i_max][j] * cos(fi) + H[j_max][j] * sin(fi);
		A[j_max][j] = H[j_max][j] * cos(fi) - H[i_max][j] * sin(fi);
		}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
		if (i != i_max && i != j_max) {
		A[i][j] = H[i][j];
			}
		}
		}
	}

	for (int i = 0; i < N; i++) {
		x[i] = A[i][i];
	}
}


void manindi::lab9()
{
double eps = 1e-3; 
double *y = new double[N];
double *yNext = new double[N]; 

	for (int i = 0; i < N; i++){
	y[i] = 1; 
	}
double delta = 0; 
double maxLambda = 0; 

	do{
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
			yNext[i] += A[i][j] * y[j]; 
			}
		}
double lambda;
	for (int i = 0; i < N; i++){
	if(fabs(y[i]) > eps && fabs(yNext[i]) > eps){
	lambda = yNext[i]/y[i];

	break;
		}
	} 
delta = fabs(lambda - maxLambda); 
maxLambda = lambda; 
		
	for (int i = 0; i < N; i++) {
	y[i] = yNext[i]; 
	}
	
	memset(yNext, 0, sizeof(double) * N);
	
	}while(delta > eps);

	cout << "maxLambda = " << maxLambda << endl;
}


std::string manindi::get_name()
{
  return "Manin D.I.";
}
