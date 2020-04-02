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

}



/**
 * Метод минимальных невязок
 */
void manindi::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void manindi::lab7()
{

}


void manindi::lab8()
{

}


void manindi::lab9()
{

}


std::string manindi::get_name()
{
  return "Manin D.I.";
}
