#include "kirdyushkindv.h"

/**
 * Введение в дисциплину
 */
void kirdyushkindv::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kirdyushkindv::lab2()
{
int max_row,i,j,k;
double sum;
double* tmp1; double tmp2;
	sum = 0;
	for (i = 0; i < N; i++)              //прямой ход
	{
		max_row = i;
		for (j = i + 1; j < N; j++)
		{
			if (fabs(A[j][i]) > fabs(A[max_row][i]))
				max_row = j;

		}
			if (max_row != i)
			{
	               tmp1 = A[i];
	               A[i]=A[max_row];
	               A[max_row] = tmp1;

	               tmp2 = b[i];
	               b[i] = b[max_row];
	               b[max_row]= tmp2;
			}

		for (j = N-1; j > i; j--)
			A[i][j] /= A[i][i];

		b[i]/=A[i][i];
		A[i][i] = 1;

		for (j = i + 1; j < N; j++)
		{
			for (k = N-1; k > i; k--)
				A[j][k] = A[j][k] - A[j][i] * A[i][k];

			b[j]=b[j]-A[j][i]*b[i];
			A[j][i] = 0;
		}
	}


	x[N-1]=b[N-1];                      //обратный ход
	for (i = N - 2; i > -1; i--)
	{
		for (j = i + 1; j < N; j++)
			sum += x[j] * A[i][j];

		x[i] = b[i] - sum;
		sum = 0;
	}

}



/**
 * Метод прогонки
 */
void kirdyushkindv::lab3()
{
double alpha[N-1], beta[N];
    alpha[0]=-A[0][1]/A[0][0];
    beta[0]=b[0]/A[0][0];
    for (int i=1; i<N-1; i++)
    {
        alpha[i]=-A[i][i+1]/(A[i][i]+A[i][i-1]*alpha[i-1]);
        beta[i]=(b[i]-A[i][i-1]*beta[i-1])/(A[i][i]+A[i][i-1]*alpha[i-1]);
    }
    beta[N-1]=(b[N-1]-A[N-1][N-2]*beta[N-2])/(A[N-1][N-1]+A[N-1][N-2]*alpha[N-2]);

    x[N-1]=beta[N-1];
    for (int i=N-2; i>=0; i--)
    {
        x[i]=alpha[i]*x[i+1]+beta[i];
    }
}



/**
 * Метод Холецкого
 */
void kirdyushkindv::lab4()
{
    double **S = new double*[N];
        for (int i=0; i<N; i++)
        {
            S[i]=new double[N];
            for(int j=0; j<N; j++)
                S[i][j]=0;
        }
    double *D = new double[N];
        if (A[0][0]>0)    D[0]=1;
        else              D[0]=-1;

        S[0][0]=sqrt(fabs(A[0][0]));

        for (int i=1; i<N; i++)
        {
            S[0][i]=A[0][i]/(S[0][0]*D[0]);
        }

        for (int i=1; i<N; i++)
        {
            double sum =0;
            for (int j=0; j<i; j++)
                sum+=S[j][i]*S[j][i]*D[j];
            if (A[i][i]-sum>=0)    D[i]=1;
            else                    D[i]=-1;

            S[i][i]=sqrt((A[i][i]-sum)*D[i]);

            for (int j=i+1; j<N; j++)
            {
                double l = 0;
                for (int k=0; k<j; k++)
                    l+=S[k][i]*S[k][j]*D[k];

                S[i][j]=(A[i][j]-l)/(S[i][i]*D[i]);
            }
        }

    double *y = new double[N];
        y[0]=b[0]/S[0][0];
        for (int i=1; i<N; i++)
        {
            double sum = 0;
            for (int j=0; j<i; j++)
                sum+=y[j]*S[j][i];
            y[i]=(b[i]-sum)/S[i][i];
        }

        x[N-1]=y[N-1]/(S[N-1][N-1]*D[N-1]);

        for (int i=N-2; i>=0; i--)
        {
            double sum =0;
            for (int j=i+1; j<N; j++)
                sum+=x[j]*S[i][j]*D[j];
            x[i]=(y[i]-sum)/(S[i][i]*D[i]);
        }
}



/**
 * Метод Якоби или Зейделя
 */
void kirdyushkindv::lab5()
{
    double eps = 1e-20;
	double *prev_x = new double[N];
	double norm = 0;

	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}

	do {
		for (int i = 0; i < N; i++)
			prev_x[i] = x[i];

		for (int i = 0; i < N; i++) {
			double result = b[i];
			for (int j = 0; j < N; j++) {
				if (i != j) {
					result -= (A[i][j] * prev_x[j]);
				}
			}

			x[i] = result / A[i][i];
		}

		norm = 0;
		for (int i = 0; i < N; i++) {
			if (abs(prev_x[i] - x[i]) > norm) {
				norm = abs(prev_x[i] - x[i]);
			}
		}
	} while (norm > eps);

	delete[] prev_x;
}




/**
 * Метод минимальных невязок
 */
void kirdyushkindv::lab6()
{
    double eps = 1e-19;
    double err;
    double Ax, Ay[N];
    double y0[N];
    double t, x_prev;
    double sum1, sum2;
    do{
        for (int i=0; i<N; i++){
            Ax=0;
            for (int j=0; j<N; j++){
                Ax+=A[i][j]*x[j];
            }
            y0[i]=b[i]-Ax;
        }
        for (int i=0; i<N; i++){
            Ay[i]=0;
            for (int j=0; j<N; j++){
                Ay[i]+=A[i][j]*y0[j];
            }
        }
        sum1=0; sum2=0;
        for (int i=0; i<N; i++){
            sum1+=y0[i]*Ay[i];
            sum2+=Ay[i]*Ay[i];
        }
        t=sum1/sum2;
        err=0;
        for(int i=0; i<N; i++){
            x_prev=x[i];
            x[i]+=t*y0[i];
            if (abs(x[i]-x_prev)>err){
                err=abs(x[i]-x_prev);
            }
        }
    }while (err>eps);
}



/**
 * Метод сопряженных градиентов
 */
void kirdyushkindv::lab7()
{
    double const eps = 1e-10;
    double r0[N], z0[N], Az[N];
    double err;

    for (int i = 0; i < N; i++) {
        double Ax = 0;
        for (int j = 0; j < N; j++) {
            Ax += A[i][j] * x[j];
        }

        r0[i] = b[i] - Ax;
        z0[i] = r0[i];
    }

    do {
        for (int i = 0; i < N; i++) {
            Az[i] = 0;
            for (int j = 0; j < N; j++) {
                Az[i] += A[i][j] * z0[j];
            }
        }

        double sum1 = 0;
        double sum2 = 0;

        for (int i = 0; i < N; i++) {
            sum1 += r0[i] * r0[i];
            sum2 += Az[i] * z0[i];
        }

        double alpha = sum1 / sum2;
        err = 0;

        for (int i = 0; i < N; i++) {
            double temp = x[i];
            x[i] = x[i] + alpha * z0[i];
            if (abs(temp - x[i]) > err) {
                err = abs(temp - x[i]);
            }
        }

        sum1 = 0;
        sum2 = 0;

        for (int i = 0; i < N; i++) {
            sum2 += r0[i] * r0[i];
            r0[i] = r0[i] - alpha * Az[i];

            sum1 += r0[i] * r0[i];
        }

        double beta = sum1 / sum2;

        for (int i = 0; i < N; i++) {
            z0[i] = r0[i] + beta * z0[i];
        }
    } while (err > eps);
}


void kirdyushkindv::lab8()
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


void kirdyushkindv::lab9()
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


std::string kirdyushkindv::get_name()
{
  return "Kirdyushkin D.V.";
}
