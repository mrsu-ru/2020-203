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
    const double eps = 1e-18;
	double* xPrev = new double[N];
	for (int i = 0; i < N; i++)
		xPrev[i] = 0;
	double* r = new double[N];
	double* rPrev = new double[N];
	double* A_xPrev = new double[N];
	double scal_r, rPrev_rPrev;
	double Ax_x;
	double alpha, beta;
	do{
		for (int i = 0; i < N; i++)
			rPrev[i] = -b[i];

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				rPrev[i] += A[i][j] * x[j];

		for (int i = 0; i < N; i++)
			xPrev[i] = rPrev[i];

		for (int i = 0; i < N; i++)
			A_xPrev[i] = 0;

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				A_xPrev[i] += A[i][j] * xPrev[j];

		rPrev_rPrev = 0;
		for (int i = 0; i < N; i++)
			rPrev_rPrev += rPrev[i] * rPrev[i];

		Ax_x = 0;
		for (int i = 0; i < N; i++)
			Ax_x += A_xPrev[i] * xPrev[i];

		alpha = rPrev_rPrev / Ax_x;

		for (int i = 0; i < N; i++)
			x[i] = x[i] + alpha * xPrev[i];

		for (int i = 0; i < N; i++)
			r[i] = rPrev[i] - alpha * A_xPrev[i];

		scal_r = 0;
		for (int i = 0; i < N; i++)
			scal_r += r[i] * r[i];

		beta = scal_r / rPrev_rPrev;

		for (int i = 0; i < N; i++)
			xPrev[i] = r[i] + beta * xPrev[i];

		for (int i = 0; i < N; i++)
			rPrev[i] = r[i];

	} while (scal_r > eps);

	delete[]xPrev;
	delete[]rPrev;
	delete[]r;
	delete[]A_xPrev;
}


void kirdyushkindv::lab8()
{

}


void kirdyushkindv::lab9()
{

}


std::string kirdyushkindv::get_name()
{
  return "Kirdyushkin D.V.";
}
