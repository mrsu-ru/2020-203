#include "zhalninrv.h"
#include "parshinad.h"

/**
 * Введение в дисциплину
 */
void parshinad::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void parshinad::lab2()
{
	double eps = 1e-12;
	for (int i = 0; i < N; i++) {
		if (fabs(A[i][i]) < eps) {
			for (int l = i + 1; l < N; l++) {
				if (fabs(A[l][i]) > eps) {
					double* tmpS = A[i];
					A[i] = A[l];
					A[l] = tmpS;
					swap(b[i], b[l]);
					break;
				}
			}
		}


		if (fabs(A[i][i] - 1.0) > eps) {
			double tmp = A[i][i];
			for (int k = i + 1; k < N; k++) {
				A[i][k] /= tmp;
			}
			A[i][i] = 1;
			b[i] /= tmp;
		}

		for (int j = i + 1; j < N; j++) {
			double tmp2 = A[j][i];
			for (int k = i; k < N + 1; k++) {
				A[j][k] -= tmp2 * A[i][k];
			}
			b[j] -= tmp2 * b[i];
		}
	}

	for (int i = 1; i < N; i++) {
		for (int j = N - i - 1; j >= 0; j--) {
			double tmp = A[j][N - i];
			b[j] -= tmp * b[N - i];
		}
	}

	for (int i = 0; i < N; i++) {
		x[i] = b[i];
	}
}



/**
 * Метод прогонки
 */
void parshinad::lab3()
{
	double *alpha = new double[N];
    double *beta = new double[N];
    
    alpha[0] = -A[0][1]/A[0][0];
    beta[0] = b[0]/A[0][0];
    
    for (int i = 1; i <= N - 1; i++){
        alpha[i] = -A[i][i+1]/(A[i][i] + A[i][i-1]*alpha[i-1]);
        beta[i] = (-A[i][i-1]*beta[i-1] + b[i])/(A[i][i]+A[i][i-1]*alpha[i-1]);
    }
    
    x[N-1] = beta[N-1];
    
    for (int i = N - 2; i >= 0; i--){
        x[i] = alpha[i]*x[i+1] + beta[i];
    }
}



/**
 * Метод Холецкого
 */
void parshinad::lab4()
{
	double** S = new double* [N];
    double* D = new double[N];

    for (int i = 0; i < N; i++)
    {
        S[i] = new double[N];
        memset(S[i], 0, sizeof(double) * N);
    }



    for (int i = 0; i < N; i++)
    {
        double sum1 = A[i][i];

        for (int k = 0; k <= i - 1; k++) {
            sum1 -= D[k] * S[k][i] * S[k][i];
        }

        D[i]=sum1/fabs(sum1);

        S[i][i] = sqrt(sum1);

        for (int j = i + 1; j < N; j++)
        {
            double sum2 = 0;
            for (int k = 0; k <= j - 1; k++) {
                sum2 += S[k][i] * S[k][j]* D[k];
            }

            S[i][j] = (A[i][j] - sum2) / (D[i] * S[i][i]);
        }
    }

    for (int i = 0; i < N; i++)
    {
        b[i] /= S[i][i];
        x[i] = b[i];

        for (int j = i + 1; j < N; j++) {
            b[j] -= b[i] * S[i][j];
        }
    }



    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            S[i][j] *= D[i];
        }
    }


    for (int i = N - 1; i >= 0; i--)
    {
        x[i] /= S[i][i];

        for (int j = i - 1; j >= 0; j--) {
            x[j] -= x[i] * S[j][i];
        }
    }
}



/**
 * Метод Якоби или Зейделя
 */
void parshinad::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void parshinad::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void parshinad::lab7()
{

}


void parshinad::lab8()
{

}


void parshinad::lab9()
{

}


std::string parshinad::get_name()
{
  return "Parshin A.D.";
}
