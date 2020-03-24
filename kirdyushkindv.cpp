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

}



/**
 * Метод Якоби или Зейделя
 */
void kirdyushkindv::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void kirdyushkindv::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void kirdyushkindv::lab7()
{

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
