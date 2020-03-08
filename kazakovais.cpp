#include "kazakovais.h"

/**
 * Введение в дисциплину
 */
void kazakovais::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kazakovais::lab2()
{
 	int max,i,j,k;
	double sum;
	sum = 0;
	
	//прямой ход
	for (i = 0; i < N; i++)
	{
		max = i;
		double* vspom; double vspom_1;

		for (j = i + 1; j < N; j++)
		{
			if (fabs(A[j][i]) > fabs(A[max][i]))
			{
				max = j;
			}
		}
			if (max != i)
			{
	               vspom = A[i];
	               A[i]=A[max];
	               A[max] = vspom;
         
	               vspom_1 = b[i];
	               b[i] = b[max];
	               b[max]= vspom_1;
			}

		for (j = N-1; j > i; j--)
		{
			A[i][j] /= A[i][i];
		}
		b[i]/=A[i][i];
		A[i][i] = 1;
		
		for (j = i + 1; j < N; j++)
		{
			for (k = N-1; k > i; k--)
			{
				A[j][k] = A[j][k] - A[j][i] * A[i][k];
			}
			b[j]=b[j]-A[j][i]*b[i];
			A[j][i] = 0;
		}
	}

	//обратный ход
	x[N-1]=b[N-1];
	for (i = N - 2; i > -1; i--)
	{
		for (j = i + 1; j < N; j++)
		{
			sum += x[j] * A[i][j];
		}
		x[i] = b[i] - sum;
		sum = 0;
	}
}



/**
 * Метод прогонки
 */
void kazakovais::lab3()
{

	//данный метод применяется для трёхдиагональных матриц
	double *alpha;
	double *beta;
	alpha = new double[N];  //массивы для прогоночных коэффициентов
	beta = new double[N];

	//прямой ход
	alpha[0] = -A[0][1] / A[0][0];
	beta[0] = b[0] / A[0][0];

	for (int i = 1; i < N - 1; i++)
	{
		alpha[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * alpha[i - 1]);
		beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / (A[i][i] + A[i][i - 1] * alpha[i - 1]);
		
	}

	beta[N - 1] = (b[N - 1] - A[N - 1][N - 2] * beta[N - 2]) / (A[N - 1][N - 1] + A[N - 1][N - 2] * alpha[N - 2]);

	//обратный ход
	x[N - 1] = beta[N - 1];

	for (int i = N - 2; i >= 0; i--)
	{
		x[i] = alpha[i] * x[i + 1] + beta[i];
	}

	
}



/**
 * Метод Холецкого (метод квадратного корня)
 */
void kazakovais::lab4()
{

	//данный метод применяется для симметричных матриц
	double *d;
	//double **s;
	double *sum; //массив вспомогательных сумм (см.лекции!!!)
	double vspom = 0;
	double *y; //вспомогательный массив для вычисления решений СЛАУ

	int i,j,k;
	
	d = new double[N];	//массив для коэффициентов диагональной матрицы
	sum = new double[N];
	y = new double[N];
	double **s = new double*[N];
	
	for (i=0; i<N;i++)
	{
		s[i] = new double[N];
	}
	
	//первый этап - поиск LU-разложения
	//начальная инициализация массива вспомогательных сумм
	for (i = 0; i < N; i++)
	{
		sum[i] = A[i][i];
	}

	if (sum[0]>0)
	{
		d[0]=1;
	}
	else 
	{
		d[0]=-1;
	}
	s[0][0]=sqrt(fabs(sum[0]));
	for (j=1;j<N;j++)
	{
		s[0][j]=A[0][j]/(d[0]*s[0][0]);
	}
	
	for (i = 1; i < N; i++)
	{
		for (k = 0; k < i; k++)
		{
			sum[i] -= d[k] * pow(s[k][i],2);
		}
		if (sum[i] > 0)
		{
			d[i] = 1;
		}
		else
		{
			d[i] = -1;
		}
		s[i][i] = sqrt(fabs(sum[i]));
		for (j = i+1; j < N; j++)
		{
			for (k = 0; k < i; k++)
			{
				vspom += d[k] * s[k][i] * s[k][j];
			}

			s[i][j] = (A[i][j]-vspom)/(d[i]*s[i][i]);
			vspom = 0;
		}
	}

	//второй этап - поиск корней
	y[0] = b[0] / s[0][0];
	for (i = 1; i < N; i++)
	{
		for (k = 0; k < i; k++)
		{
			vspom += s[k][i] * y[k];
		}
		y[i] = (b[i] - vspom) / s[i][i];
		vspom = 0;
	}
	x[N - 1] = y[N - 1] / (d[N - 1] * s[N - 1][N - 1]);
	for (i = N - 2; i >= 0; i--)
	{
		for (k = i + 1; k < N; k++)
		{
			vspom += s[i][k] * x[k];
		}
		x[i] = (y[i] - d[i] * vspom) / (d[i] * s[i][i]);
		vspom = 0;
	}
	
}



/**
*Метод простых итераций
*/
void kazakovais::lab5()
{

}



/**
 * Метод Якоби или Зейделя
 */
void kazakovais::lab6()
{

}



/**
 * Метод минимальных невязок
 */
void kazakovais::lab7()
{

}

/**
 * Метод сопряженных градиентов
 */
void kazakovais::lab8()
{

}


void kazakovais::lab9()
{

}


std::string kazakovais::get_name()
{
  return "Kazakova I.S.";
}
