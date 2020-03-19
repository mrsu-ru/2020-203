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

	alpha = new double[N];  //массивы для прогоночных коэффициентов
	beta = new double[N];
	
	alpha[0] = A[0][1]/A[0][0];		//начальные коэффициенты
	beta[0] = -b[0]/A[0][0];
	
	//прямой ход - нахождение коэффициентов
	for (int i=1;i<N;i++)
	{
		alpha[i] = A[i][i+1]/(A[i][i]-A[i][i-1]*alpha[i-1]);
		beta[i] = (A[i][i-1]*beta[i-1]-b[i])/(A[i][i]-A[i][i-1]*alpha[i-1])
	
	}
	
	
	//обратный ход - нахождение решений 
	x[N-1] = b[N-1];
	//x[N-2] = alpha[N-2]*x[N-1]+beta[N-2];
	
	for (int j= N-2;j>=0;j--)
	{
	x[j] = alpha[j]*x[j+1]+beta[j];
	}
 	
	
}



/**
 * Метод простых итераций
 */
void kazakovais::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void kazakovais::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void kazakovais::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void kazakovais::lab7()
{

}


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
