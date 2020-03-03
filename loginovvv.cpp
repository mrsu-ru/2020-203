#include "loginovvv.h"

/**
 * Введение в дисциплину
 */
void loginovvv::lab1()
{
	cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void loginovvv::lab2()
{
	//cout << "hello world!" << endl;

	/////приводим к треугольнику
	int z;
	for (int i = 0; i < N; i++)
	{
		z = i;
		//ишим наибольший по модулю 
		for (int m = i+1; m < N; m++)
			if (abs(A[m][i]) > abs(A[z][i]))
				z = m;
		if (A[z][i] == 0)
			return;
		//меняем строки 
		double* vrem;
    	double v;
		vrem = A[i];
		A[i] = A[z];
		A[z] = vrem;
    	v = b[i];
    	b[i] = b[z];
    	b[z] = v;
		//делим строку , чтоб получить единичку
    	b[i]/=A[i][i];
		for (int j = N-1; j > i; A[i][j--] /= A[i][i]);
		A[i][i] = 1;
		//вычитаем строку , чтоб получить нолики
		for (int j = i + 1; j < N; j++)
		{
            b[j]-=b[i]*A[j][i];
			for (int k = N-1; k > i; k--)
				A[j][k] -= A[i][k] * A[j][i];
			A[j][i] = 0;
		}
		
	}
	//находим решение 
	x[N - 1] = b[N-1];
	for (int i = N - 2; i >= 0; i--)
	{
		x[i] = b[i];
		for (int j = i + 1; j < N; j++)
			x[i] -= A[i][j] * x[j];
	}

}



/**
 * Метод прогонки
 */
void loginovvv::lab3()
{
	double *a = new double[N];
	double *v = new double[N];
	a[1]=-A[0][1]/A[0][0];
	v[1]=b[0]/A[0][0];
	for(int i=2; i<N; i++)
	{
		a[i]=-A[i-1][i]/(A[i-1][i-2]*a[i-1]+A[i-1][i-1]);
		v[i]=(b[i-1]-A[i-1][i-2]*v[i-1])/(A[i-1][i-2]*a[i-1]+A[i-1][i-1]);
	}
	x[N-1]=(b[N-1]-A[N-1][N-2]*v[N-1])/(A[N-1][N-1]+A[N-1][N-2]*a[N-1]);
	for(int i=N-2; i>-1; i--)
	{
		x[i]=a[i+1]*x[i+1]+v[i+1];

	}

}



/**
 * Метод простых итераций
 */
void loginovvv::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void loginovvv::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void loginovvv::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void loginovvv::lab7()
{

}


void loginovvv::lab8()
{

}


void loginovvv::lab9()
{

}


std::string loginovvv::get_name()
{
  return "Loginov V. V.";
}
