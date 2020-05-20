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
	int max, i, j, k;
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
			A[i] = A[max];
			A[max] = vspom;

			vspom_1 = b[i];
			b[i] = b[max];
			b[max] = vspom_1;
		}

		for (j = N - 1; j > i; j--)
		{
			A[i][j] /= A[i][i];
		}
		b[i] /= A[i][i];
		A[i][i] = 1;

		for (j = i + 1; j < N; j++)
		{
			for (k = N - 1; k > i; k--)
			{
				A[j][k] = A[j][k] - A[j][i] * A[i][k];
			}
			b[j] = b[j] - A[j][i] * b[i];
			A[j][i] = 0;
		}
		delete[]vspom;
	}

	//обратный ход
	x[N - 1] = b[N - 1];
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

	delete[]alpha;
	delete[]beta;
}



/**
 * Метод Холецкого (метод квадратного корня)
 */
void kazakovais::lab4()
{
	//данный метод применяется для симметричных матриц
	double *d;
	double *sum; //массив вспомогательных сумм для вычисления d[i] и s[i][i]
	double vspom = 0;
	double *y; //вспомогательный массив для вычисления решений СЛАУ

	int i, j, k;

	d = new double[N];	//массив для коэффициентов диагональной матрицы
	sum = new double[N];
	y = new double[N];
	double **s = new double*[N];

	for (i = 0; i < N; i++)
	{
		s[i] = new double[N];
	}

	//первый этап - поиск LU-разложения
	//начальная инициализация массива вспомогательных сумм
	for (i = 0; i < N; i++)
	{
		sum[i] = A[i][i];
	}

	if (sum[0] > 0)
	{
		d[0] = 1;
	}
	else
	{
		d[0] = -1;
	}
	s[0][0] = sqrt(fabs(sum[0]));
	for (j = 1; j < N; j++)
	{
		s[0][j] = A[0][j] / (d[0] * s[0][0]);
	}

	for (i = 1; i < N; i++)
	{
		for (k = 0; k < i; k++)
		{
			sum[i] -= d[k] * pow(s[k][i], 2);
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
		for (j = i + 1; j < N; j++)
		{
			for (k = 0; k < i; k++)
			{
				vspom += d[k] * s[k][i] * s[k][j];
			}

			s[i][j] = (A[i][j] - vspom) / (d[i] * s[i][i]);
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

	delete[]d;
	delete[]sum;
	delete[]y;
	delete[]s;
}



/**
* Метод Зейделя
*/
void kazakovais::lab5()
{
	//Первоначальный вид СЛАУ A*x=b
	//Приводим её к виду: x=alpha*x+beta
	double eps = 1.0e-20;  //задаваемая точность
	double **alpha = new double*[N];
	double *beta = new double[N];
	int i, j, k;

	for (i = 0; i < N; i++)
	{
		alpha[i] = new double[N];
	}

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			if (i == j)
			{
				alpha[i][j] = 0;
			}
			else
			{
				alpha[i][j] = -A[i][j] / A[i][i];
			}
		}
		beta[i] = b[i] / A[i][i];
	}

	for (i = 0; i < N; i++)
	{
		x[i] = beta[i];
	}

	double *xx = new double[N]; //массив для итерационно вычисленных значений

	for (i = 0; i < N; i++)
	{
		xx[i] = beta[i];
	}

	for (i = 0; i < N; i++)
	{

		for (k = i; k < N; k++)
		{
			xx[i] += alpha[i][k] * x[k];
		}

		for (j = 0; j < i; j++)
		{
			xx[i] += alpha[i][j] * xx[j];
		}

	}

	double Norm = abs(xx[0] - x[0]);  //норма (для определения точки выхода из цикла)

	for (i = 1; i < N; i++)
	{
		if ((abs(xx[i] - x[i])) > Norm)
		{
			Norm = abs(xx[i] - x[i]);
		}
	}

	while (Norm > eps)
	{
		for (i = 0; i < N; i++)
		{
			x[i] = xx[i];
			xx[i] = 0;
		}

		for (i = 0; i < N; i++)
		{

			for (k = i; k < N; k++)
			{
				xx[i] += alpha[i][k] * x[k];
			}

			for (j = 0; j < i; j++)
			{
				xx[i] += alpha[i][j] * xx[j];
			}


			xx[i] += beta[i];
		}

		Norm = abs(xx[0] - x[0]);

		for (i = 1; i < N; i++)
		{
			if ((abs(xx[i] - x[i])) > Norm)
			{
				Norm = abs(xx[i] - x[i]);
			}
		}
	}

	for (i = 0; i < N; i++)
	{
		x[i] = xx[i];
	}

	delete[]alpha;
	delete[]beta;
	delete[]xx;
}



/**
 * Метод минимальных невязок
 */
void kazakovais::lab6()
{
	double eps = 1e-20;
	double *x1 = new double[N];	//x1 - вспомогательный вектор для решений
	double *r = new double[N];	//r - вектор невязки
	double *Ar = new double[N];	//Ar - произведение матрицы A на вектор невязок
	double Ar2 = 0;
	double Arr = 0;
	double t = 0;	//приближение
	double norma = 0;	//норма
	int i, j;

	for (i = 0; i < N; i++)
	{
		x[i] = b[i];
		x1[i] = 0;
		Ar[i] = 0;
	}

	for (i = 0; i < N; i++)
	{
		r[i] = b[i];
		for (j = 0; j < N; j++)
		{
			r[i] -= A[i][j] * x[j];
		}
	}

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			Ar[i] += A[i][j] * r[j];
		}
	}

	for (i = 0; i < N; i++)
	{
		Ar2 += Ar[i] * Ar[i];
		Arr += Ar[i] * r[i];
	}

	t = -Arr / Ar2;

	for (i = 0; i < N; i++)
	{
		x1[i] = x[i] - t * r[i];
	}

	norma = abs(x1[0] - x[0]);
	for (i = 1; i < N; i++)
	{
		if ((abs(x1[i] - x[i])) > norma)
		{
			norma = abs(x1[i] - x[i]);
		}
	}

	while (norma > eps)
	{
		for (i = 0; i < N; i++)
		{
			x[i] = x1[i];
			x1[i] = 0;
			Ar[i] = 0;
			r[i] = 0;
		}
		Ar2 = 0;
		Arr = 0;
		t = 0;
		norma = 0;

		for (i = 0; i < N; i++)
		{
			r[i] = b[i];
			for (j = 0; j < N; j++)
			{
				r[i] -= A[i][j] * x[j];
			}
		}

		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				Ar[i] += A[i][j] * r[j];
			}
		}

		for (i = 0; i < N; i++)
		{
			Ar2 += Ar[i] * Ar[i];
			Arr += Ar[i] * r[i];
		}

		t = -Arr / Ar2;

		for (i = 0; i < N; i++)
		{
			x1[i] = x[i] - t * r[i];
		}

		norma = abs(x1[0] - x[0]);
		for (i = 1; i < N; i++)
		{
			if ((abs(x1[i] - x[i])) > norma)
			{
				norma = abs(x1[i] - x[i]);
			}
		}
	}

	for (i = 0; i < N; i++)
	{
		x[i] = x1[i];
	}

	delete[]x1;
	delete[]r;
	delete[]Ar;
}



/**
 * Метод сопряженных градиентов
 */
void kazakovais::lab7()
{
	//данный метод работает только для симметричных матриц коэффициентов
	double eps = 1e-24;
	double *r = new double[N];
	double *r1 = new double[N];
	double *z = new double[N];
	double alpha = 0;
	double beta = 0;
	double r0scal = 0;
	double r1scal = 0;
	double *Az = new double[N];
	double Azzscal = 0;
	double Norma = 0;
	double norma1 = 0;
	double norma2 = 0;
	int i, j;

	for (i = 0; i < N; i++)
	{
		x[i] = b[i];	//задаём начальное приближение
	}

	for (i = 0; i < N; i++)
	{
		r[i] = b[i];
		for (j = 0; j < N; j++)
		{
			r[i] -= A[i][j] * x[j];		//рассчитываем вектор невязки
		}
	}

	for (i = 0; i < N; i++)
	{
		z[i] = r[i];
	}

	for (i = 0; i < N; i++)
	{
		r0scal += r[i] * r[i];
		for (j = 0; j < N; j++)
		{
			Az[i] += A[i][j] * z[j];
		}
		Azzscal += Az[i] * z[i];
	}
	alpha = r0scal / Azzscal;	//коэффициент для нахождения x[i] на следующей итерации

	for (i = 0; i < N; i++)
	{
		x[i] = x[i] + alpha * z[i];
		r1[i] = r[i] - alpha * Az[i];
	}

	for (i = 0; i < N; i++)
	{
		r1scal += r1[i] * r1[i];
	}
	beta = r1scal / r0scal;

	for (i = 0; i < N; i++)
	{
		z[i] = r1[i] + beta * z[i];
	}

	norma1 = sqrt(r0scal);
	for (i = 0; i < N; i++)
	{
		norma2 += b[i] * b[i];
	}
	norma2 = sqrt(norma2);
	Norma = norma1 / norma2;

	while (Norma >= eps)
	{
		alpha = 0;
		beta = 0;
		r0scal = r1scal;
		r1scal = 0;
		for (i = 0; i < N; i++)
		{
			Az[i] = 0;
			r[i] = r1[i];
			r1[i] = 0;
		}
		Azzscal = 0;
		norma1 = 0;
		Norma = 0;

		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				Az[i] += A[i][j] * z[j];
			}
			Azzscal += Az[i] * z[i];
		}
		alpha = r0scal / Azzscal;

		for (i = 0; i < N; i++)
		{
			x[i] = x[i] + alpha * z[i];
			r1[i] = r[i] - alpha * Az[i];
		}

		for (i = 0; i < N; i++)
		{
			r1scal += r1[i] * r1[i];
		}
		beta = r1scal / r0scal;

		for (i = 0; i < N; i++)
		{
			z[i] = r1[i] + beta * z[i];
		}

		norma1 = sqrt(r0scal);
		Norma = norma1 / norma2;
	}
}


/**
*	Метод вращения для нахождения собственных значений матрицы
*/
void kazakovais::lab8()
{

}


/**
*	Нахождение наибольшего по модулю собственного значения матрицы
*/
void kazakovais::lab9()
{
	//Для данной матрицы наибольшее по модулю собственное значение вещественное и простое
	double eps = 1e-3;
	double *y0 = new double[N];	//произвольный начальный вектор
	double *y1 = new double[N]; //вектор для итерационно вычисленных значений
	double y00 = 0;
	double y10 = 0;
	double lambda_max0 = 0;
	double lambda_max1 = 0;
	double delta_lambda = 0;
	int i, j;

	for (i = 0; i < N; i++)
	{
		y0[i] = 1;
		y1[i] = 0;
	}
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			y1[i] += A[i][j] * y0[j];
		}
	}
	y00 = y0[0];

	for (i = 0; i < N; i++)
	{
		if (y1[i] != 0)
		{
			y10 = y1[i];
			break;
		}
	}

	lambda_max0 = y10 / y00;
	delta_lambda = lambda_max0;

	while (delta_lambda > eps)
	{
		for (i = 0; i < N; i++)
		{
			y0[i] = y1[i];
			y1[i] = 0;
		}
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				y1[i] += A[i][j] * y0[j];
			}
		}
		for (i = 0; i < N; i++)
		{
			if ((y1[i] != 0) && (y0[i] != 0))
			{
				y00 = y0[i];
				y10 = y1[i];
				break;
			}
		}
		lambda_max1 = y10 / y00;
		delta_lambda = fabs(lambda_max1 - lambda_max0);
		lambda_max0 = lambda_max1;
	}
	cout << "Наибольшее по модулю собственное значение матрицы A: " << lambda_max1<<endl;
}


std::string kazakovais::get_name()
{
	return "Kazakova I.S.";
}
