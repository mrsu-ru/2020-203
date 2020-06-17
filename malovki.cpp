#include "malovki.h"

using namespace std;

/**
 * Введение в дисциплину
 */
void malovki::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void malovki::lab2()
{
	float max = 0;
	int index;
	//прямой ход
	for (int i = 0; i < N; i++)
	{
		index = i;
		max = abs(A[i][i]);
		for (int k = i; k < N; k++)
		{
			if (max < abs(A[k][i]))
			{
				max = abs(A[k][i]);
				index = k;
			}
		}
		if (i != index) {
			double* vrem;
			vrem = A[index];
			A[index] = A[i];
			A[i] = vrem;
			double temp = b[index];
			b[index] = b[i];
			b[i] = temp;
		}
		for (int j = N-1; j > i; j--)
		{
			A[i][j] /= A[i][i];
		}
		b[i]/=A[i][i];
		A[i][i] = 1;
		
		for (int j = i + 1; j < N; j++)
		{
			for (int k = N-1; k > i; k--)
			{
				A[j][k] -= A[j][i] * A[i][k];
			}
			b[j]-=A[j][i]*b[i];
			A[j][i] = 0;
		}
	}

	//обратный ход
	x[N-1]=b[N-1];
	float sum = 0;
	for (int i = N - 2; i > -1; i--)
	{
		for (int j = i + 1; j < N; j++)
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
void malovki::lab3()
{
	double* alpha = new double[N];
	double* betta = new double[N];
	alpha[0] = -A[0][1]/A[0][0]; betta[0] = b[0]/A[0][0];
	for (int i = 1; i < N-2; i++)
	{
		alpha[i] = -A[i][i+1]/(A[i][i] + alpha[i-1]*A[i][i-1]);
		betta[i] = (b[i] - A[i][i-1]*betta[i-1])/(A[i][i] + alpha[i-1] * A[i][i-1]);
	}
	betta[N-1] = (b[N-1] - A[N-1][N-2]*betta[N-2]) / (A[N-1][N-1] + A[N-1][N-2]*alpha[N-2]);
	x[N-1] = betta[N-1];
	for (int i = N - 2; i >= 0; i--) 
	{
		x[i] = alpha[i]*x[i+1] + betta[i];
	}
	delete[] alpha;
	delete[] betta;
}



/**
 * Метод Холецкого
 */
void malovki::lab4()
{
	double** S = new double* [N];
	for (int i = 0; i < N; i++) {
		S[i] = new double[N];
		for(int j = 0; j < N; j++)
			S[i][j] = 0;
	}
	int* D = new int[N];
	for (int i = 0; i < N; i++)
		D[i] = 0;
	double temp;
	for(int i = 0; i < N; i++){
		temp = A[i][i];
		for (int j = 0; j < i; j++)
			temp -= D[j] * S[j][i] * S[j][i];
		D[i] = (temp > 0)? 1: -1;
		S[i][i] = sqrt(D[i] * temp);
		double nakSum;
		for (int j = i + 1; j < N; j++) {
			nakSum = 0;
			for (int k = 0; k < j; k++) 
				nakSum += D[k] * S[k][i] * S[k][j];
			S[i][j] = (A[i][j] - nakSum) / (D[i] * S[i][i]);
		}
	}
	double* y = new double[N];
	for (int i = 0; i < N; i++) {
		b[i] /= S[i][i];
		y[i] = b[i];
		for (int j = i + 1; j < N; j++)
			b[j] -= b[i] * S[i][j];
	}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			S[i][j] *= D[i];
	for (int i = N - 1; i >= 0; i--) {
		y[i] /= S[i][i];
		x[i] = y[i];
		for(int j = i - 1; j >= 0; j--)
			y[j] -= y[i] * S[j][i];
	}
	
	for (int i = 0; i < N; i++)
		delete[]S[i];
	delete[]S;
	delete[]D;
	delete[]y;
}



/**
 * Метод Якоби или Зейделя
 */
void malovki::lab5()
{
        const double eps = 1e-18;
	double sum = 0;
	double* xx = new double[N]; //вектор приближённых решений
	for (int i = 0; i < N; i++)
	{
		xx[i] = x[i];
	}
	while (true) 
	{
		label:
		for (int i = 0; i < N; i++)
		{
			sum = 0;
			for (int j = 0; j < i; j++)
				sum += A[i][j] * x[j];
			for (int j = i + 1; j < N; j++)
				sum += A[i][j] * xx[j];
			x[i] = (b[i] - sum)/A[i][i];
		}
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - xx[i]) > eps)
			{
				for (int k = 0; k < N; k++)
					xx[i] = x[i];
				goto label;
			}
		}
		break;
	}
}



/**
 * Метод минимальных невязок
 */
void malovki::lab6()
{
	const double eps = 1e-20;
	double* r = new double[N];//вектор невязок
	double rNorm = 0; //норма вектора невязок
	double* Ar = new double[N]; //произведение матрицы A на вектор невязок
	double ArNorm = 0; 
	double tau = 0; //итерационный параметр
	
	do{
		for (int i = 0; i < N; i++)
			r[i] = -b[i];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				r[i] += A[i][j] * x[j];

		for (int i = 0; i < N; i++)
			Ar[i] = 0;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				Ar[i] += A[i][j] * r[j];

		ArNorm = 0;
		for (int i = 0; i < N; i++)
			ArNorm += Ar[i] * Ar[i];

		tau = 0;
		for (int i = 0; i < N; i++)
			tau += Ar[i] * r[i];
		tau /= ArNorm;

		rNorm = 0;
		for (int i = 0; i < N; i++)
			rNorm += r[i] * r[i];

		for (int i = 0; i < N; i++)
			x[i] = x[i] - tau * r[i];

	} while (rNorm > eps);
}



/**
 * Метод сопряженных градиентов
 */
void malovki::lab7()
{
	const double eps = 1e-18;

	double* xPrev = new double[N];
	for (int i = 0; i < N; i++)
		xPrev[i] = 0;
	double* r = new double[N];
	double* rPrev = new double[N];
	double* A_xPrev = new double[N];
	double r_r, rPrev_rPrev;
	double Ax_x;
	double alpha, betta;

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

		r_r = 0;
		for (int i = 0; i < N; i++)
			r_r += r[i] * r[i];

		betta = r_r / rPrev_rPrev;

		for (int i = 0; i < N; i++)
			xPrev[i] = r[i] + betta * xPrev[i];

		for (int i = 0; i < N; i++)
			rPrev[i] = r[i];

	} while (r_r > eps);

	delete[]xPrev;
	delete[]rPrev;
	delete[]r;
	delete[]A_xPrev;
}


void malovki::lab8()
{
	int n = N;
	double** C = new double* [n];
	double eps = 1e-10, err = 0;
	int cur_i, cur_j;

	for (int i = 0; i < n; i++) {
		C[i] = new double[n];
	}

	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (i != j) {
				err += A[i][j];
			}
			C[i][j] = 0;
		}
	}

	while (err > eps) {
		//cout<<1<<endl;
		double max = abs(A[0][1]);
		//находим max в верхнем треугольнике
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				if (abs(A[i][j]) > max) {
					max = abs(A[i][j]);
					cur_i = i;
					cur_j = j;
				}
			}
		}
		double phi;
		if (A[cur_i][cur_i] != A[cur_j][cur_j]) {
			phi = 0.5 * atan(2 * A[cur_i][cur_j] / (A[cur_j][cur_j] - A[cur_i][cur_i]));
		}
		else { 
			phi = M_PI / 4; 
		}
		double s = sin(phi), c = cos(phi);

		C[cur_i][cur_j] = c * c * A[cur_i][cur_i] - 2 * s * c * A[cur_i][cur_j] + s * s * A[cur_j][cur_j];
		C[cur_j][cur_j] = s * s * A[cur_i][cur_i] + 2 * s * c * A[cur_i][cur_j] + c * c * A[cur_j][cur_j];
		C[cur_i][cur_j] = (c * c - s * s) * A[cur_i][cur_j] + s * c * (A[cur_j][cur_j] - A[cur_i][cur_i]);
		C[cur_j][cur_i] = C[cur_i][cur_j];

		for (int k = 0; k < n; k++) {
			if (k != cur_i && k != cur_j) {
				C[cur_i][k] = c * A[cur_i][k] - s * c * A[cur_j][k];
				C[k][cur_i] = C[cur_i][k];
				C[cur_j][k] = s * A[cur_i][k] + c * A[cur_j][k];
				C[k][cur_j] = C[cur_j][k];
			}

			for (int l = 0; l < n; l++)
				if (l != cur_j && k != cur_i)
					C[k][l] = A[k][l];
		}

		err = 0;
		for (int i = 0; i < n; i++)
			for (int j = i + 1; j < n; j++)
				if (i != j) err += C[i][j];


		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				A[i][j] = C[i][j];
	}
	cout << "СЗ Mat(A):" << endl;
	for (int i = 0; i < n; i++)
		cout << A[i][i] << endl;
	cout << "))";
}


void malovki::lab9()
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


std::string malovki::get_name()
{
  return "Malov K. I.";
}
