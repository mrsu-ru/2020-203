#include "garinma.h"

/**
 * Введение в дисциплину
 */
void garinma::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void garinma::lab2()
{
	double tmp;
for (int k = 0; k < N; k++) { 
  int m=k;
    for(int i=k+1;i<N;i++)
      if(abs(A[i][k]) > abs(A[m][k]))
        m=i;
  for(int i=0;i<N;i++)
  std::swap(A[k][i],A[m][i]);
  std::swap(b[k],b[m]);


tmp = A[k][k];
for (int j = 0; j < N; j++)    //Прямой ход
  A[k][j] = A[k][j] / tmp;
b[k] =b[k]/tmp;

for (int i = k + 1; i < N; i++) {
  tmp = A[i][k];
  for (int j = 0; j < N; j++) {
    A[i][j] =A[i][j]- A[k][j] * tmp;
  }
b[i] =b[i]- b[k] * tmp;
}
}

for (int k = N - 1; k > 0; k--)  //Обратный ход
{
  for (int i = k - 1; i >= 0; i--)
    {
       tmp = A[i][k];
       for (int j = 0; j < N; j++)
           A[i][j] =A[i][j]- A[k][j] * tmp;
       b[i] =b[i] - b[k] * tmp;
    }
}
for(int i=0; i<N; i++)
x[i]=b[i];
}








/**
 * Метод прогонки
 */
void garinma::lab3()
{
  int n1=N-1;
double *alfa = new double[N];
double *beta = new double[N];
double y = A[0][0];
alfa[0] = -A[0][1] / y;
beta[0] = b[0] / y;
for (int i = 1; i < n1; i++) {
y = A[i][i] + A[i][i - 1] * alfa[i - 1];
alfa[i] = -A[i][i + 1] / y;
beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
}

x[n1] = (b[n1] - A[n1][n1-1] * beta[n1-1]) / (A[n1][n1] + A[n1][n1-1] * alfa[n1-1]);
for (int i = n1-1; i >= 0; i--) {
x[i] = alfa[i] * x[i + 1] + beta[i];
}

}



/**
 * Метод Холецкого
 */
void garinma::lab4()
{
  double **S = new double*[N];
    double **ST = new double*[N];
    double *D = new double[N];
    for(int i=0; i<N; i++)
    {
       S[i] = new double[N];
       ST[i] = new double[N];
       D[i]=0;
       for(int j=0; j<N; j++)
       {
           S[i][j]=0;
       }
    }
    //вычисляем D и S
    for(int i=0; i<N; i++)
    {
        double isum = 0;
        for(int k=0; k<i; k++)
                isum += std::abs(S[k][i]) * std::abs(S[k][i]) * D[k];

        if((A[i][i] - isum) > 0)
            D[i] = 1;
        else
            D[i] = -1;

        S[i][i] = sqrt(std::abs(A[i][i] - isum));

        for(int j=i+1; j<N; j++)
        {
            double sum=0;
            for(int k=0; k<i; k++)
                sum += S[k][i] * S[k][j] * D[k];
            S[i][j] = (A[i][j] - sum) / (S[i][i] * D[i]);
        }
    }
    //транспонируем S
    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            ST[i][j] = S[j][i];
    //умножаем ST на D
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            ST[i][j] *= D[j];
    //ST*D * y = b
    double *y = new double[N];
    for(int i=0; i<N; i++)
    {
        double s=0;
        for(int j=0; j<i; j++)
            s += y[j] * ST[i][j];
        y[i]=(b[i]-s)/ST[i][i];
    }
    //S * x = y
    for(int i=N-1; i>=0; i--)
    {
        double s=0;
        for(int j=i+1; j<N; j++)
            s += x[j] * S[i][j];
        x[i]=(y[i] - s)/S[i][i];
    }

    for(int i=0; i<N; i++)
    {
        delete [] S[i];
        delete [] ST[i];
    }
    delete [] S;
    delete [] ST;
    delete [] D;
    delete [] y;
}





/**
 * Метод Якоби или Зейделя
 */
void garinma::lab5()
{

    double *new_x = new double[N],
    eps = 1.e-10;
    bool c;

  for (int i = 0; i < N; i++)
        x[i] = 1;

    do
    {
        c = false;
        for (int i = 0; i < N; i++)
        {
            new_x[i] = b[i];
            for (int j = 0; j < N; j++)
            {
                if (i == j) continue;
                new_x[i] -= A[i][j]*x[j];
            }

            new_x[i] /= A[i][i];
            if (!c) c = (fabs(new_x[i] - x[i]) > eps);
            x[i] = new_x[i];
        }

    }while(c);

    delete[] new_x;
}



/**
 * Метод минимальных невязок
 */
void garinma::lab6()
{
double *new_x = new double[N],
*r = new double[N],
eps = 1.e-10;

    for (int i = 0; i < N; i++)
        x[i] = 0;

    do
    {
        for (int i = 0; i < N; i++)
        {
            r[i] = b[i];
            for (int j = 0; j < N; j++)
                r[i] -= A[i][j] * x[j];
        }

        double tau, P = 0, Q = 0, t;
        for (int i = 0; i < N; i++)
        {
            t = 0;
            for (int j = 0; j < N; j++)
                t += A[i][j] * r[j];


            P += r[i] * t;
            Q += t * t;
        }

        tau = P / Q;
        for (int i = 0; i < N; i++)
            new_x[i] = x[i] + tau * r[i];

        double maxdif = 0;
        for (int i = 0; i < N; i++)
        {
            if (fabs(x[i] - new_x[i]) > maxdif)
				maxdif = fabs(x[i] - new_x[i]);
            x[i] = new_x[i];
        }

        if (maxdif < eps) break;
    }while(true);

    delete[] new_x;
    delete[] r;
}




/**
 * Метод сопряженных градиентов
 */
void garinma::lab7()
{
	double Del, s, sAbs;
    double eps = 1.e-10;

	double *K = new double[N];
	double *L = new double[N];
	double *M = new double[N];
	double *xrez = new double[N];//итерационные решения


	//задание начального приближения
	for (int i = 0; i<N; i++){
		xrez[i] = 0;
	}


	do {
		//нахождение скалярного произведения матрицы системы и вектора приближенного решения
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * xrez[j];
		}

		//нахождение градиента
		for (int i = 0; i < N; i++) {
			L[i] = K[i] - b[i];
		}

		//скалярное произведение матрицы системы и градиента
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * L[j];
		}


		for (int i = 0; i < N; i++) {
			M[i] = 0;
			for (int j = 0; j < N; j++) {
				M[i] += A[i][j] * K[j];
			}
		}

		s = 0;
		sAbs = 0;

		//нахождение величины смещения по направлению градиента(скалярного шага)
		for (int i = 0; i < N; i++) {
			s += K[i] * L[i];
			sAbs += M[i] * K[i];
		}
		if (s == sAbs)
			s = 1;
		else
			s = s / sAbs;
		//записываем новое приближенное решение
		for (int i = 0; i < N; i++)
			x[i] = xrez[i] - s*L[i];

		//проверка на уменьшение погрешности
		Del = abs(x[0] - xrez[0]);

		for (int i = 0; i < N; i++) {
			if (abs(x[i] - xrez[i])>Del)
				Del = abs(x[i] - xrez[i]);
				xrez[i] = x[i];
		}
	} while (eps < Del);
}


void garinma::lab8()
{
	double eps = 1.e-10;

	double **T = new double*[N];//промежуточная матрица
	for (int i = 0; i < N; i++) {
		T[i] = new double[N];
	}

//}
	while (true) {  //вычисляем норму и находим наибольший не диагональный элемент
		int im = 0;
		int jm = 0;
		double norma = 0.;
		for (int i = 0; i < N - 1; i++) {
			for (int j = i + 1; j < N; j++) {
				norma += (A[i][j] * A[i][j]);
				if (abs(A[i][j]) > abs(A[im][jm])) {
					im = i;
					jm = j;
				}
			}
		}

		if (sqrt(norma) < eps) { //след-но матрица диагональная
			break;
		}
         //значения переменных матрицы вращения
		double f = .5 * atan(2. * A[im][jm] / (A[im][im] - A[jm][jm]));
		double c = cos(f);
		double s = sin(f);

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				T[i][j] = A[i][j];
			}
		}

		for (int k = 0; k < N; k++) {  //новое приближение
			T[k][im] = A[k][im] * c + A[k][jm] * s;
			T[k][jm] = A[k][jm] * c - A[k][im] * s;
		}

		for (int k = 0; k < N; k++) {
			A[im][k] = T[im][k] * c + T[jm][k] * s;
			A[jm][k] = T[jm][k] * c - T[im][k] * s;
		}

		for (int i = 0; i < N; i++) {
			if ((i != im) && (i != jm)) {
				for (int j = 0; j < N; j++) {
					A[i][j] = T[i][j];
				}
			}
		}
	}

   // Запиываем в вектор решений вектор собственных значений матрицы
	for (int i = 0; i < N; i++) {
		x[i] = A[i][i];
	}

	for (int i = 0; i < N; i++) {
		delete[] T[i];
	}
    delete[] T;

}


void garinma::lab9()
{
double * Y = new double[N];//первый вектор приближения
    double * M = new double[N];//второй вектор приближения
    double maxL, L, sum;
    double EPS = 1e-15;

    //первичное приближение начального вектора
    for (int i = 0; i < N; i++)
        Y[i] = 0;
    Y[0] = 1;

    do{
        sum = 0;
        //нахождение скалярного произведения векторов приближения
        for (int i = 0; i < N; i++)
            sum += Y[i] * Y[i];

        L = sqrt(sum);//норма вектора приближения

        //построение последовательности векторов
        for (int i = 0; i < N; i++)
        {
            M[i] = 0;
            for (int j = 0; j < N; j++)
                M[i] += A[i][j] * Y[j] / L;
        }
        sum = 0;

        //сравнение нормы полученного вектора с заданной погрешностью
        for (int i = 0; i < N; i++)
            sum += M[i] * M[i];
        maxL = sqrt(sum);

        for (int i = 0; i<N; i++)
            Y[i] = M[i];
    } while (abs(maxL - L)>EPS);
    x[0]=maxL;
}


std::string garinma::get_name()
{
  return "Garin M.A.";
}
