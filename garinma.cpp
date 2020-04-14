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

}



/**
 * Метод сопряженных градиентов
 */
void garinma::lab7()
{

}


void garinma::lab8()
{

}


void garinma::lab9()
{

}


std::string garinma::get_name()
{
  return "Garin M.A.";
}
