#include "landyshevav.h"

/**
 * Введение в дисциплину
 */
void landyshevav::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void landyshevav::lab2()
{
    double p;
    int maxn;

    for (int k = 0; k < N - 1; k++) {
        maxn = k;
        for (int i = k + 1; i < N; i++)
            if (abs(A[i][k]) > abs(A[maxn][k])) maxn = i; ///  Выбор главного элемента
        std::swap(A[maxn], A[k]);///Меняем строки местами
        std::swap(b[maxn], b[k]);

        for (int i = k + 1; i < N; i++) {
            p = A[i][k] / A[k][k];
            for (int j = k; j < N; j++)
                A[i][j] -= p * A[k][j];
            b[i] -= p * b[k];
        }
    }

    for (int i = 0; i < N; i++) {
        x[i] = b[i];
    }

    for (int i = N - 1; i >= 0; i--) {//Перебираем массив
        for (int j = i + 1; j < N; j++)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];

    }

}






/**
 * Метод прогонки
 */
void landyshevav::lab3()
{

        double* al = new double[N];
        double* bl = new double[N];
        al[0] = -A[0][1] / A[0][0]; bl[0] = b[0] / A[0][0];
        for (int i = 1; i < N - 2; i++)
        {
            al[i] = -A[i][i + 1] / (A[i][i] + al[i - 1] * A[i][i - 1]);
            bl[i] = (b[i] - A[i][i - 1] * bl[i - 1]) / (A[i][i] + al[i - 1] * A[i][i - 1]);
        }
        bl[N - 1] = (b[N - 1] - A[N - 1][N - 2] * bl[N - 2]) / (A[N - 1][N - 1] + A[N - 1][N - 2] * al[N - 2]);
        x[N - 1] = bl[N - 1];
        for (int i = N - 2; i >= 0; i--)
        {
            x[i] = al[i] * x[i + 1] + bl[i];
        }
        delete[] al;
        delete[] bl;
 
}



/**
 * Метод Холецкого
 */
void landyshevav::lab4()
{
    double* d = new double[N];
    double** S = new double* [N];
    for (int i = 0; i < N; i++) {
        S[i] = new double[N];
    }

    if (A[0][0] > 0) d[0] = 1;
    else d[0] = -1;
    S[0][0] = sqrt(abs(A[0][0]));

    for (int i = 1; i < N; i++) {
        S[0][i] = A[0][i] / (S[0][0] * d[0]);
    }


    double sum_d = 0;
    for (int i = 1; i < N; i++) {//матрица S
        for (int k = 0; k < i; k++) {
            sum_d += pow(S[k][i], 2) * d[k];
        }
        if ((A[i][i] - sum_d) > 0) {
            d[i] = 1;
        }
        else {
            d[i] = -1;
        }
        S[i][i] = sqrt(d[i] * (A[i][i] - sum_d));
        sum_d = 0;
        double sum_s = 0;
        for (int j = i + 1; j < N; j++) {
            for (int k = 0; k < j; k++) {
                sum_s += d[k] * S[k][i] * S[k][j];
            }
            S[i][j] = (A[i][j] - sum_s) / (d[i] * S[i][i]);
            sum_s = 0;
        }
    }



    double* y = new double[N];
    y[0] = b[0] / S[0][0];

    double sum_s = 0; //Решение уравнения S^t*y=b
    for (int i = 1; i < N; i++) {
        for (int j = 0; j < i; j++) {
            sum_s += S[j][i] * y[j];
        }
        y[i] = (b[i] - sum_s) / S[i][i];
        sum_s = 0;
    }



    x[N - 1] = y[N - 1] / (S[N - 1][N - 1] * d[N - 1]);//Решение уравнения (SD)*x=y

    double sum_sDx = 0;
    for (int i = N - 2; i >= 0; i--) {
        for (int k = i + 1; k < N; k++) {
            sum_sDx += S[i][k] * x[k];
        }
        x[i] = (y[i] - sum_sDx) / (S[i][i] * d[i]);
        sum_sDx = 0;
    }
}



/**
 * Метод Якоби или Зейделя
 */
void landyshevav::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void landyshevav::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void landyshevav::lab7()
{

}


void landyshevav::lab8()
{

}


void landyshevav::lab9()
{

}


std::string landyshevav::get_name()
{
  return "Landyshev A.V.";
}
