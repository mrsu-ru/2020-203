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
    double eps = 1e-20;
    for (int i = 0; i < N; i++) {
        x[i] = 0;
    }

    double* prev_x = new double[N];

    double norma = 0;
    do {
        for (int i = 0; i < N; i++) {
            prev_x[i] = x[i];
        }

        for (int i = 0; i < N; i++) {
            double result = b[i];
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    result -= (A[i][j] * prev_x[j]);
                }
            }

            x[i] = result / A[i][i];
        }

        norma = 0;
        for (int i = 0; i < N; i++) {
            if (abs(prev_x[i] - x[i]) > norma) {
                norma = abs(prev_x[i] - x[i]);
            }
        }
    } while (norma > eps);

    delete[] prev_x;
}



/**
 * Метод минимальных невязок
 */
void landyshevav::lab6()
{
    double eps = 1e-15;

    double* prevX = new double[N];
    double* r = new double[N];

    while (true) {

        for (int i = 0; i < N; i++)
            prevX[i] = x[i];

        for (int i = 0; i < N; i++) {
            r[i] = b[i];

            for (int j = 0; j < N; j++) {
                r[i] -= A[i][j] * x[j];
            }
        }

        double tau = 0;
        double denomTau = 0;

        for (int i = 0; i < N; i++) {
            double Ar = 0;

            for (int j = 0; j < N; j++) {
                Ar += A[i][j] * r[j];
            }

            tau += Ar * r[i];
            denomTau += Ar * Ar;
        }

        tau /= denomTau;

        for (int i = 0; i < N; i++) {
            x[i] = prevX[i] + tau * r[i];
        }

        double maxErr = abs(x[0] - prevX[0]);
        for (int i = 1; i < N; i++)
            if (abs(x[i] - prevX[i]) > maxErr)
                maxErr = abs(x[i] - prevX[i]);

        if (maxErr < eps)
            break;

    }

    delete[] prevX;
    delete[] r;
}



/**
 * Метод сопряженных градиентов
 */
void landyshevav::lab7()
{
    double Del, s, sAbs;
    double eps = 1.e-10;

    double* K = new double[N];
    double* L = new double[N];
    double* M = new double[N];
    double* xrez = new double[N];


    for (int i = 0; i < N; i++) {
        xrez[i] = 0;
    }


    do {
       
        for (int i = 0; i < N; i++) {
            K[i] = 0;
            for (int j = 0; j < N; j++)
                K[i] += A[i][j] * xrez[j];
        }

        for (int i = 0; i < N; i++) {
            L[i] = K[i] - b[i];
        }


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

        for (int i = 0; i < N; i++) {
            s += K[i] * L[i];
            sAbs += M[i] * K[i];
        }
        if (s == sAbs)
            s = 1;
        else
            s = s / sAbs;

        for (int i = 0; i < N; i++)
            x[i] = xrez[i] - s * L[i];


        Del = abs(x[0] - xrez[0]);

        for (int i = 0; i < N; i++) {
            if (abs(x[i] - xrez[i]) > Del)
                Del = abs(x[i] - xrez[i]);
            xrez[i] = x[i];
        }
    } while (eps < Del);
}


void landyshevav::lab8()
{
    double eps = 1.e-10;
    double Ja[N][N];
    while (true)
    {
        int i_max = 0;
        int j_max = 0;
        for (int i = 0; i < N - 1; i++)
        {
            for (int j = i + 1; j < N; j++)
            {
                if (fabs(A[i][j]) > fabs(A[i_max][j_max]))
                {
                    i_max = i;
                    j_max = j;
                }
            }
        }

        if (fabs(A[i_max][j_max]) < eps)
        {
            break;
        }

        double fi = 0.5 * atan(2 * A[i_max][j_max] / (A[i_max][i_max] - A[j_max][j_max]));
        double cos_fi = cos(fi);
        double sin_fi = sin(fi);

        for (int k = 0; k < N; k++)
        {
            Ja[k][i_max] = A[k][i_max] * cos_fi + A[k][j_max] * sin_fi;
            Ja[k][j_max] = A[k][j_max] * cos_fi - A[k][i_max] * sin_fi;
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if ((i != i_max) && (j != j_max))
                {
                    Ja[i][j] = A[i][j];
                }
            }
        }

        for (int k = 0; k < N; k++)
        {
            A[i_max][k] = Ja[i_max][k] * cos_fi + Ja[j_max][k] * sin_fi;
            A[j_max][k] = Ja[j_max][k] * cos_fi - Ja[i_max][k] * sin_fi;
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if ((i != i_max) && (i != j_max))
                {
                    A[i][j] = Ja[i][j];
                }
            }
        }
    }

    for (int i = 0; i < N; i++)
    {
        x[i] = A[i][i];
    }
}



void landyshevav::lab9()
    {

        double eps = 1.0e-3;
        double* yk, * yk_1;
        double del_k1 = 0, del_k = 1;
        int n = N;

        yk_1 = new double[n];
        yk = new double[n];

        for (int i = 0; i < n; i++)
            yk_1[i] = 1;

        while (fabs(del_k - del_k1) > eps) {
            del_k1 = del_k;

            for (int i = 0; i < n; i++) {
                double s = 0;
                for (int j = 0; j < n; j++) {
                    s += A[i][j] * yk_1[j];
                }
                yk[i] = s;
            }

            for (int i = 0; i < n; i++) {
                if (yk[0] != 0 && yk_1[0] != 0) {
                    del_k = yk[i] / yk_1[i];
                }
            }

            for (int i = 0; i < n; i++)
                yk_1[i] = yk[i];

        }
        cout << "Answer by LabWork 9: " << del_k << endl;

    }



std::string landyshevav::get_name()
{
  return "Landyshev A.V.";
}
