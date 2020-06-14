#include "kvashninka.h"

/**
 * Введение в дисциплину
 */
void kvashninka::lab1()
{
    cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kvashninka::lab2()
{
    int *sw1 = new int[N];
    for (int i = 0; i < N; i++)
    {
        sw1[i] = i;
    }
    for (int i = 0; i < N; i++)
    {
        int Maxi = i;
        for (int j = i + 1; j < N; j++)
        {
            if (abs(A[j][i]) > abs(A[Maxi][i]))
            {
                Maxi = j;
            }
        }
        if (Maxi != i)
        {
            swap(A[Maxi], A[i]);
            swap(b[Maxi], b[i]);
            swap(sw1[Maxi], sw1[i]);
        }
        b[i] /= A[i][i];
        for (int j = N - 1; j >= i; j--)
        {
            A[i][j] /= A[i][i];
        }

        for (int j = 0; j < N; j++)
        {
            if (j != i)
            {
                b[j] -= A[j][i] * b[i];
                for (int k = N - 1; k >= i; k--)
                {
                    A[j][k] -= A[j][i] * A[i][k];
                }

            }

        }
    }
    for (int i = 0; i < N; i++)
    {
        x[i] = b[i];
    }
    for (int i = 0; i < N; i++)
    {
        if (sw1[i] != i)
        {
            swap(A[sw1[i]], A[i]);
            swap(b[sw1[i]], b[i]);
        }
    }
}



/**
 * Метод прогонки
 */
void kvashninka::lab3()
{
    double alpha[N - 1], beta[N];
    alpha[0] = - A[0][1] / A[0][0];
    beta[0] = b[0] / A[0][0];

    for (int i = 1; i < N; i++)
    {
        double y = A[i][i] + A[i][i - 1] * alpha[i - 1];
        alpha[i] = - A[i][i + 1] / y;
        beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
    }

    x[N - 1] = beta[N - 1];
    for (int i = N - 2; i >= 0; i--)
    {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
}



/**
 * Метод Холецкого
 */
void kvashninka::lab4()
{
    double** S = new double* [N];
    double* D = new double[N];

    for (int i = 0; i < N; i++)
    {
        S[i] = new double[N];
        memset(S[i], 0, sizeof(double) * N);
    }



    for (int i = 0; i < N; i++)
    {
        double sum1 = A[i][i];

        for (int k = 0; k <= i - 1; k++)
        {
            sum1 -= D[k] * S[k][i] * S[k][i];
        }

        D[i]=sum1/fabs(sum1);

        S[i][i] = sqrt(sum1);

        for (int j = i + 1; j < N; j++)
        {
            double sum2 = 0;
            for (int k = 0; k <= j - 1; k++)
            {
                sum2 += S[k][i] * S[k][j]* D[k];
            }

            S[i][j] = (A[i][j] - sum2) / (D[i] * S[i][i]);
        }
    }

    for (int i = 0; i < N; i++)
    {
        b[i] /= S[i][i];
        x[i] = b[i];

        for (int j = i + 1; j < N; j++)
        {
            b[j] -= b[i] * S[i][j];
        }
    }



    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            S[i][j] *= D[i];
        }
    }


    for (int i = N - 1; i >= 0; i--)
    {
        x[i] /= S[i][i];

        for (int j = i - 1; j >= 0; j--)
        {
            x[j] -= x[i] * S[j][i];
        }
    }
}



/**
 * Метод Якоби или Зейделя
 */
void kvashninka::lab5()
{
    double eps = 1e-35;
    double* z = new double[N];

    for (int i = 0; i < N; i++)
    {
        z[i] = b[i];
        for (int j = 0; j < N; j++)
        {
            z[i] -= A[i][j] * x[j];
        }
    }

    double norma = 0;

    for (int i = 0; i < N; i++)
    {
        norma += z[i] * z[i];
    }

    while (norma > eps)
    {
        for (int i = 0; i < N; i++)
        {
            double sum1 = 0;
            double sum2 = 0;
            for (int j = 0; j < i; j++)
            {
                sum1 += A[i][j] * x[j];
            }

            for (int j = i + 1; j < N; j++)
            {
                sum2 += A[i][j] * x[j];
            }

            x[i] = (b[i] - sum1 - sum2) / A[i][i];
        }

        for (int i = 0; i < N; i++)
        {
            z[i] = b[i];
            for (int j = 0; j < N; j++)
            {
                z[i] -= A[i][j] * x[j];
            }
        }

        norma = 0;

        for (int i = 0; i < N; i++)
        {
            norma += z[i] * z[i];
        }
    }
}



/**
 * Метод минимальных невязок
 */
void kvashninka::lab6()
{
    double eps = 1e-35;
    double *r = new double[N];

    for (int i = 0; i < N; i++)
    {
        r[i] = -b[i];
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            r[i] += A[i][j] * x[j];
        }
    }

    double rNorm = 0;
    for (int i = 0; i < N; i++)
    {
        rNorm += r[i] * r[i];
    }

    double *Ar = new double[N];
    memset(Ar, 0, sizeof(double) * N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            Ar[i] += A[i][j] * r[j];
        }
    }

    double ArNorm = 0;
    for (int i = 0; i < N; i++)
    {
        ArNorm += Ar[i] * Ar[i];
    }

    double t = 0;
    for (int i = 0; i < N; i++)
    {
        t += Ar[i] * r[i];
    }
    t /= ArNorm;

    while (rNorm > eps)
    {
        for (int i = 0; i < N; i++)
        {
            x[i] -= t * r[i];
        }

        for (int i = 0; i < N; i++)
        {
            r[i] = -b[i];
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                r[i] += A[i][j] * x[j];
            }
        }

        memset(Ar, 0, sizeof(double) * N);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Ar[i] += A[i][j] * r[j];
            }
        }

        ArNorm = 0;
        for (int i = 0; i < N; i++)
        {
            ArNorm += Ar[i] * Ar[i];
        }

        t = 0;
        for (int i = 0; i < N; i++)
        {
            t += Ar[i] * r[i];
        }
        t /= ArNorm;

        rNorm = 0;
        for (int i = 0; i < N; i++)
        {
            rNorm += r[i] * r[i];
        }
    }
}




/**
 * Метод сопряженных градиентов
 */
void kvashninka::lab7()
{

}


void kvashninka::lab8()
{

}


void kvashninka::lab9()
{

}


std::string kvashninka::get_name()
{
    return "Kvashnin K.A.";
}
