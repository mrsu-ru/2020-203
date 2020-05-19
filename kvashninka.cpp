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

}



/**
 * Метод простых итераций
 */
void kvashninka::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void kvashninka::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void kvashninka::lab6()
{

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
