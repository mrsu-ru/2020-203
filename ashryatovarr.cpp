#include "ashryatovarr.h"

/**
 * Введение в дисциплину
 */
void ashryatovarr::lab1()
{
    cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void ashryatovarr::lab2()
{

    for (int k = 0; k < N; k++)
    {
        int idmax = 0;
        for (int i=0; i<N; i++)
        {
            if(abs(A[i][k]) > abs(A[idmax][k]))
                idmax = i;
        }

        for (int i = 0; i < N; i++)
        {
            swap(A[idmax][i], A[k][i]);
        }
        swap(b[idmax], b[k]);

        for (int i = k + 1; i < N; i++)
        {
            double tmp = A[i][k]/A[k][k];
            for (int j = k; j < N; j++)
            {
                A[i][j] -= A[k][j]*tmp;
            }
            b[i] -= b[k]*tmp;
        }
    }

    for(int i = 0; i<N; i++)
    {
        x[i]=b[i];
    }

    for (int i = N-1; i >= 0; i--)
    {
        for (int j = i+1; j < N; j++)
        {
            x[i] -= A[i][j]*x[j];
        }
        x[i] /= A[i][i];
    }

}



/**
 * Метод прогонки
 */
void ashryatovarr::lab3()
{

}



/**
 * Метод простых итераций
 */
void ashryatovarr::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void ashryatovarr::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void ashryatovarr::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void ashryatovarr::lab7()
{

}


void ashryatovarr::lab8()
{

}


void ashryatovarr::lab9()
{

}


std::string ashryatovarr::get_name()
{
    return "Ashryatova R.R.";
}
