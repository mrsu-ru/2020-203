#include "gorbunovaa.h"

/**
 * Ââåäåíèå â äèñöèïëèíó
 */
void gorbunovaa::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Ìåòîä Ãàóññà ñ âûáîðîì ãëàâíîãî ýëåìåíòà
 */
void gorbunovaa::lab2()
{
	setlocale(LC_ALL, "Russian");
	double* x;
	int i, j, m, k;
	double aa, bb;
	// Для матрицы нужена переменная 'a'
	x = (double *)malloc(n * sizeof(double));

	for (k = 0; k < n; k++) //Поиск максимального элемента в первом столбце
	{
		aa = abs(a[k][k]);
		i = k;
		for (m = k + 1; m < n; m++)
			if (abs(a[m][k]) > aa)
			{
				i = m;
				aa = abs(a[m][k]);
			}

		if (aa == 0)   //проверка на нулевой элемент
		{
			cout << "Система не имеет решений" << endl;
		}

		if (i != k)  //  перестановка i-ой строки, содержащей главный элемент k-ой строки
		{
			for (j = k; j < n + 1; j++)
			{
				bb = a[k][j];
				a[k][j] = a[i][j];
				a[i][j] = bb;
			}
		}
		aa = a[k][k];//преобразование k-ой строки (Вычисление масштабирующих множителей)
		a[k][k] = 1;
		for (j = k + 1; j < n; j++)
			a[k][j] = a[k][j] / aa;
		for (i = k + 1; i < n; i++)//преобразование строк с помощью k-ой строки
		{
			bb = a[i][k];
			a[i][k] = 0;
			if (bb != 0)
				for (j = k + 1; j < n + 1; j++)
					a[i][j] = a[i][j] - bb * a[k][j];
		}
	}

	for (i = n; i > 0; i--)   //Нахождение решений СЛАУ
	{
		x[i] = 0;
		aa = a[i][n + 1];
		for (j = n; j > i; j--)
			aa = aa - a[i][j] * x[j];
		x[i] = aa;
	}

	cout << "Решение системы:" << endl;  //вывод решений
	for (i = 1; i < n + 1; i++)
	{
		cout << "x[" << i << "]=" << x[i];
		cout << endl;
	}
	system("PAUSE");
}



/**
 * Ìåòîä ïðîãîíêè
 */
void gorbunovaa::lab3()
{

}



/**
 * Ìåòîä ïðîñòûõ èòåðàöèé
 */
void gorbunovaa::lab4()
{

}



/**
 * Ìåòîä ßêîáè èëè Çåéäåëÿ
 */
void gorbunovaa::lab5()
{

}



/**
 * Ìåòîä ìèíèìàëüíûõ íåâÿçîê
 */
void gorbunovaa::lab6()
{

}



/**
 * Ìåòîä ñîïðÿæåííûõ ãðàäèåíòîâ
 */
void gorbunovaa::lab7()
{

}


void gorbunovaa::lab8()
{

}


void gorbunovaa::lab9()
{

}


std::string gorbunovaa::get_name()
{
  return "Gorbunov A. A.";
}