#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <windows.h>
#include "pch.h"

using namespace std;
class Worker
{
	string SecName, titul, publishing;
	int year, zp;

public:
	Worker();
	void showWorker();
	string get_SecName()
	{
		return SecName;
	}
	string get_titul()
	{
		return titul;
	}
	int get_year()
	{
		return year;
	}
	int get_zp()
	{
		return zp;
	}
};

Worker::Worker()
{
	cout << "Фамилия и инициалы:  ";
	cin >> SecName;
	cout << "Должность:  ";
	cin >> titul;
	cout << "Год поступления на работу:  ";
	cin >> year;
	cout << "Зарплата:  ";
	cin >> zp;

}
void Worker::showWorker()
{
	cout << "Фамилия и инициалы:  " << SecName << endl;
	cout << "Должность:  " << titul << endl;
	cout << "Год поступления на работу:  " << year << endl;
	cout << "Зарплата:  " << zp << endl;
}
void creature(Worker *obj, int n)
{
	for (int i = 0; i < n; i++)
	{
		cout << "=====================================\n";
		cout << "Запись №" << i + 1 << ":\n";
		obj[i].showWorker();
	}
	cout << "=====================================\n";
}
void search_zp(Worker *obj, int n, int zp)
{
	for (int i = 0; i < n; i++)
	{
		if (zp < obj[i].get_zp())
		{
			cout << "Фамилия и инициалы : ' " << obj[i].get_SecName() << " '" << endl;
		}
	}
}
void search_year(Worker *obj, int n, int year)
{
	for (int i = 0; i < n; i++)
	{
		if (year < obj[i].get_year())
			cout << "Фамилия и инициалы : ' " << obj[i].get_SecName() << " '" << endl;
	}
}
void search_titul(Worker *obj, int n, string titul)
{
	for (int i = 0; i < n; i++)
	{
		if (titul == obj[i].get_titul())
			cout << "Фамилия и инициалы : ' " << obj[i].get_SecName() << " '" << endl;
	}
}
int main()
{
	setlocale(LC_ALL, "rus");
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	int n = 0, v = 1;
	int year, zp;
	string titul;
	cout << "Кол-во записей -> ";
	cin >> n;
	cout << "Создание списка:\n";
	Worker* obj = new Worker[n];
	cout << "1. Вывести список;\n"
		"2. Список работников: стаж работы (больше числа);\n"
		"3. Список работников: зарплата (больше числа);\n"
		"4. Список работников: должность;\n"
		"0. Выход.\n";
	do
	{
		cout << " -> ";
		cin >> v;
		switch (v)
		{
		case 1:
			cout << "Список:\n";
			creature(obj, n);
			break;
		case 2:
			cout << "Список по стажу работы\n";
			cout << "Введите стаж:  ";
			cin >> year;
			search_year(obj, n, year);
			break;
		case 3:
			cout << "Список по зарплате\n";
			cout << "Введите издательство:  ";
			cin >> zp;
			search_zp(obj, n, zp);
			break;
		case 4:
			cout << "Список по даолжности\n";
			cout << "Введите год:  ";
			cin >> titul;
			search_titul(obj, n, titul);
			break;
		default:
			return 0;
		}
	} while (v);
	return 0;
}
