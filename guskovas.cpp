#include "guskovas.h"

/**
 * Введение в дисциплину
 */
void guskovas::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void guskovas::lab2()//int n, int m, double e, double** arr, double* x) {
{

	for(int i = 0; i < N; i++){
		int indexMax = i; 
		for (int j = i + 1; j < N; j++) { 
			if (abs(A[j][i]) > abs(A[indexMax][i])) indexMax = j; 
		} 
		if (indexMax != i) { 
			for (int k = 0; k < N; k++) { 
				swap(A[i][k], A[indexMax][k]); 
			} 
		} 
	}

	//ТУДА
	int n = N;

	for(int k = 0; k < n-1; k++){
		for(int i = k+1; i < n; i++){
			double c = A[i][k] / A[k][k];
			//A[i][k] = 0;
			for(int j = k+1; j < n; j++){
				A[i][j] -= A[k][j] * c; 
			}
			b[i] -= b[k] * c;
		}
	}

	for (int i = 0; i < n; i++){
			x[i] = b[i];
		}

	//Обратно

	for(int i = n-1; i >= 0; i--){
		double S = 0;
		for(int j = i+1; j < n; j++){
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
	}
	
}



/**
 * Метод прогонки
 */
void guskovas::lab3()//N, A, b, x
{
	double *ALFA = new double[N];
	double *BETA = new double[N];

	// ТУДА
	//	находим ALFA и BETA
	ALFA[0] = -A[0][1] / A[0][0];
	BETA[0] = b[0] / A[0][0];
	for (int i = 1; i < N; i++) {
		ALFA[i] = -A[i][i + 1] 
					/ 
					(
						A[i][i] + 
						A[i][i - 1] * ALFA[i - 1]
					);
		BETA[i] =   (b[i] - A[i][i - 1] * BETA[i - 1]) 
					/ 
					(A[i][i] + A[i][i - 1] * ALFA[i - 1]);
	}

	// Обратная прогонка
	//после общего вида формул подставляем иксы 
	x[N - 1] = BETA[N - 1];
	for (int i = N - 2; i >= 0; i--) {
		x[i] = ALFA[i] * x[i + 1] + BETA[i];
	}
}



/**
 * Метод Холецкого
 */
void guskovas::lab4()
{
	double **S = new double* [N];  
		for (int i = 0; i < N; ++i) S[i] = new double[N];
		S[0][0] = sqrt( fabs(A[0][0]) );
    double n; //
    double *y = new double[N];
    double *D = new double[N];
    	A[0][0] > 0 ?
		D[0] = 1
		:
		D[0] = -1;


    for (int i = 1; i < N; ++i) S[0][i] = A[0][i] / ( D[0]*S[0][0] );
  
	    for (int i = 1; i < N; ++i){
	    	n = 0;
	    	for (int j=0; j<i; j++) n += D[j] * S[j][i] * S[j][i];

		    A[i][i] - n >= 0 ?
				D[i] = 1
		  		:
				D[i] = -1;

		S[i][i] = sqrt( D[i]*(A[i][i] - n) );
		  
	    for (int j= i+1; j < N; ++j) {
	    	double l = 0;
		    for (int k = 0; k < j; ++k) l += D[k] * S[k][i] * S[k][j];

		    S[i][j] = 
		    ( A[i][j]-l )
		    / 
		    ( D[i]*S[i][i] );
		}

	}
		  
	
	for (int i = 0; i < N; ++i){
		n=0;
		if(i == 0){
			y[0] = b[0] / S[0][0];
			++i;
		}

		for (int j=0; j < i; ++j){
			n += y[j]*S[j][i];
		}
		y[i] = (b[i] - n) / S[i][i];
	}

	x[N-1] = y[N-1] / (D[N-1] * S[N-1][N-1]);
	 
	for (int i = N-2; i >= 0; --i){
		n = 0;
		for (int j = i+1; j < N; ++j) n += x[j] * D[j] * S[i][j];
		
		x[i]=
		(y[i] - n)
		/
		( D[i]*S[i][i] );
	}	  
}



/**
 * Метод Якоби или Зейделя
 */
void guskovas::lab5()//якоби
{
	double *f = new double[N];

	double norma = 0;//error
	do {
		for (int i = 0; i < N; i++) f[i] = x[i];

		for (int i = 0; i < N; i++) {
			double result = b[i];

			for (int j = 0; j < N; j++) result = i != j ? 
					result - (A[i][j] * f[j]) 
					:
					result
			;

			x[i] = result / A[i][i];
		}

		norma = 0;
		for (int i = 0; i < N; i++) norma = abs(f[i] - x[i]) > norma ? 
			abs(f[i] - x[i]) 
			: 
			norma
		;

	} while (norma > 1e-20);

}



/**
 * Метод минимальных невязок
 */
void guskovas::lab6()
{

	double *prevX = new double[N];
	double *y = new double[N];
	double tau = 0.0, err = 0.0, Ay = 0.0, denom = 0.0;

	do{

		for (int i = 0; i < N; i++) {
			y[i] = b[i];
			for (int j = 0; j < N; j++) y[i] -= A[i][j] * prevX[j]; //y^n
		}

		tau = 0.0; denom = 0.0;

		for (int i = 0; i < N; i++) {
			Ay = 0.0;

			for (int j = 0; j < N; j++) Ay += A[i][j] * y[j];

			tau += Ay * y[i]; denom += Ay * Ay;
		}
		tau /= denom; //t^n

		for (int i = 0; i < N; i++) x[i] = prevX[i] + (tau * y[i]);
		
		err = 0.0;
		for (int i = 0; i < N; i++) if (abs(x[i] - prevX[i]) > err) err = abs(x[i] - prevX[i]);
			
		for (int i = 0; i < N; i++) prevX[i] = x[i];
	}
	while(err > 1e-20);

}



/**
 * Метод сопряженных градиентов
 */
void guskovas::lab7()
{

}


void guskovas::lab8()
{

}


void guskovas::lab9()
{

}


std::string guskovas::get_name()
{
  return "Guskov A.S.";
}
