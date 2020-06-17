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
		for (int i = 0; i < N; i++) if ( abs(x[i] - prevX[i]) > err ){ err = abs( x[i] - prevX[i] ); break; };
			
		for (int i = 0; i < N; i++) prevX[i] = x[i];
	}
	while(err > 1e-20);

}



/**
 * Метод сопряженных градиентов
 */
void guskovas::lab7()
{

	double* prevX = new double[N];
	double* prevR = new double[N];
	double* r = new double[N];
	double* z = new double[N];

	for (int i = 0; i < N; i++) {
		r[i] = b[i];
		z[i] = b[i];
	}

	double err = 0;
	double a = 0;
	double det = 0;
	double beta = 0;
	double Az = 0;

	do {

		for (int i = 0; i < N; i++) {
			prevR[i] = r[i];
			prevX[i] = x[i];
		}

		a = 0; det = 0;

		for (int i = 0; i < N; i++) {
			Az = 0;
			
			for (int j = 0; j < N; j++) Az += (A[i][j] * z[j]);
			
			a += (prevR[i] * prevR[i]); det += (Az * z[i]);
		}
		a /= det;

		for (int i = 0; i < N; i++)	x[i] = prevX[i] + (a * z[i]);
		

		err = abs(x[0] - prevX[0]);

		for (int i = 1; i < N; i++)
			if (abs(x[i] - prevX[i]) > err)
				err = abs(x[i] - prevX[i]);

		for (int i = 0; i < N; i++) {
			Az = 0;

			for (int j = 0; j < N; j++) Az += (A[i][j] * z[j]);

			r[i] = prevR[i] - (a * Az);
		}

		beta = 0; det = 0;

		for (int i = 0; i < N; i++) {
			beta += (r[i] * r[i]);
			det += (prevR[i] * prevR[i]);
		}
		beta /= det;

		for (int i = 0; i < N; i++) z[i] = r[i] + (beta * z[i]);

	}while( !(err < 1e-20) );

}


void guskovas::lab8()
{
  double err = 0;
  double** C = new double*[N];

  for (int i = 0; i < N; i++){
  	for (int j = i+1; j < N; j++) err += (i != j)? A[i][j] * A[i][j] : 0;
  	C[i] = new double[N];
   	for (int j= i+1; j < N; j++) C[i][j] = 0; // вычисляем ошибку errc
   }
	  
  while (sqrt(err) > 1e-20){

	int mI = 0, mJ = 1;

  	for (int i = 0; i < N; i++)
      for (int j = i+1; j<N; j++)
  	    if ( abs(A[i][j]) > abs(A[mI][mJ]) ){
  	    	mI = i; mJ = j;	// arg max |Aij|
		}


	double a = (A[mI][mI] == A[mJ][mJ])? // (11)
					M_PI/4 : 
					a = 0.5 * ( atan( 2*A[mI][mJ] / (A[mI][mI]-A[mJ][mJ]) ) );
	
	double c = cos(a);
	double s = sin(a);
 
 	// Найдем значения элементов матрицы С (10)
	C[mI][mI] = pow(c, 2)*A[mI][mI] - 2*s*c*A[mI][mJ] + pow(s, 2)*A[mJ][mJ];
	C[mJ][mJ] = pow(s, 2)*A[mI][mI] + 2*s*c*A[mI][mJ] + pow(c, 2)*A[mJ][mJ];
	C[mI][mJ] = (pow(c, 2) - pow(s, 2))*A[mI][mJ] + s*c*(A[mJ][mJ] - A[mI][mI]);
	C[mJ][mI] = C[mI][mJ];
	
	for (int k = 0; k < N; k++){
	  if (k != mI  &&  k != mJ){
	  	C[mI][k] = c*A[mI][k] - s*c*A[mJ][k];
	  	C[k][mI] = C[mI][k];
	  	C[mJ][k] = s*A[mI][k] + c*A[mJ][k];
	  	C[k][mJ] = C[mJ][k];
	  } 
	  for (int l = 0; l < N; l++)
	    if (k != mI && k != mJ && l != mI && l != mJ) C[k][l] = A[k][l];	  
	}
	
	err = 0;
	for (int i = 0; i < N; i++)
      for (int j = i+1; j < N; j++)
  	    if (i != j) err += C[i][j] * C[i][j]; // Вычисляем errc

	for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) A[i][j] = C[i][j]; // A := C

  }
  
  for (int i = 0; i < N; i++) x[i] = A[i][i];

}


void guskovas::lab9()
{
	double *xk = new double[N]; //Yk
	double L = 0; // лямбда
	double newL = 0; //следующее приближение
    x[0] = 1; // x = y0,   a1 != 0     (a1=1)

    do
    {
        newL = 0;
        for (int i = 0; i < N; i++)
        {
            xk[i] = 0;
            for (int j = 0; j < N; j++) xk[i] += A[i][j] * x[j]; // Yk+1 = Ak*Yk = Ak*x (2)
            newL += x[i] * xk[i]; // 
        }

        if (fabs(newL - L) < 1e-20) break; /////////////////////BREAK

        L = newL;
        double n = 0;
        for (int i = 0; i < N; i++) n += xk[i] * xk[i];
        n = sqrt(n);

        for (int i = 0; i < N; i++) x[i] = xk[i] / n; // yk+1=yk/||yk||

    } while (1);

    cout<<"RESULT: "<< L << endl;

}


std::string guskovas::get_name()
{
  return "Guskov A.S.";
}
