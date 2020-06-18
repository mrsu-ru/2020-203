#include "zevaykinae.h"

/**
 * Введение в дисциплину
 */
void zevaykinae::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void zevaykinae::lab2()
{

	for (int i = 0; i < N; i++) 
		x[i] = b[i]; 
	long double m; 
	for (int k = 0; k < N - 1; k++) { 
		for (int i = k + 1; i < N; i++) { 
			m = A[i][k] / A[k][k]; 	
			for (int j = k; j < N; j++) { 
				A[i][j] = A[i][j] - m * A[k][j]; 
			} 
			x[i] = x[i] - m * x[k]; 
		} 
	} 
	for (int i = N - 1; i >= 0; i--) { 
		for (int j = i + 1; j < N; j++) 
		x[i] = x[i] - A[i][j] * x[j]; 
		x[i] = x[i] / A[i][i]; 
	} 
}



/**
 * Метод прогонки
 */
void zevaykinae::lab3()
{
	double *alpha = new double [N];
	double *beta = new double [N];

	alpha[0] = -A[0][1]/A[0][0];
	beta[0] = b[0]/A[0][0];

	for(int i=1; i<N; i++) //здесь определяются прогоночные коэффициенты
	{
		alpha[i] = A[i][i+1]/(-A[i][i] - A[i][i-1]*alpha[i-1]);
		beta[i] = (-b[i] + A[i][i-1]*beta[i-1])/(-A[i][i] - A[i][i-1]*alpha[i-1]);
	}

	x[N-1] = beta[N-1];
	for(int i=N-2; i>=0; i--) //решение
		x[i] = alpha[i]*x[i+1] + beta[i];

	delete [] alpha;
	delete [] beta;
}



/**
 * Метод Холецкого
 */
void zevaykinae::lab4()
{
	double** S = new double* [N];
	for (int i = 0; i < N; i++) {
		S[i] = new double[N];
		for(int j = 0; j < N; j++)
			S[i][j] = 0;
	}
	int* D = new int[N];
	for (int i = 0; i < N; i++)
		D[i] = 0;
	double temp;
	for(int i = 0; i < N; i++){
		temp = A[i][i];
		for (int j = 0; j < i; j++)
			temp -= D[j] * S[j][i] * S[j][i];
		D[i] = (temp > 0)? 1: -1;
		S[i][i] = sqrt(D[i] * temp);
		double nakSum;
		for (int j = i + 1; j < N; j++) {
			nakSum = 0;
			for (int k = 0; k < j; k++) 
				nakSum += D[k] * S[k][i] * S[k][j];
			S[i][j] = (A[i][j] - nakSum) / (D[i] * S[i][i]);
		}
	}
	double* y = new double[N];
	for (int i = 0; i < N; i++) {
		b[i] /= S[i][i];
		y[i] = b[i];
		for (int j = i + 1; j < N; j++)
			b[j] -= b[i] * S[i][j];
	}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			S[i][j] *= D[i];
	for (int i = N - 1; i >= 0; i--) {
		y[i] /= S[i][i];
		x[i] = y[i];
		for(int j = i - 1; j >= 0; j--)
			y[j] -= y[i] * S[j][i];
	}
	
	for (int i = 0; i < N; i++)
		delete[]S[i];
	delete[]S;
	delete[]D;
	delete[]y;
}



/**
 * Метод Якоби
 */
void zevaykinae::lab5()
{
	double eps = 1e-14;
	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double *p_x = new double[N];
	double norma = 0;
	do {
		for (int i = 0; i < N; i++) {
			p_x[i] = x[i];
		}
		for (int i = 0; i < N; i++) {
			double sum = 0;
			for (int j = 0; j < N; j++) {
				if (i != j) {
					sum += (A[i][j] * p_x[j]);
				}
			}

			x[i] = (b[i] - sum) / A[i][i];
		}
		norma = 0;
		for (int i = 0; i < N; i++) {
			if ((p_x[i] - x[i]) > norma) {
				norma = abs(p_x[i] - x[i]);
			}
		}
	} while (norma > eps);

	delete[] p_x;
}



/**
 * Метод минимальных невязок
 */
void zevaykinae::lab6()
{
 	double eps = 1e-15;

	// Задаем вектор значений x на предыдущий итерации
	double *px = new double[N];
	for (int i = 0; i < N; i++) {
		px[i] = 0.0;
	}

	// Задаем вектор невязок
	double *r = new double[N];

	int iteration = 0;
	while (true) {
		iteration++;

	// Рассчитываем вектор невязки
		for (int i = 0; i < N; i++) {
			r[i] = b[i];

			for (int j = 0; j < N; j++) {
				r[i] -= (A[i][j] * px[j]);
			}
		}

	// Рассчитываем итерационный параметр
		double t = 0.0;
		double temp = 0.0;
		for (int i = 0; i < N; i++) {
			double Ar = 0.0;

			for (int j = 0; j < N; j++) {
				Ar += (A[i][j] * r[j]);
			}

			t += (Ar * r[i]);
			temp += (Ar * Ar);
		}
		t /= temp;

	// Рассчитывается новое приближение к вектору x
		for (int i = 0; i < N; i++) {
			x[i] = px[i] + t * r[i];
		}

	// Посчитаем максимальную по модулю погрешность
		double err = 0.0;
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - px[i]) > err) {
				err = abs(x[i] - px[i]);
			}
		}

	// При достижении необходимой точности завершаем процесс
		if (err < eps) {
			break;
		}

	// Текущее зачение итерации представим как предыдущее
		for (int i = 0; i < N; i++) {
			px[i] = x[i];
		}
	}

	cout << "Number of iterations : " << iteration << endl;

	delete[] px;
	delete[] r;
}



/**
 * Метод сопряженных градиентов
 */
void zevaykinae::lab7()
{
	for (int i=0; i<N; i++)
           	   x[i]=0;
      	 double *r=new double[N];
       	for (int i=0; i<N; i++){
       	     r[i]=b[i];
        	 for (int j=0; j<N; j++)
         	   r[i]-=A[i][j]*x[j];
       	}
       	double *z=new double[N];
        	 for (int i=0; i<N; i++)
         	   z[i]=r[i];
       	double eps=10e-16;
       	double var=0;
       	double alpha = 0;
       	for(;;){
        	    double differ=0, sum1=0, sum2=0, vec=0;
          	   for (int i=0; i<N; i++){
              	  	vec=0;
               	 for (int k=0; k<N; k++)
               	      vec+=A[i][k]*z[k];
               	 sum1+=r[i]*r[i];
               	 sum2+=vec*z[i];
           	 }
            	alpha=sum1/sum2;
            	for (int i=0; i<N; i++){
           	         var=x[i];
            		x[i]+=alpha*z[i];
             	   	if(abs(x[i]-var)>differ)
               	     		differ=abs(x[i]-var);
             	}
           	if(differ<eps) break;
            	sum2=0;
            	sum1=0;

             	for (int i=0; i<N; i++){
              	    	vec=0;
               	 	sum2+=r[i]*r[i];
                	for (int j=0; j<N; j++)
                 		vec+=A[i][j]*z[j];
                	r[i]-=alpha*vec;
                 	sum1+=r[i]*r[i];
             	 }

           	 for (int i=0; i<N; i++){
              	 	sum1+=r[i]*r[i];
                 	sum2+=vec*z[i];
           	 }

            for (int i=0; i<N; i++)
           	 z[i]=r[i]+sum1*z[i]/sum2;
     	}
      	delete[] r;
      	delete[] z;
}


void zevaykinae::lab8()
{
     double **H = new double*[N], eps = 1.e-10;
     for (int i = 0; i < N; i++) H[i] = new double[N];

    do
    {
        double n = 0;
        int i_max = 0, j_max = 1;
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++)
            {
                if (fabs(A[i][j]) > fabs(A[i_max][j_max]))
                {
                    i_max = i;
                    j_max = j;
                }

                n += A[i][j] * A[i][j];
            }

        if (sqrt(n) < eps) break;

        double fi = 0.5 * atan(2 * A[i_max][j_max] / (A[i_max][i_max] - A[j_max][j_max]));
        for (int i = 0; i < N; i++)
        {
            H[i][i_max] = A[i][i_max] * cos(fi) + A[i][j_max] * sin(fi);
            H[i][j_max] = A[i][j_max] * cos(fi) - A[i][i_max] * sin(fi);
        }

        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                if (j != i_max && j != j_max) H[i][j] = A[i][j];

        for (int j = 0; j < N; j++)
        {
            A[i_max][j] = H[i_max][j] * cos(fi) + H[j_max][j] * sin(fi);
            A[j_max][j] = H[j_max][j] * cos(fi) - H[i_max][j] * sin(fi);
        }

        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                if (i != i_max && i != j_max) A[i][j] = H[i][j];

    }while(true);

    for (int i = 0; i < N; i++) x[i] = A[i][i];

    for (int i = 0; i < N; i++) delete[] H[i];
    delete[] H;
}


void zevaykinae::lab9()
{
        for (int i=0; i<N; i++)
            x[i]=0;
            x[0]=1;
       double *y=new double[N];
       double eps=1e-15;
       double prev_l;
       double l = 0;
       double sum;

       for(;;){
            sum = 0;
            prev_l = l;
            l = 0;
            for (int i = 0; i<N; i++){
                y[i] = 0;
                for (int k=0; k<N; k++)
                     y[i] += A[i][k]*x[k];
                l += y[i]*x[i];
                sum += y[i]*y[i];
            }
           sum = pow(sum,0.5);
             for (int i=0; i<N; i++)
                x[i] = y[i]/sum;
             if(abs(l - prev_l)<eps) break;
            }
            x[0] = l;
        delete[] y;
}


std::string zevaykinae::get_name()
{
  return "Zevaykin A.E.";
}
