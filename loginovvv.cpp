#include "loginovvv.h"

/**
 * Введение в дисциплину
 */
void loginovvv::lab1()
{
	cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void loginovvv::lab2()
{
	//cout << "hello world!" << endl;

	/////приводим к треугольнику
	int z;
	for (int i = 0; i < N; i++)
	{
		z = i;
		//ишим наибольший по модулю 
		for (int m = i+1; m < N; m++)
			if (abs(A[m][i]) > abs(A[z][i]))
				z = m;
		if (A[z][i] == 0)
			return;
		//меняем строки 
		double* vrem;
    	double v;
		vrem = A[i];
		A[i] = A[z];
		A[z] = vrem;
    	v = b[i];
    	b[i] = b[z];
    	b[z] = v;
		//делим строку , чтоб получить единичку
    	b[i]/=A[i][i];
		for (int j = N-1; j > i; A[i][j--] /= A[i][i]);
		A[i][i] = 1;
		//вычитаем строку , чтоб получить нолики
		for (int j = i + 1; j < N; j++)
		{
            b[j]-=b[i]*A[j][i];
			for (int k = N-1; k > i; k--)
				A[j][k] -= A[i][k] * A[j][i];
			A[j][i] = 0;
		}
		
	}
	//находим решение 
	x[N - 1] = b[N-1];
	for (int i = N - 2; i >= 0; i--)
	{
		x[i] = b[i];
		for (int j = i + 1; j < N; j++)
			x[i] -= A[i][j] * x[j];
	}

}



/**
 * Метод прогонки
 */
void loginovvv::lab3()
{
	double *a = new double[N];
	double *v = new double[N];
	a[1]=-A[0][1]/A[0][0];
	v[1]=b[0]/A[0][0];
	for(int i=2; i<N; i++)
	{
		a[i]=-A[i-1][i]/(A[i-1][i-2]*a[i-1]+A[i-1][i-1]);
		v[i]=(b[i-1]-A[i-1][i-2]*v[i-1])/(A[i-1][i-2]*a[i-1]+A[i-1][i-1]);
	}
	x[N-1]=(b[N-1]-A[N-1][N-2]*v[N-1])/(A[N-1][N-1]+A[N-1][N-2]*a[N-1]);
	for(int i=N-2; i>-1; i--)
	{
		x[i]=a[i+1]*x[i+1]+v[i+1];

	}

}



/**
 * Метод Холецкого
 */
void loginovvv::lab4()
{
	double **S = new double*[N];
    for (int i = 0;  i < N; i++) {
        S[i] = new double[N];
		for(int j=0; j<N; j++)
			S[i][j]=0.0;
    }
	double *d = new double[N];


	for(int i=0; i<N; i++)
	{
		double temp=A[i][i];
		for(int e=0; e<=i-1; e++){
			 temp-=d[e]*S[e][i]*S[e][i];
		}
		d[i] = signbit(temp) == false ? 1 : -1;
		S[i][i]=sqrt(temp*d[i]);
	

		for(int j=i+1; j<N; j++)
		{
			double t=0;
			for(int k=0; k<=j-1; k++)
				t+=d[k]*S[k][i]*S[k][j];
			 S[i][j] = (A[i][j]-t)/(d[i]*S[i][i]);
		}
	}
	
	
	for(int i=0; i<N; i++){
		b[i]/=S[i][i];
		for(int j=i+1; j<N; j++)
			b[j]-=b[i]*S[i][j];
	}

	for(int i=0; i<N; i++)
		for(int j=i; j<N; j++)
			S[i][j]*=d[i];
	
	for(int i=N-1; i>=0; i--){
		b[i]/=S[i][i];
		for(int j=i-1; j>=0; j--)
			b[j]-=b[i]*S[j][i];
	}

	x=b;
}



/**
 * Метод Якоби или Зейделя
 */
void loginovvv::lab5()
{
	const double eps = 1e-20;
    double *f = new double[N];
	double norma=0;

    do{
		for (int i = 0; i < N; i++) 
			f[i] = x[i];

    	for (int i = 0; i < N; i++)
    	{
      		double sum1 = 0;
      		double sum2 = 0;
      		for (int j = 0; j < i; j++)
        		sum1 += A[i][j] * x[j];

      		for (int j = i + 1; j < N; j++)
        		sum2 += A[i][j] * x[j];
      
      		x[i] = (b[i] - sum1 - sum2) / A[i][i];
    	}

		norma = 0;
    	for (int i = 0; i < N; i++) 
			if (abs(f[i] - x[i]) > norma)
				norma = abs(f[i] - x[i]);
			
  } while (norma > eps);

  delete []f;
}



/**
 * Метод минимальных невязок
 */
void loginovvv::lab6()
{
	const double eps = 1e-20;
	double *r = new double[N]; 
	double *y = new double[N];
	double norma;

	do{
		for (int i = 0; i < N; i++) 
			y[i] = x[i];
			
		//вычисляем невязку
		for(int i=0; i<N; i++){
			r[i]=-b[i];
			for(int j=0; j<N; j++)
				r[i]+=A[i][j]*x[j];
		}

		//найдем итерационный параметр
		double tau=0;
		double Ar_norm=0;
		for (int i = 0; i < N; i++) {
			double Ar=0;	
        	for (int j = 0; j < N; j++)
            	Ar += A[i][j] * r[j];

			tau+= Ar*r[i];
			Ar_norm+=Ar*Ar;
		}
		tau/=Ar_norm;

		for(int i=0; i<N; i++)
			x[i] = y[i] - tau * r[i];

		norma = 0;
    	for (int i = 0; i < N; i++) 
			if (abs(y[i] - x[i]) > norma)
				norma = abs(y[i] - x[i]);
	}	while(norma>eps);
}



/**
 * Метод сопряженных градиентов
 */
void loginovvv::lab7()
{
	const double eps = 1e-21;
	double *r = new double[N];
	double *z = new double[N];
	double *Az = new double[N];
	double alfa, betta, tmp, norma, b_norm=0;

	for(int i=0; i<N; i++){
		b_norm += b[i]*b[i]; 
		r[i] = b[i];
		for(int j=0; j<N; j++)
			r[i]-=A[i][j]*x[j];
		z[i]=r[i];
	}

	do{
		for(int i=0; i<N; i++){
			Az[i]=0;
			for(int j=0; j<N; j++)
				Az[i]+=A[i][j]*z[j];
		}
		alfa=0; tmp=0;
		for(int i=0; i<N; i++){
			alfa+=r[i]*r[i];
			tmp+=Az[i]*r[i];
		}
		alfa/=tmp;

		betta=0; tmp=0;
    	for(int i=0; i<N; i++){
			tmp+=r[i]*r[i];
			x[i]+=alfa*z[i];
			r[i]-=alfa*Az[i];
			betta+=r[i]*r[i];
		}
		norma=sqrt(betta/b_norm);
		betta/=tmp;
		
		for(int i=0; i<N; i++)
			z[i]=r[i]+betta*z[i];
	}while(norma>eps);
}


void loginovvv::lab8()
{
	double err = 0;
	double** C = new double*[N];
	double eps = 1e-20;

	for (int i = 0; i < N; i++){
		C[i] = new double[N];
	}


	for (int i = 0; i < N; i++){
		for (int j = i+1; j < N; j++){
    		if(i!=j){ 
      			err+=A[i][j]* A[i][j];
   			}
    		C[i][j] = 0;
  		}
	}

	while (sqrt(err) > eps){
		int mI = 0, mJ = 1;
		for (int i = 0; i < N; i++)
			for (int j = i+1; j<N; j++)
				if ( abs(A[i][j]) > abs(A[mI][mJ]) ){
  	    			mI = i; mJ = j;	
		}

	double phi;
  	if(A[mI][mI]!= A[mJ][mJ]){
    	phi = 0.5*atan(2*A[mI][mJ]/(A[mI][mI] - A[mJ][mJ]));
 	}else{phi = M_PI/4;}
	double c = cos(phi), s = sin(phi);
 
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
	    	if (k != mI && k != mJ && l != mI && l != mJ) 
				C[k][l] = A[k][l];	  
	}
	
	err = 0;
	for (int i = 0; i < N; i++)
      for (int j = i+1; j < N; j++)
  	    if (i != j) 
			err += C[i][j] * C[i][j]; 

	for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) 
	  	A[i][j] = C[i][j]; 

  }
  
  	for (int i = 0; i < N; i++) 
  		x[i] = A[i][i];

}



void loginovvv::lab9()
{
	double eps = 1.0e-3;
	double *y = new double[N];
    double *y_next = new double[N];
    double l = 0, l_next = 1;
    

    for (int i=0; i<N; i++)
  		y[i]=1;
  
    while(fabs(l_next - l) > eps){
  		l = l_next;
  	
		for(int i=0; i<N; i++){
	    	double s = 0;
        	for(int j=0; j<N; j++){
        		s += A[i][j]*y[j];
        	}
	    y_next[i]=s;
    	}
    
		for (int i=0; i<N; i++){
			if (y_next[0] != 0 && y[0] != 0){
				l_next = y_next[i]/y[i];
			}
		}
	
		for (int i=0; i<N; i++)
			y[i]=y_next[i];
   	
    }
    cout<<"RESULT: " << l_next << endl;
}


std::string loginovvv::get_name()
{
  return "Loginov V. V.";
}
