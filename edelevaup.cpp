//#include "zhalninrv.h"
#include "edelevaup.h"

/**
 * Введение в дисциплину
 */
void edelevaup::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void edelevaup::lab2()
{

	for (int i = 0; i < N; i++) {

		int max = i;
		for (int j = i + 1; j < N; j++) {
			if (abs(A[j][i]) > abs(A[max][i]))
				max = j;
		}
		double *s;
		double k;
		if (max != i) {
			s = A[i];
			A[i] = A[max];
			A[max] = s;
			k = b[i];
			b[i] = b[max];
			b[max] = k;
		}
		
		b[i] /= A[i][i];
		for (int j = N - 1; j > i; j--)
			A[i][j] /= A[i][i];
		A[i][i] = 1;

		for (int j = i + 1; j < N; j++) {

			for (int k = N; k > i; k--) {
				A[j][k] -= A[i][k] * A[j][i];
			}
			b[j] -= b[i] * A[j][i];

			A[j][i] = 0;
		}


	}

	for (int i = N - 1; i > 0; i--) {
		double s = 0, l = 0;
		for (int j = i; j < N; j++) {
			s += b[j ] * A[i - 1][j];
			A[i - 1][j] = 0;
			
		}
		b[i - 1] -= s;
		
	}

	for (int i = 0; i < N; i++) {
		x[i] = b[i];
		
		}
	

}



/**
 * Метод прогонки
 */
void edelevaup::lab3()
{
	double *al = new double[N];
    double *bet = new double[N];
    al[0]=A[0][1]/-A[0][0]; 
	bet[0]=b[0]/A[0][0];
	for( int i=1; i<N; i++){
		al[i]=A[i][i+1]/(-A[i][i]-A[i][i-1]*al[i-1]);
		bet[i]=(A[i][i-1]*bet[i-1]-b[i])/(-A[i][i]-A[i][i-1]*al[i-1]);
}
 
  x[N-1]=bet[N-1];
  for (int i=N-2; i>=0; i--){
	  x[i]=al[i]*x[i+1]+bet[i];
  }
 
}
/**
 * Метод Холецкого (метод квадратного корня)
 */
void edelevaup::lab4()
{
  double *D = new double[N];
  double **S = new double* [N];  
   double *y = new double[N];
   //double l=0, temp=0;
  for (int i=0; i<N; i++){
	  S[i]=new double[N];
	  for(int j=0; j<N; j++)
		S[i][j]=0;
  } 
   if (A[0][0]>0)
		  D[0]=1;
	  else 
		  D[0]=-1;
	  S[0][0]=sqrt(fabs(A[0][0]));
	  
   for (int i=1; i<N; i++){
   S[0][i]=A[0][i]/(D[0]*S[0][0]);
   }
   
	  for (int i=1; i<N; i++){
	  //temp = A[i][i];
	  double temp =0;
      for (int j=0; j<i; j++){
		  temp+=D[j]*S[j][i]*S[j][i];
	  }
	  if (A[i][i]-temp>=0)
		  D[i]=1;
	  else 
		  D[i]=-1;
	  S[i][i]=sqrt(D[i]*(A[i][i]-temp));
	  
	  for (int j=i+1; j<N; j++) {
		//  l =A[i][j];
			  double l =0;

	   for (int k=0; k<j; k++){
	   l+=D[k]*S[k][i]*S[k][j];
	   }
	   
		  S[i][j]=(A[i][j]-l)/(D[i]*S[i][i]);
	   }
		  }
		  
		  y[0]=b[0]/S[0][0];
		  for (int i=1; i<N; i++){
			  	  double temp =0;

			  for (int j=0; j<i; j++){
				  temp+=y[j]*S[j][i];
			  }
			  y[i]=(b[i]-temp)/S[i][i];
		  }		
	 x[N-1]=y[N-1]/(D[N-1]*S[N-1][N-1]);
	 
for (int i=N-2; i>=0; i--){
		  double temp =0;
	for (int j=i+1; j<N; j++){
		temp+=x[j]*D[j]*S[i][j];
	}
	x[i]=(y[i]-temp)/(D[i]*S[i][i]);
}
			  
}





/**
 * Метод Якоби или Зейделя
 */
void edelevaup::lab5 ()
{
  double *f = new double[N];
  double e = pow(10, -30);
  double temp =0;
   
   do{
		for (int i = 0; i < N; i++)
		{
			f[i] = x[i];
		}
     
      for (int i = 0; i < N; i++) {

	  double sum1 = 0, sum2 = 0;

	  for (int j = 0; j < i; j++)

	    sum1 += A[i][j] * x[j];

	  for (int j = i + 1; j < N; j++)

	    sum2 += A[i][j] * x[j];

	  x[i] = (b[i] - sum1 - sum2) / A[i][i];
	}
	
	temp =0;
	  for (int i = 0; i < N; i++)
    {

      if (abs(x[i] - f[i]) > temp)
	      temp = abs (x[i] - f[i]);
    }
} while (temp>e);
	
}


/**
 * Метод минимальных невязок
 */
void edelevaup::lab6()
{
double eps = 1e-19;
	double *x1 = new double[N];	
	double *r = new double[N];	
	double *Ar = new double[N];	
	double norma =0;

	

	for (int i = 0; i < N; i++)
	{
		x[i] = b[i];
		x1[i]=0;

	}
	do{
	for (int i=0; i<N; i++){ 
		double t=0;
		for (int j=0; j<N; j++){
			t+=A[i][j]*x[j];
		}
		r[i]=t-b[i];
	}
	for (int i=0; i<N; i++){ 
		double t=0;
		for (int j=0; j<N; j++){
			t+=A[i][j]*r[j];
		}
		Ar[i]=t;
	}
	double s1=0, s2 =0;
	for (int i=0;i<N;i++){
  		s1 += r[i]*Ar[i];
		s2 += Ar[i]*Ar[i];
  	}
	double tau=s1/s2;
	for(int i=0;i<N;i++){
  		x[i]=x[i]-tau*r[i];
		  
	  }
	norma =0;
	for (int i=0;i<N;i++){
		norma+=(x1[i]-x[i])*(x1[i]-x[i]);
  		x1[i]=x[i];
  	}
	} while (sqrt(norma) > eps);
	
}


/**
 * Метод сопряженных градиентов
 */
void edelevaup::lab7()
{
  int i,j; double e = 1.0e-21;
 double *c, *r, *z, *r1;
 
  z = new double[N];
  c = new double[N];
  r = new double[N];
  r1 = new double[N];
  
  for (int i=0;i<N;i++){
  	x[i]=b[i];
  }
  
	for (i=0;i<N;i++){
  	  double s = 0;
  	  for (j=0;j<N;j++){
  		 s += A[i][j]*x[j];
	   }
	  r[i]=-s+b[i];
	  z[i]=-s+b[i];
      }
  
  double xxx = 1;
  double k=0;
  while (xxx>e){
	double al1=0,al2=0;
	k++;

  	for (i=0;i<N;i++){
  		al1 += r[i]*r[i];
  	}
  	
  	for (i=0;i<N;i++){
  	  double s = 0;
  	  for (j=0;j<N;j++){
  		 s += A[i][j]*z[j];
	   }
	  c[i]=s;
      }

    for (i=0;i<N;i++){
  		al2 += c[i]*z[i];
  	}
   double al=al1/al2;
   cout<<al<<endl;
  
   for (i=0;i<N;i++){
  	  x[i]+=al*z[i];
    }
   
   for (i=0;i<N;i++){
  	  r1[i] = r[i] - al*c[i];
    }
    double bt1=0;
  	for (i=0;i<N;i++){
  		bt1 += r1[i]*r1[i];
  	}
  	double bt=bt1/al1;
  	double u1=0,u2=0;
  	for (i=0;i<N;i++){
  	  z[i] = r1[i] + bt*z[i];
  	  r[i] = r1[i];
  	  u1+=r1[i]*r1[i];
  	  u2+=b[i]*b[i];
    }
    xxx = sqrt(u1)/sqrt(u2);
    cout<<xxx<<endl;
  
  }
}
// Метод вращения для нахождения собственных значений матрицы 

void edelevaup::lab8()
{
    double eps = 1.0e-6;
    double errc = 0;
    double c = 0, s = 0, alpha = 0;
    int max_i = 0, max_j = 1;
    double** C;

    C = new double*[N];
 
		for (int i = 0; i < N; i++)
	{
	     C[i] = new double[N];
        for (int j = i+1; j < N; j++)
        {
          
            if (i != j){
                errc += A[i][j]*A[i][j];
            }
            	C[i][j] = 0;
        }
    }
    while (sqrt(errc) > eps)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = i + 1; j < N; j++)
            {
                if (abs(A[i][j]) >= abs(A[max_i][max_j]))
                {

                    max_i = i;
                    max_j = j;
                }
            }
        }
        if (A[max_i][max_i] == A[max_j][max_j])
            alpha = M_PI / 4;
		else alpha = 0.5*atan((2 * A[max_i][max_j]) / (A[max_j][max_j] - A[max_i][max_i]));
        c = cos(alpha);
        s = sin(alpha);

        C[max_i][max_i]= c*c*A[max_i][max_i] - 2*s*c*A[max_i][max_j]+s*s*A[max_j][max_j];
        C[max_j][max_j]= s*s*A[max_i][max_i]+2*s*c*A[max_i][max_j]+c*c*A[max_j][max_j];
        C[max_i][max_j]= (c*c-s*s)*A[max_i][max_j]+s*c*(A[max_i][max_i]-A[max_j][max_j]);
        C[max_j][max_i]= C[max_i][max_j];

        for (int k = 0; k < N; k++)
        {
            if (k != max_i && k != max_j)
            {
                C[max_i][k] = c * A[max_i][k] - s * c * A[max_j][k];
                C[k][max_i] = C[max_i][k];
                C[max_j][k] = s * A[max_i][k] + c * A[max_j][k];
                C[k][max_j] = C[max_j][k];
            }
            for (int l = 0; l < N; l++)
            {
                if (k != max_i && k != max_j && l != max_i && l != max_j)
                    C[k][l] = A[k][l];
            }
        }
    
    errc = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            if (i != j)
            {
                errc += C[i][j] * C[i][j];
            }
        }
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i][j] = C[i][j];
        }
    }

    for (int i = 0; i < N; i++)
        x[i] = A[i][i];
}
}

//Наибольшее по модулю собственное значение матрицы
void edelevaup::lab9()
{
  double eps = 1.0e-3;
  double *yk = new double[N];	
  double *yk_1 = new double[N];
  double y0 = 0;
  double y1 = 0;
  double lambda_0 = 0;
  double lambda_1 = 0;
  double delta= 0;
  for ( int i = 0; i < N; i++)
	{
		yk_1[i] = 1;
		yk[i] = 0;
		
	}
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			yk[i] += A[i][j] * yk_1[j];
		}
	}
	y0 = yk_1[0];

	for (int i = 0; i < N; i++)
	{
		if (yk[i] != 0)
		{
			y1 = yk[i];
			break;
		}
	}

	lambda_0 = y1 / y0;
	delta = lambda_0;
	while (delta > eps)
	{
		for (int i = 0; i < N; i++)
		{
			yk_1[i] = yk[i];
			yk[i] = 0;
		}
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				yk[i] += A[i][j] * yk_1[j];
			}
		}
		for (int i = 0; i < N; i++)
		{
			if ((yk[i] != 0) && (yk_1[i] != 0))
			{
				y0 = yk_1[i];
				y1 = yk[i];
				break;
			}
		}
		lambda_1 = y1 / y0;
		delta = fabs(lambda_1 - lambda_0);
		lambda_0 = lambda_1;
	}

	cout << "Наибольшее по модулю СЗ матрицы A: " << lambda_1<<endl;
}


std::string edelevaup::get_name()
{
  return "Edeleva U.P.";
}
