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
 * Метод простых итераций
 */
void edelevaup::lab5()
{

}



/**
 * Метод Якоби или Зейделя
 */
void
edelevaup::lab6 ()
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
void edelevaup::lab7()
{

}



/**
 * Метод сопряженных градиентов
 */
void edelevaup::lab8()
{

}


void edelevaup::lab9()
{

}


/*void edelevaup::lab10()
{

}
*/

std::string edelevaup::get_name()
{
  return "Edeleva U.P.";
}
