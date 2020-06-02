#include "golovatyukam.h"

/**
 * Введение в дисциплину
 */
void golovatyukam::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void golovatyukam::lab2()
{
	int m = N + 1;
	double eps =0.000001;
	
	double** arr = new double*[N];
	for (int i = 0; i < N; i++) {
		arr[i] = new double[N + 1];
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			arr[i][j] = A[i][j];
		}
	}

	for (int i = 0; i < N; i++) {
		arr[i][N] = b[i];
	}

	//cout << "\nПрямой ход: " << endl;
	for (int i = 0; i < N; i++) {

		//////////////////////////////////////////////////////////////
		//проаверяем наличие решения
		int max = i;
		for (int j = i + 1; j < N; j++) {
			if (abs(arr[j][i]) > abs(arr[max][i]))
				max = j;
		}
		if (max != i) {
			for (int k = 0; k < N + 1; k++) {
				swap(arr[i][k], arr[max][k]);
			}
		}

		if (abs(arr[i][i]) <= eps) {
			cout << "Данная система не имеет решений " << endl;
			exit(0);
		}
		/////////////////////////////////////////////////////////////

		double ved = arr[i][i];
		for (int j = i; j < N + 1; j++) {
			arr[i][j] /= ved;
		}
		for (int j = i + 1; j < N; j++) {
			ved = arr[j][i];
			arr[j][i] = 0;
			for (int k = i + 1; k < N + 1; k++)
				arr[j][k] -= arr[i][k] * ved;
		}
	}

	//cout << "\n\nОбратный ход: " << endl;

	for (int i = N - 1; i > 0; i--) {
		double s = 0;
		for (int j = i + 1; j < m; j++) {
			double ved = arr[i - 1][j - 1];
			arr[i - 1][j - 1] = 0;
			s += arr[j - 1][m - 1] * ved;
		}
		arr[i - 1][m - 1] = arr[i - 1][m - 1] - s;
	}
	
	for (int i = 0; i < N; i++) {
		x[i] = arr[i][N];
	}

}



/**
 * Метод прогонки
 */
void golovatyukam::lab3()
{
	int n =N;
	
	double *alfa = new double[n];
	double *betta = new double[n];
	double *F = new double[n];

	alfa[0] = -A[0][1] / A[0][0];
	betta[0] = b[0] / A[0][0];
	F[0] = -b[0];

	for (int i = 1; i < n - 1; i++) {
		F[i] = -1 * b[i];

		alfa[i] = -A[i][i + 1] / (A[i][i - 1] * alfa[i - 1] + A[i][i]);
		betta[i] = (-A[i][i - 1] * betta[i - 1] - F[i]) / (A[i][i - 1] * alfa[i - 1] + A[i][i]);
	}

	betta[n - 1] = (-A[n - 1][n - 2] * betta[n - 2] + b[n - 1]) / (A[n - 1][n - 2] * alfa[n - 2] + A[n - 1][n - 1]);

	x[n - 1] = betta[n - 1];

	for (int i = n - 2; i >= 0; i--) {
		x[i] = betta[i] + alfa[i] * x[i + 1];
	}
}



/**
 * Метод Холецкого
 */
void golovatyukam::lab4()
{
	int n = N;
	double *D = new double[n];

	if (A[0][0] > 0) D[0] = 1;
	else D[0] = -1;


	double **S = new double*[n];
	for (int i = 0; i < n; i++) {
		S[i] = new double[n];
	}
	S[0][0] = sqrt(fabs(A[0][0]));

	//s[0][j] - первая строка 
	for (int j = 1; j < n; j++) 
		S[0][j] = A[0][j] / (D[0] * S[0][0]);
	

	for (int i = 1; i < n; i++) {
		
		//--------------- find d[i][i] && s[i][i] ------------------------------
		double sum = 0;
		for (int k = 0; k < i; k++) {
			sum += D[k] * pow(S[k][i], 2);
		}
		D[i] = copysign(1, A[i][i] - sum);
		S[i][i] = sqrt(fabs(A[i][i] - sum));
		
		
		//-------------find s[i][j]---------------------------------- 
		for (int j = i + 1; j < n; j++) {
			double tmps = 0;
			for (int k = 0; k < j; k++) {
				tmps += D[k] * S[k][i] * S[k][j];
			}
			S[i][j] = (A[i][j] - tmps) / (D[i] * S[i][i]);
		}

	}

	double* y = new double[n];
	y[0] = b[0]/S[0][0];

	for (int i = 1; i < n; i++) {
    		double sum = 0;
    		for (int j = 0; j < i; j++) {
     			 sum += y[j]*S[j][i];
    		}
    		y[i] = (b[i] - sum) / S[i][i];
 	 }

  x[n-1] = y[n-1]/(S[n-1][n-1]*D[n-1]);

  for(int i = n-2; i>=0; i--){
    double sum = 0;
      for(int j = i+1; j<n; j++){
        sum+=x[j]*D[j]*S[i][j];
      }
      x[i] = (y[i] - sum)/(D[i]*S[i][i]);
  }
}



/**
 * Метод Якоби или Зейделя
 */
void golovatyukam::lab5()
{
int n = N;	
double eps = 1e-69, norma;
double *y = new double [n];
  do{
    for(int i = 0; i<n; i++)
      y[i]=x[i];
    
	norma = 0;
	
    for(int i = 0; i<n; i++){
      double sum1 = 0, sum2 = 0;
      for(int j = i + 1; j < n; j++)
      sum1 += A[i][j]*x[j];
    
      for(int j = i-1; j>= 0; j--)
      sum2 += A[i][j]*x[j];
    
      x[i] = (b[i] - sum1 - sum2)/A[i][i];
    }
  
    for (int i = 0; i < n; i++){
      norma += abs(x[i] - y[i]);
    }
  } while (norma>eps);
}



/**
 * Метод минимальных невязок
 */
void golovatyukam::lab6()
{
	int n = N;

	double* F = new double[n];
	double* r = new double[n];
	double* alfa = new double[n];
	double a, norma, k = 0;
	double eps = 1e-19;

	do {
		//r  - вектор невязки, F - функционал 
		for (int i = 0; i < n; i++) {
			double tmp = 0;
			for (int j = 0; j < n; j++) {
				tmp += A[i][j] * x[j];
			}
			r[i] = tmp - b[i];
			F[i] = 2 * r[i];
		}
		//alfa
		double* Ar = new double[n];
		for (int i = 0; i < n; i++) {
			double tmps = 0;
			for (int j = 0; j < n; j++) {
				tmps += A[i][j] * r[j];
			}
			Ar[i] = tmps;
		}
		double ts1 = 0, ts2 = 0;
		for (int i = 0; i < n; i++) {
			ts1 += abs(Ar[i] * r[i]);
			ts2 += abs(Ar[i] * Ar[i]);
		}
		a = ts1 / (2 * ts2);

		//x[k+1]
		double*y = new double[n];
		for (int i = 0; i < n; i++) {
			y[i] = x[i];
		}
		for (int i = 0; i < n; i++) {
			x[i] = x[i] - a * F[i];
		}

		norma = 0;
		for (int i = 0; i < n; i++) {
			norma += (y[i] - x[i])*(y[i] - x[i]);
		}
		//cout << k++ << " " << norma << endl;

	} while (sqrt(norma) > eps);
}



/**
 * Метод сопряженных градиентов
 */
void golovatyukam::lab7()
{

  int n = N;
  double* r = new double[n];
  double* r0 = new double[n];
  double  eps = 1e-21, norma_b, norma;
  double* Az = new double[n];
  double* z = new double[n];

  //for(int i = 0; i<n; i++)
  //x[i] = b[i];

  for(int i = 0; i<n; i++)
  norma_b += b[i]*b[i]; 

  //r0 = (Ах0 - b)
  for (int i = 0; i< n; i++){
    double tmp=0;
    for(int j = 0; j< n; j++){
      tmp +=A[i][j]*x[j]; 
    }
    r[i] = b[i] - tmp;
    z[i] = r[i];
  }

  double alfa, betta;
  //Az
  for (int i = 0; i< n; i++){
    double tmp=0;
    for(int j = 0; j< n; j++){
      tmp +=A[i][j]*z[j]; 
    }
    Az[i] = tmp;
  }
  //alfa
  double ts1 = 0, ts2 = 0;
    for (int i = 0; i<n; i++){
      ts1 += abs(r[i]*r[i]);
      ts2 += abs(Az[i]*z[i]); 
    }
  alfa = ts1/ts2;

  //запоминаем r[i]
  for(int i=0; i<n; i++){
    r0[i] = r[i];
  }

  for(int i=0; i<n; i++){
    x[i] = x[i] + alfa*z[i];
    r[i] = r[i] - alfa*Az[i];
  }

   for (int i = 0; i<n; i++){
      ts1 += abs(r[i]*r[i]);
      ts2 += abs(r0[i]*r0[i]); //abs(As[i]*r[i])
    }
  betta = ts1/ts2;

  for(int i = 0; i<n; i++){
    z[i] = r[i] + betta*z[i];
  }
  double tr = 0;
  for(int i = 0; i<n; i++){
    tr +=r[i]*r[i];
  }

  norma = sqrt(tr/norma_b);

do{
  //Az
  for (int i = 0; i< n; i++){
    double tmp=0;
    for(int j = 0; j< n; j++){
      tmp +=A[i][j]*z[j]; 
    }
    Az[i] = tmp;
  }
  //tao / alfa
  double ts1 = 0, ts2 = 0;
    for (int i = 0; i<n; i++){
      ts1 += abs(r[i]*r[i]);
      ts2 += abs(Az[i]*z[i]); //abs(As[i]*r[i])
    }
  alfa = ts1/ts2;

  //запоминаем r[i]
  for(int i=0; i<n; i++){
    r0[i] = r[i];
  }

  for(int i=0; i<n; i++){
    x[i] = x[i] + alfa*z[i];
    r[i] = r[i] - alfa*Az[i];
  }

   for (int i = 0; i<n; i++){
      ts1 += abs(r[i]*r[i]);
      ts2 += abs(r0[i]*r0[i]); //abs(As[i]*r[i])
    }
  betta = ts1/ts2;

  for(int i = 0; i<n; i++){
    z[i] = r[i] + betta*z[i];
  }
  double tr = 0;
  for(int i = 0; i<n; i++){
    tr +=r[i]*r[i];
  }

  norma = sqrt(tr/norma_b);

}while(norma>=eps);
 
}


void golovatyukam::lab8(){
int n = N;
double** C = new double*[n];
double eps = 1e-10,err = 0;
int cur_i, cur_j;

for (int i = 0; i < n; i++){
	C[i] = new double[n];
}

//вычисляем ошибку
for (int i = 0; i < n; i++){
  for (int j = i+1; j < n; j++){
    if(i!=j){ 
      err+=A[i][j];
    }
    C[i][j] = 0;
  }
}

while(err>eps){
  //cout<<1<<endl;
  double max = abs(A[0][1]);
  //находим max в верхнем треугольнике
  for (int i = 0; i < n; i++){
    for (int j = i+1; j < n; j++){
      if(abs(A[i][j]) > max){
        max = abs(A[i][j]);
        cur_i = i;
        cur_j = j;
      }
    }
  }
  double phi;
  if(A[cur_i][cur_i]!= A[cur_j][cur_j]){
    phi = 0.5*atan(2*A[cur_i][cur_j]/(A[cur_j][cur_j] - A[cur_i][cur_i]));
  }else{phi = M_PI/4;}
  double s = sin(phi), c = cos(phi);

  //используем ф-лы (10)
  C[cur_i][cur_j] = c*c*A[cur_i][cur_i] - 2*s*c*A[cur_i][cur_j] + s*s*A[cur_j][cur_j];
  C[cur_j][cur_j]	= s*s*A[cur_i][cur_i] + 2*s*c*A[cur_i][cur_j] + c*c*A[cur_j][cur_j];
	C[cur_i][cur_j] = (c*c - s*s)*A[cur_i][cur_j] + s*c*(A[cur_j][cur_j] - A[cur_i][cur_i]);
	C[cur_j][cur_i] = C[cur_i][cur_j];

  for (int k=0; k<n; k++){
	  if (k != cur_i  &&  k != cur_j){
	  	C[cur_i][k] = c*A[cur_i][k] - s*c*A[cur_j][k];
	  	C[k][cur_i] = C[cur_i][k];
	  	C[cur_j][k] = s*A[cur_i][k] + c*A[cur_j][k];//s*c??
	  	C[k][cur_j] = C[cur_j][k];
	  }

    for(int l = 0; l<n; l++)
      if(l!= cur_j && k!=cur_i)
         C[k][l] = A[k][l];
  }

  err = 0;
  for (int i = 0; i < n; i++)
    for (int j = i+1; j < n; j++)
      if(i!=j) err+=C[i][j];


  for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
	      A[i][j] = C[i][j];
}
    cout<<"СЗ Mat(A):"<<endl;
	for (int i=0; i<n; i++)
		cout<<A[i][i]<<endl;
    cout<<"))";

}


void golovatyukam::lab9(){
  int n = N;
  double *y = new double[n];
  double *y_next = new double[n];

  double eps =1e-2, lyambda = 1, lyambda_next = 0;

  for(int i = 0; i<n; i++)
    y[i] = b[i];//y0 

  do{
  for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
        y_next[i] += A[i][j]*y[j];
      }
  }
  lyambda = lyambda_next;

  for(int i=0; i<n; i++){
    if(y[i]!= 0 && y_next[i] != 0){
      lyambda_next = y_next[i]/y[i];
      break;
    }
  }
  for (int i=0; i<n; i++)
   	y[i]=y_next[i];
  }while(fabs(lyambda_next - lyambda)>eps);

cout<<"Answer: "<<lyambda_next << endl;
}



std::string golovatyukam::get_name()
{
  return "Golovatyukam A.M.";
}
