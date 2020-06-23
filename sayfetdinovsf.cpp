#include "sayfetdinovsf.h"

/**
 * �������� � ����������
 */
void sayfetdinovsf::lab1()
{
  cout << "hello world!" << endl;
}


/**
 * ����� ������ � ������� �������� ��������
 */
void sayfetdinovsf::lab2()
{
int i,j,k,max;

  for (i=0; i<N; b[i]=-b[i], i++)

  for (i=0; i<N; i++) {
      max=i;
      double* s; double s1;
      for (j=i+1; j<N; j++) if (A[j][i]*A[j][i]>A[max][i]*A[max][i]) max=j;
      if (max!=i) {
         s=A[i]; s1=b[i];
         A[i]=A[max]; b[i]=b[max];
         A[max]=s; b[max]=s1;
      }
      
      for (j=N-1; j>i; A[i][j]/=A[i][i], j--);
      b[i]/=A[i][i];
      A[i][i]=1;
      
      for (j=N-1; j>i; j--) {
          for (k=i+1; k<N; A[j][k]-=A[i][k]*A[j][i], k++);
          b[j]-=b[i]*A[j][i];
          A[j][i]=0;
      }
      
      
  }
  
  double sum;
  x[N-1]=b[N-1];
  for (i=N-2;i>=0;i--){
    sum=0;
    for (j=i+1;j<N;j++) sum+=x[j]*A[i][j];
    x[i]=b[i]-sum;
  }
}



/**
 * ����� ��������
 */
void sayfetdinovsf::lab3()
{
	double *p = new double[N];
	double *q = new double[N];

	// ������ ���
	p[0] = -A[0][1] / A[0][0];
	q[0] = b[0] / A[0][0];
	for (int i = 1; i < N; i++) {
		p[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * p[i - 1]);
		q[i] = (b[i] - A[i][i - 1] * q[i - 1]) / (A[i][i] + A[i][i - 1] * p[i - 1]);
	}

	// �������� ���
	x[N - 1] = q[N - 1];
	for (int i = N - 2; i >= 0; i--) {
		x[i] = p[i] * x[i + 1] + q[i];
	}

	delete[] p;
	delete[] q;
}



/**
 * Метод Холецкого
 */
void sayfetdinovsf::lab4()
{
double eps = 1e-15;
double t = 1e-5;
for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double x1;
	double *xr = new double[N];
	int step = 0;

	do {
		step++;
		for (int i = 0; i < N; i++) {
			xr[i] = x[i];
			for (int k = 0; k < N; k++)
				xr[i] -= t*A[i][k] * x[k];
			xr[i] += t * b[i];

		}
		x1 = 0.;
		for (int i = 0; i < N; i++) {
			x1 += (xr[i]-x[i])*(xr[i]-x[i]);
		}

		for (int i = 0; i < N; i++) {
			x[i] = xr[i];
		}
	} while (sqrt(x1)>eps);
}



/**
 * ����� ����� ��� �������
 */
void sayfetdinovsf::lab5()
{
//����� �����
	double *oldx = new double[N]; 
	for (int i=0; i<N; i++) { 
		x[i]=0;} // ��������� ������� ������ 
	double Err=0.0; 
	double eps=1e-20; // �����������
	int k=0; 
	do { 
		k++; 
		Err=0.0; 
		for(int i=0; i<N; i++) 
			oldx[i]=x[i]; // ���������� ������� 
			for(int i=0; i<N; i++) 
			{ 
				double s=0;
				for(int j=0; j<i; j++) 
					s += A[i][j] * oldx[j]; //��� ������� ����������
				for(int j=i+1; j<N; j++) 
					s += A[i][j] * oldx[j]; //��� ������� ����������
				x[i]=(b[i] - s)/A[i][i]; // ����������� ����� ������� 
			}			 
			Err=std::abs(oldx[0]-x[0]); 
			for(int i=0; i<N; i++) 
			{ 
				if(std::abs(oldx[i]-x[i]) > Err) 
				Err = std::abs(oldx[i]-x[i]);//������������ ������� ����� ���������� �������� � �������. 
			} 
	} while(Err >= eps); 
delete [] oldx; 
}



/**
 * ����� ����������� �������
 */
void sayfetdinovsf::lab6()
{
double eps = 1e-17;    
    double *z, *Az;
    int i,j;
    z = new double[N];
    Az = new double[N];
    for (i=0; i<N; i++) x[i]=0;
    for (i=0; i<N; i++) z[i]=b[i];
	
    for (i=0; i<N; i++)
	   for (j=0; j<N; j++) 
	      z[i]-=A[i][j]*x[j];

    double Z=0;
    double t=0;
    double scalar=0;
	
    for (i=0; i<N; i++)
	    Z += z[i]*z[i];
    
    while (Z > eps*eps){
	    for (i=0; i<N; i++){
	        Az[i] = 0;
	        for (j=0; j<N; j++)	Az[i] += A[i][j]*z[j];
	    }
		
	    scalar = 0;
	    for (i=0; i<N; i++) scalar += Az[i]*z[i];
	    t = -scalar;
	    scalar = 0;
	    for (i=0; i<N; i++) scalar += Az[i]*Az[i];
	    t/=scalar;
	
	    for (i=0; i<N; i++) x[i] = x[i] - t*z[i];
			
	    for (i=0; i<N; i++)  z[i] = b[i];
            for (i=0; i<N; i++)
	            for (j=0; j<N; j++) 
		            z[i] -= A[i][j]*x[j];		
		
	    Z=0;
	    for (i=0; i<N; i++)  Z += z[i]*z[i];
	}
	
	delete[] Az;
	delete[] z;
}



/**
 * ����� ����������� ����������
 */
void sayfetdinovsf::lab7()
{
	double eps = 1e-20;
	double* prevX = new double[N];
	double* prevR = new double[N];
	double* r = new double[N];
	double* z = new double[N];
	for (int i = 0; i < N; i++) {
		r[i] = b[i];
		z[i] = r[i];
	}

	while (true) {

		for (int i = 0; i < N; i++) {
			prevR[i] = r[i];
			prevX[i] = x[i];
		}

		double alpha = 0, denAlpha = 0;

		for (int i = 0; i < N; i++) {
			double Az = 0;
			for (int j = 0; j < N; j++) {
				Az += A[i][j] * z[j];
			}
			alpha += prevR[i] * prevR[i];
			denAlpha += Az * z[i];
		}
		alpha /= denAlpha;

		for (int i = 0; i < N; i++) {
			x[i] = prevX[i] + alpha * z[i];
		}

		double maxErr = abs(x[0] - prevX[0]);
		for (int i = 1; i < N; i++)
			if (abs(x[i] - prevX[i]) > maxErr)
				maxErr = abs(x[i] - prevX[i]);

		if (maxErr < eps)
			break;

		for (int i = 0; i < N; i++) {
			double Az = 0;

			for (int j = 0; j < N; j++) {
				Az += A[i][j] * z[j];
			}

			r[i] = prevR[i] - alpha * Az;
		}

		double beta = 0, denBeta = 0;
		for (int i = 0; i < N; i++) {
			beta += r[i] * r[i];
			denBeta += prevR[i] * prevR[i];
		}
		beta /= denBeta;

		for (int i = 0; i < N; i++) {
			z[i] = r[i] + beta * z[i];
		}
	}

	delete[] prevX;
	delete[] r;
	delete[] prevR;
	delete[] z;
}


void sayfetdinovsf::lab8()
{

}


void sayfetdinovsf::lab9()
{

}


std::string sayfetdinovsf::get_name()
{
  return "sayfetdinovsf";
}
