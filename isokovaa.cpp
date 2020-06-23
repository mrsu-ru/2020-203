#include "isokovaa.h"

/**
 * Введение в дисциплину
 */
void isokovaa::lab1()
{
  cout << "hello world!" << endl;
}

/**
 * Метод Гаусса с выбором главного элемента
 */
void isokovaa::lab2()
{
int *sw1=new int[N];
for (int i=0; i<N; i++){
sw1[i]=i;
}
for (int i=0; i<N; i++){
int Maximum=i;
for (int j=i+1; j <N; j++){
    if (abs(A[j][i])>abs(A[Maximum][i])) {
        Maximum=j;
}
}
if (Maximum!=i){
    swap(A[Maximum],A[i]);
    swap(b[Maximum],b[i]);
    swap(sw1[Maximum],sw1[i]);
}
b[i]/=A[i][i];
for (int j=N-1; j>=i; j--){
    A[i][j]/=A[i][i];
}
for (int j=0; j<N; j++){
    if (j!=i){
        b[j]-=A[j][i]*b[i];
        for (int k=N-1; k>=i; k--){
            A[j][k]-=A[j][i]*A[i][k];
}
}
}
}
for (int i=0; i<N; i++){
    x[i]=b[i];
}
for (int i=0; i<N; i++){
    if (sw1[i]!=i){
        swap(A[sw1[i]],A[i]);
        swap(b[sw1[i]],b[i]);
}
}
}

/**
 * Метод прогонки
 */
void isokovaa::lab3()
{
double *tg=new double[N];
double *vg=new double[N];
tg[0]=-A[0][1]/A[0][0];
vg[0]=b[0]/A[0][0];
for (int i=1; i<N; i++){
tg[i]=-A[i][i+1]/(A[i][i]+A[i][i-1]*tg[i-1]);
vg[i]=(b[i]-A[i][i-1]*vg[i-1])/(A[i][i]+A[i][i-1]*tg[i-1]);
}
x[N-1]=vg[N-1];
for (int i=N-2; i>=0; i--){
x[i]=tg[i]*x[i+1]+vg[i];
}
}

/**
 * Метод Халецкого
 */
void isokovaa::lab4()
{
double *d;
double *sum;
double v=0;
double *m;
int i,j,k;
d=new double[N];
sum=new double[N];
m=new double[N];
double **s=new double*[N];
for (i=0; i<N;i++){
s[i] = new double[N];
}
for (i = 0; i < N; i++){
sum[i] = A[i][i];
}
if (sum[0]>0){
    d[0]=1;
}
else {
    d[0]=-1;
}
s[0][0]=sqrt(fabs(sum[0]));
for (j=1;j<N;j++){
s[0][j]=A[0][j]/(d[0]*s[0][0]);
}
for (i=1; i<N; i++){
for (k=0; k<i; k++){
sum[i]-=d[k]*pow(s[k][i],2);
}
if (sum[i]>0){
    d[i]=1;
}
else{
d[i]=-1;
}
s[i][i] = sqrt(fabs(sum[i]));
for (j=i+1; j<N; j++){
for (k=0; k<i; k++){
v+=d[k]*s[k][i]*s[k][j];
}
s[i][j]=(A[i][j]-v)/(d[i]*s[i][i]);
v=0;
}
}
m[0]=b[0]/s[0][0];
for (i=1; i<N; i++){
for (k=0; k<i; k++){
v+=s[k][i]*m[k];
}
m[i]=(b[i]-v)/s[i][i];
v=0;
}
x[N-1] = m[N-1]/(d[N-1]*s[N-1][N-1]);
for (i=N-2; i>=0; i--){
for (k=i+1; k<N; k++){
v+=s[i][k]*x[k];
}
x[i]=(m[i]-d[i]*v)/(d[i]*s[i][i]);
v=0;
}
}

/**
 * Метод Якоби или Зейделя
 */
void isokovaa::lab5()
{
double epsil=1e-20;
for (int i=0; i<N; i++){
x[i]=0;
}
double *prev=new double[N];
double nr=0;
do{
for (int i=0; i<N; i++){
prev[i]=x[i];
}
for (int i=0; i<N; i++){
double result=b[i];
for (int j=0; j<N; j++){
if (i!=j){
result-=(A[i][j]*prev[j]);
}
}
x[i]=result/A[i][i];
}
nr=0;
for (int i=0; i<N; i++){
if (abs(prev[i]-x[i])>nr){
nr=abs(prev[i]-x[i]);
}
}
}
while (nr>epsil);
delete[] prev;
}

/**
 * Метод минимальных невязок
 */
void isokovaa::lab6()
{
double eps = 1e-17;
double *z, *Az;
int i,j;
z = new double[N];
Az = new double[N];
for (i=0; i<N; i++) {
        x[i]=0;
}
for (i=0; i<N; i++) {
        z[i]=b[i];
}
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
for (j=0; j<N; j++) Az[i] += A[i][j]*z[j];
}
scalar = 0;
for (i=0; i<N; i++) scalar += Az[i]*z[i];
t = -scalar;
scalar = 0;
for (i=0; i<N; i++) scalar += Az[i]*Az[i];
t/=scalar;
for (i=0; i<N; i++) x[i] = x[i] - t*z[i];
for (i=0; i<N; i++) z[i] = b[i];
for (i=0; i<N; i++)
for (j=0; j<N; j++)
z[i] -= A[i][j]*x[j];
Z=0;
for (i=0; i<N; i++) Z += z[i]*z[i];
}
delete[] Az;
delete[] z;
}



/**
 * Метод сопряженных градиентов
 */
void isokovaa::lab7()
{
// double *z;
// double *r;
// double *Az;
// double alpha, beta, eps = 1e-17;
// int i,j, sum = 1, k=0;
// z = new double[N];
// r = new double[N];
// Az = new double[N];
// for (i=0; i<N; i++){
// r[i] = b[i];
// x[i] = 0;
// z[i] = r[i];
// }
// sum = scal(r, r, N);
// do {
// for (i=0; i<N; i++){
// Az[i] = 0;
// for (j=0; j<N; j++) Az[i] += A[i][j] * z[j];
// }
// alpha = scal(r, r, N)/scal(Az, z, N);
// for (i=0; i<N; i++) {
//         x[i] = x[i] + alpha * z[i];
// }
// for (i=0; i<N; i++) r[i] = r[i] - alpha * Az[i];
// beta = sum/scal(r, r, N);
// for (i=0; i<N; i++) z[i] = r[i] + beta * z[i];
// sum = scal(r, r, N);
// k++;
// } while (scal(r,r,N) > eps*eps);
}


void isokovaa::lab8()
{
double err = 0;
double eps = 1e-6;
double max;
double alpha;
double c,s;
int i, j, I, J;
for (int i = 0; i < N; i++)
for (int j = i + 1; j < N; j++)
err += A[i][j] * A[i][j] + A[j][i] * A[j][i];
double **C;
C = new double*[N];
for (i=0; i<N; i++){
C[i] = new double[N];
}
while (err > eps){
alpha = 0;
max = abs(A[0][1]);
I = 0;
J = 1;
for (i = 0; i<N; i++)
for (j = i+1; j<N; j++)
if (abs(A[i][j]) > max){
max = abs(A[i][j]);
I = i; J = j;
}
if (abs(A[I][I] - A[J][J]) > eps) alpha = atan(2*A[I][J] / (A[J][J] - A[I][I])) / 2;
else alpha = atan(1);
c = cos(alpha);
s = sin(alpha);
for (i=0; i<N; i++)
for (j=0; j<N; j++)
if (i!=I && j!=I && i!=J && j!=J) C[i][j] = A[i][j];
C[I][I] = c*c*A[I][I] - 2*s*c*A[I][J] + s*s*A[J][J];
C[J][J] = s*s*A[I][I] + 2*s*c*A[I][J] + c*c*A[J][J];
C[I][J] = (c*c - s*s)*A[I][J] + s*c*(A[I][I] - A[J][J]);
C[J][I] = C[I][J];
for (int k=0; k<N; k++)
if (k!=I && k!=J){
C[I][k] = c*A[I][k] - s*A[J][k];
C[k][I] = C[I][k];
C[J][k] = s*A[I][k] + c*A[J][k];
C[k][J] = C[I][k];
}
for (i=0; i<N; i++)
for (j=0; j<N; j++)
A[i][j] = C[i][j];
err = 0;
for (i = 0; i < N; i++)
for (j = i + 1; j < N; j++)
err += C[i][j] * C[i][j] + C[j][i] * C[j][i];
}
for (i=0; i<N; i++) x[i] = A[i][i];
for (i = 0; i < N; i++) delete []C[i];
delete []C;
}


void isokovaa::lab9()
{
double eps = 1e-3;
double *yk = new double[N];
double *yk1 = new double[N];
int i,j;
for (int i = 0; i < N; i++) yk[i] = 1;
double err = 0;
double lambda1 = 0;
double lambda2 = 0;
do{
for (i = 0; i < N; i++)
for (j = 0; j < N; j++)
yk1[i] += A[i][j] * yk[j];
for (i = 0; i < N; i++)
if(fabs(yk[i]) > eps && fabs(yk1[i]) > eps){
lambda1 = yk1[i]/yk[i];
break;
}
err = fabs(lambda1 - lambda2);
lambda2 = lambda1;
for (i = 0; i < N; i++) {
yk[i] = yk1[i];
yk1[i] = 0;
}
}
while(err > eps);
cout<<"max lambda = "<<lambda2<<endl;
}


std::string isokovaa::get_name()
{
  return "Isokov A.A.";
}
