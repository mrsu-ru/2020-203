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
 * Метод Холецкого
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

}



/**
 * Метод сопряженных градиентов
 */
void isokovaa::lab7()
{

}


void isokovaa::lab8()
{

}


void isokovaa::lab9()
{

}


std::string isokovaa::get_name()
{
  return "Isokov A.A.";
}
