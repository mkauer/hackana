//
//
//  Library to calculate 90% limit, using O. Helene formula
//  Author Vladimir Vasiliev, 2002
//  Usage
//  SolveHelen(CL,mean_background, observation)
//
//

#include <stdlib.h>
//#include <conio.h>
#include <math.h>
//#include <mem.h>
//#include <process.h>
//#include <alloc>
#include <ctype.h>
#include <string>
#include "helen.hpp"


#define maxnumber  1000
#define e_step 0.001
#define eps 0.00001
#define f_eps 0.001
#define MAXITERATIONS 1000
#define infinity 1000000

#define STDDEV 1.0
#define CL90   1.9
#define HELEN90 2.30
#define EXPERIMENT  1
#define THEORETICAL 2
#define FIT         3
#define BOUND 		4
#define INIFILE     5

const long double C1 = 0.5*log(2.*(long double)3.14159265358979);
const long double B1 = 1./(12.);
const long double B2 = -1./(12*30);
const long double B3 = 1./(30*42);
const long double B4 = -1./(42*30);
const long double B5 = 5./(56*66);


FILE   *in_exper,*in_f1,*in_f2,*in_theory,*in_str;
long double a1,a2,dpa1,dpa2,dma1,dma2,cla1,cla2,F;
long double *exper,*f1,*f2,*theory;
long double scale1=1,scale2=1,thscale = 1;
int N,nmin = 0,nmax = 0;
long double himin = 0;
long double cl90;
char **pp;

long double sgn(long double f){
	if(f<0) return -1;
	if(f>0) return 1;
	return 1;
}
long double value_cl90(){
	long double K=0;
	long double a=1.6889;
	long double b=0.260035;
	long double c=0.0047;
	long double cl;
	int i;
	for(i=nmin;i<nmax;i++){
		K +=exper[i];
	}
	cl = a+b*pow(M_E,-c*K);
	return(cl);
}
long double R(long double x1,long double x2)
{
	return fabs(x1-x2);
}
long double pow(long double x,long double y)
{
	if(x<=0) return 0;
	long double dummyx=expl(1000);
	return expl(y*logl(x));
}
long double l_factorial(long double n)
{
	if(n<= 0) return(1.);
	else
	return (C1 - n + (n+0.5)*logl(n)+(B1+(B2+(B3+(B4+B5/n)/n)/n)/n)/n );
}
long double i_factorial(int n)
{
	long double i,f=1;
	for(i=1;i<=n;i++){
		f*=i;
	}
	return f;
}
long double factorial(int x)
{
	if(x<10) return logl(i_factorial(x));
	else return l_factorial(x);
}
long double FHelen(long double N,long double mb,int n0)
{
	long double dx,dy;
	int i;
	long double f1=0,f2=0;
	for (i = 0;i<=n0;i++){
		if(mb+N!=0)
			f1+= expl(i*logl(mb+N) - factorial(i));
		if(mb!=0)
			f2+= expl(i*logl(mb) - factorial(i));
	}
	if(mb==0) f1=expl(-N)*f1;
	else 	  f1 = expl(-N)*f1/f2;
	return f1;
}
long double SolveHelen(long double e,long double mb,int n0)
{
	int i;
	long double x0,x1,x2,f1,f2,f0;
	e=1-e;
	x1=1.;
	x2=5*n0+10.;
	x0 = (x1+x2)/2;
//	if(n0==0) return HELEN90;
	f1=FHelen(x1,mb,n0);
	f1-=e;
	do{
		x2*=2;
		f2=FHelen(x2,mb,n0)-e;
	}while(f1*f2>0);
	f0=FHelen(x0,mb,n0)-e;
	do{
		if(f1>0 && f0 <0){
			x2 = x0;
			f2 = f0;
		}
		if(f0>=0 && f2<0){
			x1=x0;
			f1=f0;
		}
		x0=(x1+x2)/2;
		f0=FHelen(x0,mb,n0)-e;
	}while(fabsl(x1-x2)>eps);
	return (x1+x2)/2;
}
/*main(int np ,char** par){
	float mb,e;
	float res;
	int n0;
	printf("\r\n\r\nHelen calculator. v1.1 by Vasilyev V.\r\n");
	printf("Calculates limits by O.Helen's formula. \r\n\r\n");
	printf("Input background (mb): ");
	scanf( "%f",&mb);
	printf("\r\nInput n. of events (n0) : ");
	scanf("%i",&n0);
	printf("\r\nInput CL value (e): ");
	scanf("%f",&e);
	res=SolveHelen(e,mb,n0);
	printf("\r\n\r\n Signal s < %f with CL=%f",res,e);

}
*/
