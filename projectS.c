#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
int main(int argc, char *argv[])
{
/* Initialize MPI */
int step,i,j,l,L;
double min,J,JJ;

/* random number generator initialization */
srandom(0);
int I = 100;
int K = 4;

/* data generation */
int *y = (int *) malloc(I*sizeof(int));  
int *x = (int *) malloc(I*sizeof(int));  
double *cx = (double *) malloc(2*K*sizeof(double));  
double *cy = (double *) malloc(2*K*sizeof(double));  
double *cxr = (double *) calloc(2*K,sizeof(double));  
double *cyr = (double *) calloc(2*K,sizeof(double));  
int *c = (int *) malloc(I*sizeof(int));  
time_t start,end; 
time(&start);
for (i = 0; i < I; i++)
x[i] = random()%64 + (i%4)*64;
for (i = 0; i < I; i++)
y[i] = random()%128 + (i%2)*128;
for (i = 0; i < K; i++)
{	l = random()%I ;
	cx[i] = x[l] ;
	cy[i] = y[l] ;
}
/* K Means */
step = 0 ;
while((step < 2) || (JJ-J>5))
{	step++;
	for(i=0 ; i<I ; i++)
		{
 		min = (x[i]-cx[0])*(x[i]-cx[0])+(y[i]-cy[0])*(y[i]-cy[0]);
		c[i] = 0;
		for (j=1 ; j<K ; j++)
			
		if((x[i]-cx[j])*(x[i]-cx[j])+(y[i]-cy[j])*(y[i]-cy[j]) <min)
		{
			c[i] = j;
			min = (x[i]-cx[j])*(x[i]-cx[j])+(y[i]-cy[j])*(y[i]-cy[j]);
		}
		}
		
	for (j=1 ; j<K ; j++)
		{	min = 0;
			cx[j]=0;
			cy[j]=0;
			for(i=0 ; i<I ; i++)
			{		
				if(c[i]==j)
				{
					cx[j]=cx[j]+x[i];
					cy[j]=cy[j]+y[i];
					min++ ;
				}
			}
			if (min != 0)
				{cx[j]=cx[j]/min;
				 cy[j]=cy[j]/min;
				}
			cx[j+K]=min;
			cy[j+K]=min;
		}

	JJ=J;
	J=0;

	for (i=0;i<I;i++)
		{
		J=J+(x[i]-cx[c[i]])*(x[i]-cx[c[i]])+(y[i]-cy[c[i]])*(y[i]-cy[c[i]]);
		}
}

time(&end);
printf("%f", (double)(end-start));
/* Printing the result */
for (j=0; j<K;j++)
	{
	printf("%f,%f \n",cx[j],cy[j]);
	}

/* Wrinting the final solution */
/*FILE *fp;
int s = 0 ;
fp = fopen("project.txt","w");
for(j=0;j<K;j++)
	fprintf(fp, "%f,%f,",cx[j],cy[j]);
for (i=0;i<I;i++)
	{fprintf(fp, "%d,%d,%d,", x[i],y[i],c[i]);}
fclose(fp); */


free(x); 
free(y);
free(cx); 
free(cy);
free(cxr); 
free(cyr);
free(c);
exit(0);}