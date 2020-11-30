#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
int main(int argc, char *argv[])
{
/* Initialize MPI */
int p,P,step,i,j;
double min,J,JJ,Jg;
MPI_Status status ;
MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &P);
MPI_Comm_rank(MPI_COMM_WORLD, &p);

/* random number generator initialization */
srandom(p);
int N = 1000;
int I = (int)((N+P-p-1)/P); 
int K = 2;
int D = (int)(log(P)/log(2));
int DD = pow(2,D);
/* data generation */
int *y = (int *) malloc(I*sizeof(int));  
int *x = (int *) malloc(I*sizeof(int));  
double *cx = (double *) malloc(2*K*sizeof(double));  
double *cy = (double *) malloc(2*K*sizeof(double));  
double *cxr = (double *) calloc(2*K,sizeof(double));  
double *cyr = (double *) calloc(2*K,sizeof(double));  
int *c = (int *) malloc(I*sizeof(int));  

for (i = 0; i < I; i++)
x[i] = random()%64 + (i%4)*64;
for (i = 0; i < I; i++)
y[i] = random()%128 + (i%2)*128;
for (i = 0; i < K; i++)
cx[i] = (random()%256) ;
for (i = 0; i < K; i++)
cy[i] = (random()%256) ;
/* K Means */
step = 0 ;
while((step < 2) || (JJ-J>20))
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

	if(p>=DD)
		{
		 MPI_Send(&cx[0],2*K,MPI_DOUBLE,p^(int)(DD),0,MPI_COMM_WORLD) ;
		 MPI_Send(&cy[0],2*K,MPI_DOUBLE,p^(int)(DD),0,MPI_COMM_WORLD) ;
		}
	if(p<P-DD)
		{
		MPI_Recv(&cxr[0],2*K,MPI_DOUBLE,(int)(p^(int)(DD)),0,MPI_COMM_WORLD,&status) ;
		MPI_Recv(&cyr[0],2*K,MPI_DOUBLE,p^(int)(DD),0,MPI_COMM_WORLD,&status) ;
		for (j=0; j<K;j++)
			{
				if ((int)(cx[j+K]) != (int)(cxr[j+K]))
				{
				cx[j]= cx[j+K]*(cx[j]/(cx[j+K]+cxr[j+K]))+cxr[j+K]*(cxr[j]/(cx[j+K]+cxr[j+K]));
				cy[j]= cy[j+K]*(cy[j]/(cy[j+K]+cyr[j+K]))+cyr[j+K]*(cyr[j]/(cy[j+K]+cyr[j+K]));
				}
				else 
				{
				cx[j] = 0.5*(cx[j]+cxr[j]);
				cy[j] = 0.5*(cy[j]+cyr[j]);
				}		
			}
		}
	if(p<DD){
	for(i=0;i<D;i++)
		{if (p>(int)(p^(int)(pow(2,i))))
		{
			MPI_Recv(&cxr[0],2*K,MPI_DOUBLE,(int)(p^(int)(pow(2,i))),0,MPI_COMM_WORLD,&status) ;
			MPI_Send(&cx[0],2*K,MPI_DOUBLE,p^(int)(pow(2,i)),0,MPI_COMM_WORLD) ;
			MPI_Recv(&cyr[0],2*K,MPI_DOUBLE,p^(int)(pow(2,i)),0,MPI_COMM_WORLD,&status) ;
			MPI_Send(&cy[0],2*K,MPI_DOUBLE,p^(int)(pow(2,i)),0,MPI_COMM_WORLD) ;
			for (j=0; j<K;j++)
			{
				if ((int)(cx[j+K]) != (int)(cxr[j+K]))
				{
				cx[j]= cx[j+K]*(cx[j]/(cx[j+K]+cxr[j+K]))+cxr[j+K]*(cxr[j]/(cx[j+K]+cxr[j+K]));
				cy[j]= cy[j+K]*(cy[j]/(cy[j+K]+cyr[j+K]))+cyr[j+K]*(cyr[j]/(cy[j+K]+cyr[j+K]));
				}
				else 
				{
				cx[j] = 0.5*(cx[j]+cxr[j]);
				cy[j] = 0.5*(cy[j]+cyr[j]);
				}		
			}
		}
		else 
		{
			MPI_Send(&cx[0],2*K,MPI_DOUBLE,(int)(p^(int)(pow(2,i))),0,MPI_COMM_WORLD) ;
			MPI_Recv(&cxr[0],2*K,MPI_DOUBLE,p^(int)(pow(2,i)),0,MPI_COMM_WORLD,&status) ;
			MPI_Send(&cy[0],2*K,MPI_DOUBLE,p^(int)(pow(2,i)),0,MPI_COMM_WORLD) ;
			MPI_Recv(&cyr[0],2*K,MPI_DOUBLE,p^(int)(pow(2,i)),0,MPI_COMM_WORLD,&status) ;
			for (j=0; j<K;j++)
			{   if ((int)(cx[j+K]) != (int)(cxr[j+K]))
				{
				cx[j]= cx[j+K]*(cx[j]/(cx[j+K]+cxr[j+K]))+cxr[j+K]*(cxr[j]/(cx[j+K]+cxr[j+K]));
				cy[j]= cy[j+K]*(cy[j]/(cy[j+K]+cyr[j+K]))+cyr[j+K]*(cyr[j]/(cy[j+K]+cyr[j+K]));
				}
				else 
				{
				cx[j] = 0.5*(cx[j]+cxr[j]);
				cy[j] = 0.5*(cy[j]+cyr[j]);
				}
			}
			}
		}
	}
	if(p<P-DD) 
		{MPI_Send(&cx[0],2*K,MPI_DOUBLE,p^(int)(DD),0,MPI_COMM_WORLD) ;
		 MPI_Send(&cy[0],2*K,MPI_DOUBLE,p^(int)(DD),0,MPI_COMM_WORLD) ;
		}
	if(p>=DD)
		{
		MPI_Recv(&cxr[0],2*K,MPI_DOUBLE,(int)(p^(int)(DD)),0,MPI_COMM_WORLD,&status) ;
		MPI_Recv(&cyr[0],2*K,MPI_DOUBLE,p^(int)(DD),0,MPI_COMM_WORLD,&status) ;
		for (j=0; j<K;j++)
			{	cx[j]= cxr[j];
				cy[j]= cyr[j];
			}
		}

	JJ=J;
	J=0;
 	for (i=0;i<I;i++)
		{
		J=J+(x[i]-cx[c[i]])*(x[i]-cx[c[i]])+(y[i]-cy[c[i]])*(y[i]-cy[c[i]]);
		}

	MPI_Allreduce(&J,&Jg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 
	J = Jg ;
}

/* Printing the result */
	for (j=0; j<K;j++)
	{
	printf("%d : %f,%f \n",p,cx[j],cy[j]);
	}

/* Wrinting the final solution */
/*if (p==0) 
	{FILE *fp;
	int s = 0 ;
	fp = fopen("project.txt","w");
	for(j=0;j<K;j++)
		fprintf(fp, "%f,%f,",cx[j],cy[j]);
	for (i=0;i<I;i++)
	{fprintf(fp, "%d,%d,%d,", x[i],y[i],c[i]);}
	fclose(fp);
	MPI_Send (&s,1,MPI_INT,p+1,0,MPI_COMM_WORLD);} //Sending a signal to the next process
else 
	{FILE *fp;
	int s = 0 ;
	MPI_Recv (&s,1,MPI_INT,p-1,0,MPI_COMM_WORLD,&status); //Receiving the signal of the provious process
	fp = fopen("project.txt","a");
	for (i=0;i<I;i++)
	{fprintf(fp, "%d,%d,%d,", x[i],y[i],c[i]);}
	fclose(fp);
	if(p!=P-1){
	MPI_Send (&s,1,MPI_INT,p+1,0,MPI_COMM_WORLD);}} //Sending a signal to the next process
*/
free(x); 
free(y);
free(cx); 
free(cy);
free(cxr); 
free(cyr);
free(c);
MPI_Finalize();
exit(0);}