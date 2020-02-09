#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Initialization.h"
#include "PrintToFile.h"
#include "MatrixOperation.h"
#include <time.h>
#include <omp.h>

int main(int argc,char *argv[])
{
	int i,j,k;
	double START_WTIME,END_WTIME;
	double START_CLOCK,END_CLOCK;
	double Time,WTime;
	double Iter_Running_Time,Total_Running_Time;

	char *Directory1,*Directory2;
	Directory1 = "Array.txt";
	Directory2 = "Matrix.txt";

	//************ OpenMP Setup ***********//
	int Num_Threads,Num_Threads_Max;
	#pragma omp parallel
	{
		Num_Threads_Max = omp_get_num_threads();
	}
	printf("Maximum number of threads = %d\n",Num_Threads_Max);

	Num_Threads = 12;
	omp_set_num_threads(Num_Threads);
	printf("Number of threads chosen = %d\n",Num_Threads);
	printf("\n");

	//************ Memory Allocations ***********//
	double *x,*b,*r,*r_new,*p,*Ax;
	double **Matrix;
	double Error,Error0 = 1.0e-6;
	double alpha,beta,Sum1,Sum2;
	int n = 20000;
	int Bandwidth = 5;

	Matrix = (double**)malloc(n*sizeof(double*));
	for(i=0;i<n;i++)
	{
		*(Matrix+i) = (double*)malloc((2*Bandwidth-1)*sizeof(double));
	}

	x = (double*)malloc(n*sizeof(double));
	b = (double*)malloc(n*sizeof(double));
	r = (double*)malloc(n*sizeof(double));
	r_new = (double*)malloc(n*sizeof(double));
	p = (double*)malloc(n*sizeof(double));
	Ax = (double*)malloc(n*sizeof(double));

	InitializeMatrix(Matrix,n,2*Bandwidth-1,0.0);
	InitializeArray(x,n,0.0);
	InitializeArray(b,n,1.0);
	MatrixDefinition_Banded(Matrix,n,Bandwidth);

	//************ CG ***********//
	//Initialization//
	Total_Running_Time = 0.0;
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		#pragma omp for schedule(static)
		for(int i=0;i<n;i++)
		{InitializeArray_Parallel(tid,i,x,0.0);}

		#pragma omp for schedule(static) reduction(+:Error)
		for(int i=0;i<n;i++)
		{
			MatrixMultiply_Parallel(tid,i,Matrix,x,Ax,n,Bandwidth);
			r[i] = b[i] - Ax[i];
			p[i] = r[i];
			Error += r[i] * r[i];
		}
	}
	k = 0;
	Error = sqrt(Error);
	//Iteration//
	while(Error>Error0)
	{
		START_CLOCK = clock();
		Sum1 = 0.0;
		Sum2 = 0.0;
		Error = 0.0;
		#pragma omp parallel default(none) shared(Matrix,p,Ax,r,r_new,x,b,n,Bandwidth,alpha,beta,Error,Sum1,Sum2)
		{
			int tid = omp_get_thread_num();

			#pragma omp for schedule(static) reduction(+:Sum1,Sum2)
			for(int i=0;i<n;i++)
			{MatrixMultiply_Parallel(tid,i,Matrix,p,Ax,n,Bandwidth); Sum1 += r[i]*r[i]; Sum2 += p[i]*Ax[i];}

			#pragma omp master
			{alpha = Sum1/Sum2;}
			#pragma omp barrier

			#pragma omp for schedule(static) reduction(+:Error)
			for(int i=0;i<n;i++)
			{x[i] += alpha*p[i]; r_new[i] = r[i]-alpha*Ax[i]; Error += r_new[i]*r_new[i];}

			#pragma omp master
			{beta = Error/Sum1; Error = sqrt(Error);}
			#pragma omp barrier

			#pragma omp for schedule(static)
			for(int i=0;i<n;i++)
			{p[i] = r_new[i]+beta*p[i]; r[i] = r_new[i];}
		}
		END_CLOCK = clock();
		Iter_Running_Time = (double)(END_CLOCK - START_CLOCK)/CLOCKS_PER_SEC;
		Total_Running_Time += Iter_Running_Time;
		k++;
		if(k % 200 == 0) {printf("Iteration number = %d, Error = %12E, Iteration time = %12.6lf, Total time = %12.6lf\n",k,Error,Iter_Running_Time,Total_Running_Time);}
	}
	printf("Iteration number = %d, Error = %12E, Total time = %12.6lf\n",k,Error,Total_Running_Time);
	printf("\n");
	system("pause");

	return 0;
}
