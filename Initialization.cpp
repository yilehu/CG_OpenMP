#include <stdio.h>

void InitializeArray(double *Array,int m,double InitialValue)
{
	int i;
	for(i=0;i<m;i++)
	{
		Array[i] = InitialValue;
	}
}

void InitializeMatrix(double **Matrix,int m,int n,double InitialValue)
{
	int i,j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			Matrix[i][j] = InitialValue;
		}
	}
}

void InitializeArray_Parallel(int tid,int i,double *Array,double InitialValue)
{
	Array[i] = InitialValue;
	//printf("Thread = %d, i = %d\n",tid,i);
}