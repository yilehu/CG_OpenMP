#include <stdio.h>

void MatrixDefinition_Banded(double **Matrix,int m,int Bandwidth)
{
	int i,j,k,tid;
	for(i=0;i<m;i++)
	{
		for(j=0;j<Bandwidth;j++)
		{
			if(i>Bandwidth-j-2)
			{
				Matrix[i][j] = j+1;
				
			}
			if(i<m-(Bandwidth-1-j))
			{
				Matrix[i][2*Bandwidth-2 - j] = j+1;
			}
		}
	}
}

void MatrixMultiply_Banded(double **A,double *x,double *b,int m,int n,int Bandwidth)
{
	int i,j;
	for(i=0;i<Bandwidth-1;i++)
	{
		b[i] = 0.0;
		for(j=Bandwidth-1-i;j<n;j++)
		{
			b[i] += A[i][j]*x[i-(Bandwidth-1)+j];
		}
	}
	for(i=Bandwidth-1;i<m-Bandwidth+1;i++)
	{
		b[i] = 0.0;
		for(j=0;j<n;j++)
		{
			b[i] += A[i][j]*x[i-(Bandwidth-1)+j];
			
		}
	}
	for(i=m-Bandwidth+1;i<m;i++)
	{
		b[i] = 0.0;
		for(j=0;j<Bandwidth-1+m-i;j++)
		{
			b[i] += A[i][j]*x[i-(Bandwidth-1)+j];
		}
	}
}

double Dotproduct(double *a,double *b,int n)
{
	double temp = 0.0;
	for(int i=0;i<n;i++)
	{
		temp += a[i]*b[i];
	}
	return temp;
}

void MatrixMultiply_Parallel(int tid,int i,double **A,double *x,double *b,int m,int Bandwidth)
{
	int j,j_start,j_end;
	if(i>=0 && i<Bandwidth-1)
	{
		j_start = Bandwidth-1-i; j_end = 2*Bandwidth-1;
	}
	else if(i>=Bandwidth-1 && i<m-Bandwidth+1)
	{
		j_start = 0; j_end = 2*Bandwidth-1;
	}
	else if(i>=m-Bandwidth+1 && i<m)
	{
		j_start = 0; j_end = Bandwidth-1+m-i;
	}
	b[i] = 0.0;
	for(j=j_start;j<j_end;j++)
	{
		b[i] += A[i][j]*x[i-(Bandwidth-1)+j];
	}
}