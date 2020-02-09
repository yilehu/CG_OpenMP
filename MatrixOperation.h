#ifndef _MATRIXOPERATION_ 
#define _MATRIXOPERATION_ 

void MatrixDefinition_Banded(double **Matrix,int m,int Bandwidth);

void MatrixMultiply_Banded(double **A,double *x,double *b,int m,int n,int Bandwidth);

double Dotproduct(double *a,double *b,int n);

void MatrixMultiply_Parallel(int tid,int i,double **A,double *x,double *b,int m,int Bandwidth);

#endif