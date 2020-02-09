#ifndef _INITIALIZATION_ 
#define _INITIALIZATION_ 

void InitializeMatrix(double **Matrix,int m,int n,double InitialValue);

void InitializeArray(double *Array,int m,double InitialValue);

void InitializeArray_Parallel(int tid,int i,double *Array,double InitialValue);

#endif