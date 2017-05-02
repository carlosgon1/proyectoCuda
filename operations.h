#ifndef _OPERATIONS_H_
#define _OPERATIONS_H_
#define G 6.67e-11

extern "C" void calAllForcesGravitational_CPU(float *F,float *dist,float *m,float *x,int SIZE,int numVars);
extern "C" void clearSpacePoints(float *matrix,int sizeR,int sizeC);
extern "C" void clearMass(float *vector,int size);
extern "C" int *SortSumAllColumns(float *matrix,int size);
int *bubblesort(float *sumatoria,int *index,int size);

#endif 
