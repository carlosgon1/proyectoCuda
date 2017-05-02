#ifndef _UTILITYARRAY_H_
#define _UTILITYARRAY_H_
extern "C" void printArray(int *array,int size);				//PRINT ARRAY TYPE INT OF size "SIZE"
extern "C" void printArrayfloat(float *array,int size);			//PRINT ARRAY FLOAT TYPE INT OF size "SIZE"
extern "C" void printMatrix(float *matrix, int sizeRow, int sizeCol);	//PRINT MATRIX TYPE FLOAT OF size "SIZEROW x SIZECOL
extern "C" void initializateArray(float *array,int size,int val);		//INITIALIZE ARRAY TO MULT
#endif 
