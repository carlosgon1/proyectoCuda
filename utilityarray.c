#include <stdio.h>
void printArray(int *array, int size)
{
	int i;
	printf("\n ************************ARRAY*************************\n");
	for (i = 0; i<size; i++)
	{
		printf("%i ", array[i]);
	}
	printf("\n");
	printf("\n");
}

void printArrayfloat(float *array, int size)
{
	printf("\n ************************ARRAY*************************\n");
	int i;
	for (i = 0; i<size; i++)
	{
		printf("%.3e 	", array[i]);
	}
	printf("\n");
	printf("\n");
}


void printMatrix(float *matrix, int sizeRow, int sizeCol)
{
	int i, j;
	printf("\n ************************MATRIX**************************\n");
	for (i = 0; i<sizeRow; i++)
	{
		for (j = 0; j<sizeCol; j++)
		{
			printf("%0.3e ", matrix[i*sizeCol+j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("\n");

}

void initializateArray(float *array,int size,int val)
{
	int i;
	for(i=0;i<size;i++)
	{
		array[i]=1*val;
	}
}
