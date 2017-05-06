#include <math.h>
#include <stdlib.h>
#define G 6.67e-11

int *bubblesort(float *sumatoria,int *index,int size)
{
	int i,j;
	for(i=0;i<size;i++)
	{
		for(j=0;j<size-1;j++)
		{
			if(sumatoria[j]<sumatoria[j+1])
			{
				 
				float aux1=sumatoria[j];
				sumatoria[j]=sumatoria[j+1];
				sumatoria[j+1]=aux1;
				///////ORDEN IDX//////////
				int aux2=index[j];
				index[j]=index[j+1];
				index[j+1]=aux2;
			}			
		}		
	}
	free(sumatoria);
	return index;
}

void clearSpacePoints(float *points,int sizeR,int sizeC)
{
	int j;
	int idy = 0;
	int idx=0;

	for(j=0;j<sizeR;j++)
	{
		if(points[j*sizeC] != -2)
		{
			for (idx = 0; idx<sizeC; idx++)
			{
				points[idy*sizeC+idx] = points[j*sizeC+idx];
			}
			idy++;
		}
	}
}
 
 
/*DELETE MASS WHEN THIS IS -1 AND RESIZE VECTOR*/
void clearMass(float *vector,int size)
{
	int i,idx=0;
	for(i=0;i<size;i++)
	{
		if(vector[i]!=-1)
		{
			vector[idx]=vector[i];
			idx++;
		}	
	}
}

int *SortSumAllColumns(float *matrix,int size)
{
	float *vector;
	vector=(float *)malloc(size*sizeof(float));

	int *indexSort;
	indexSort = (int *)malloc(size * sizeof(int));

	int i,j;
	for(i=0;i<size;i++)
	{
		vector[i]=0;
		for(j=0;j<size;j++)
		{
			vector[i]+=matrix[j*size+i];						//SUMAS LAS COLUMANAS DE LA MATRIZ Y AGREGA EN UN VECTOR
		}
		indexSort[i]=i;
	}
	return bubblesort(vector,indexSort,size);
}



void calAllForcesGravitational_CPU(float *F,float *dist,float *m,float *x,int SIZE,int numVars)
{
	int i,j,k;
	float dSum;
	for(i=0;i<SIZE;i++)
		{
			for(j=0;j<SIZE;j++)
			{
				if(i != j)
				{
					dSum = 0;
					for(k=0;k<numVars;k++)
					{
						dSum += pow( x[i*numVars+k] - x[j*numVars+k],2);
					}
					dist[i * SIZE + j] = sqrt(dSum);
					F[ i * SIZE + j ] = (G*m[i]*m[j]) / dSum;
				
				}
				else
				{
					dist[i*SIZE + j] =3;
				}
			}
		}

}



