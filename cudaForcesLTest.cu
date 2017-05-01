#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cuda.h>

#define G 6.67e-11
#define R 30
#define C 3
typedef struct
{
	float minD;				//GUARDA EL MINIMO DEL VECTOR.
	int md;					//GUARDA EL INDICE DONDE ENCONTRO EL MINIMO.
}Particules;
void checkCUDAError(const char *msg);
__global__ void kernel(float *F,float *points,float *dist,float *mass,int M,int numVar);
void calAllForcesGravitational_GPU(float *h_Forces,float *h_dist,float *h_mass,float *h_spacePoints,int M,int N);
void calAllForcesGravitational_CPU(float *F,float *dist,float *m,float *x,int SIZE,int numVars);
void clearSpacePoints(float *matrix,int sizeR,int sizeC);
void clearMass(float *vector,int size);
int *bubblesort(float *sumatoria,int *index,int size);
 int *SortSumAllColumns(float *matrix,int size);
Particules min(float *vector,int SIZE);

void printArray(int *array,int size);
void printArrayfloat(float *array,int size);
void printMatrix(float *matrix, int size, int sizeX);
 

void initializateArray(float *array,int size,int mult);
float x1[R][C]={   {0.1506, 0.0898, 0.0229},
					{0.2145, 0.1316, 0.2570},
					{0.1903, 0.1984, 0.2603},
					{0.2757, 0.0004, 0.2116},
					{0.0309, 0.1581, 0.0190},
					{0.2663, 0.2360, 0.2723},
					{0.1939, 0.0978, 0.1605},
					{0.2810, 0.2015, 0.0249},
					{0.1758, 0.0134, 0.1723},
					{0.1659, 0.2623, 0.1682},
					{0.5780, 0.5930, 0.5035},
					{0.5192, 0.3544, 0.3773},
					{0.4037, 0.5772, 0.6048},
					{0.4352, 0.4619, 0.4602},
					{0.5692, 0.4652, 0.5358},
					{0.4934, 0.5750, 0.4121},
					{0.5949, 0.5094, 0.5401},
					{0.5965, 0.5426, 0.5782},
				    {0.5244, 0.5962, 0.3561},
				    {0.4081, 0.4082, 0.5684},
				    {0.9133, 0.8158, 0.8304},
				    {0.9095, 1.0000, 0.7354},
				    {0.8366, 0.8548, 0.7556},
				    {0.8147, 0.8896, 0.9289},
				    {0.9720, 0.8373, 0.9794},
				    {0.8591, 0.8378, 0.7319},
		   		    {0.9793, 0.9635, 0.7751},
					{0.7467, 0.8534, 0.8549}, 
				    {0.8930, 0.9876, 0.9589},
				    {0.8329, 0.9730, 0.7375}
			    };
 
 /*
int main()
{	
	int SIZE = R;					 
	int numVars = C;
	size_t sizeM = SIZE * sizeof(float);
	size_t sizeN = numVars * sizeof(float);
	size_t sizeMxM = SIZE * SIZE* sizeof(float);
	int i,j,k;
	int INTERACTIONS_EXIST = 1;
	int iters=0;
	
	float radius = 0.5;
	float modifier = 0.9;
	float mass = 100;

	float *x;
	x=(float*)malloc(sizeM*sizeN);

	 /////////////////////////
	for(i=0;i<SIZE;i++)
	{
		for(j=0;j<numVars;j++)
		{
			x[i*numVars+j]=x1[i][j];
		}
	}
	//CREA ARREGLO  CON TODOS LOS PESOS DE LAS MASAS
	float *m;
	m=(float*)malloc(sizeM);						
	initializateArray(m,SIZE,mass);

	int *index;

	float *F;
	F = (float*)malloc(sizeMxM);
	
	float *dist;
	dist = (float*)malloc(sizeMxM);

	float d, distp, t;
	int op1, op2, a, b;
	
	while(INTERACTIONS_EXIST ==1)
	{
		 
		INTERACTIONS_EXIST = 0;
		iters = iters + 1;
		initializateArray(F, SIZE*SIZE, 0);
		initializateArray(dist, SIZE*SIZE, 0);
		//Calculate all interacting forces in the system
		//calAllForcesGravitational_CPU(F,dist,m,x,SIZE,numVars);
   		calAllForcesGravitational_GPU(F,dist,m,x,SIZE,numVars);
		//reorder particule Data.
		printMatrix(F,SIZE,SIZE);
		getchar();
		index=SortSumAllColumns(F,SIZE);
		printArray(index,SIZE);
		//Unificate Particules
		int newSize=SIZE;
		for(i=0;i<SIZE;i++)
		{
			int idx=index[i];
			Particules p=min(dist+(idx*SIZE),SIZE);
			op1=idx;
			op2=p.md;
			//printf("i=%d 	m1=%.4f m2=%4f  \n",i,m[op1],m[op2]);
			if(m[op1]!=-1 && m[op2]!=-1 && op1!=op2 && p.minD<radius)
			{
				//if(ii==2)
				printf("i=%d 	m1=%.4f m2=%4f  op1=%d  op2=%d \n",i+1,m[op1],m[op2],op1,op2);
				//printf("i=%d  ",op1);  
				INTERACTIONS_EXIST=1;
				if(m[op1] >= m[op2])
				{
					a=op1;
					b=op2;
				}
				else
				{
					a=op2;
					b=op1;
				}
				m[a]=m[a]+m[b];
				//%move particulas according to masses
				float dSum;
		        printf("\nma=%.2f \n",m[a]);
				for(k=0;k<numVars;k++)
				{
					dSum= dSum + pow((x[a*numVars+k]-x[b*numVars + k]),2);
				}
				d=sqrt(dSum);
				distp=(m[b] / (m[a]+m[b]))*d;
				t=distp/d;

				//sum2Vector_MulxT
				for (j = 0; j<numVars; j++)
				{
					x[a*numVars+j] = x[a*numVars + j] + t*(x[b*numVars + j] - x[a*numVars + j]);
				}

				m[b]=-1;
				x[b*numVars] = -2;
				newSize--;
			}
			 
		}
		free(index);

		//CLEAN UP ELEMENT FROM MASS AND SPACE POINTS
		clearMass(m,SIZE);
		clearSpacePoints(x,SIZE,numVars);

		printf("valores de masas\n",SIZE);
		printArrayfloat(m,newSize);
		printMatrix(x,newSize,numVars);
		radius=radius*modifier;
		SIZE=newSize;
	}
//printArray(index,SIZE);
	//printMatrix(F,SIZE);
	//printf("\n\n\n");
	//printMatrix(dist,SIZE);
	free(m);
	free(x);
	free(F);
	free(dist);
	 
	return 0;
}
*/
 
int main()
{  
    int i,k,j;
    int J=R;
    size_t sizeM = J * sizeof(float);           // para representar tamaños cuando se usa el sizeof()
     size_t sizeN= C* sizeof(float); 
   int SIZE=J;

  
	int INTERACTIONS_EXIST = 1;
	int iters=0;
	
	float radius = 0.5;
	float modifier = 0.9;
	 


    float *x;
    x=(float*)malloc(sizeM*sizeN);

     /////////////////////////
    for(i=0;i<R;i++)
    {
        for(j=0;j<C;j++)
        {
            x[i*C+j]=x1[i][j];
        }
    }


    float *F;
    F = (float*)malloc(sizeM*sizeM);  

    float *dist;
    dist= (float*)malloc(sizeM*sizeM);  

    float *m;
    m= (float*)malloc(sizeM); 

    int *index;

    initializateArray(m,J,100);
    /////////////////////////////////////

    float d, distp, t;
	int op1, op2, a, b;
	
	while(INTERACTIONS_EXIST ==1)
	{
		 
		INTERACTIONS_EXIST = 0;
		iters = iters + 1;
		initializateArray(F, SIZE*SIZE, 0);
		initializateArray(dist, SIZE*SIZE, 0);
		//Calculate all interacting forces in the system
		//calAllForcesGravitational_CPU(F,dist,m,x,SIZE,numVars);
   		calAllForcesGravitational_GPU(F,dist,m,x,SIZE,C);
		//reorder particule Data.
		printMatrix(F,SIZE,SIZE);
		getchar();
		index=SortSumAllColumns(F,SIZE);
		printArray(index,SIZE);
		//Unificate Particules
		int newSize=SIZE;
		for(i=0;i<SIZE;i++)
		{
			int idx=index[i];
			Particules p=min(dist+(idx*SIZE),SIZE);
			op1=idx;
			op2=p.md;
			//printf("i=%d 	m1=%.4f m2=%4f  \n",i,m[op1],m[op2]);
			if(m[op1]!=-1 && m[op2]!=-1 && op1!=op2 && p.minD<radius)
			{
				//if(ii==2)
				printf("i=%d 	m1=%.4f m2=%4f  op1=%d  op2=%d \n",i+1,m[op1],m[op2],op1,op2);
				//printf("i=%d  ",op1);  
				INTERACTIONS_EXIST=1;
				if(m[op1] >= m[op2])
				{
					a=op1;
					b=op2;
				}
				else
				{
					a=op2;
					b=op1;
				}
				m[a]=m[a]+m[b];
				//%move particulas according to masses
				float dSum;
		        printf("\nma=%.2f \n",m[a]);
				for(k=0;k<C;k++)
				{
					dSum= dSum + pow((x[a*C+k]-x[b*C + k]),2);
				}
				d=sqrt(dSum);
				distp=(m[b] / (m[a]+m[b]))*d;
				t=distp/d;

				//sum2Vector_MulxT
				for (j = 0; j<C; j++)
				{
					x[a*C+j] = x[a*C + j] + t*(x[b*C + j] - x[a*C + j]);
				}

				m[b]=-1;
				x[b*C] = -2;
				newSize--;
			}
			 
		}
		free(index);

		//CLEAN UP ELEMENT FROM MASS AND SPACE POINTS
		clearMass(m,SIZE);
		clearSpacePoints(x,SIZE,C);

		printf("valores de masas\n",SIZE);
		printArrayfloat(m,newSize);
		printMatrix(x,newSize,C);
		radius=radius*modifier;
		SIZE=newSize;
	}
//printArray(index,SIZE);
	//printMatrix(F,SIZE);
	//printf("\n\n\n");
	//printMatrix(dist,SIZE);
	free(m);
	free(x);
	free(F);
	free(dist);
	 
	return 0;
}


//////////////////////////////////////////
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


 

void initializateArray(float *array,int size,int mult)
{
	int i;
	for(i=0;i<size;i++)
	{
		array[i]=1*mult;
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

Particules min(float *vector,int size)
 {
 	Particules p;
 	p.minD=0;
 	p.md=0;
 	int i;
 	for(i=0;i<size;i++)
	{
		if(p.minD==0)
		{
			p.minD=vector[i];
			p.md=i;
		}
		else if(vector[i]<p.minD)
		{
			p.minD=vector[i];	
			p.md=i;
		}
	}
	return p;
 }




/*
***********************************
*ONLY PRINTS OF MATRIX OR VECTORS.
***********************************
*/
 
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
		printf("%.4e 	", array[i]);
	}
	printf("\n");
	printf("\n");
}


void printMatrix(float *matrix, int size, int sizeX)
{
	int i, j;
	printf("\n ************************ARRY LIKE MATRIX*************************\n");
	for (i = 0; i<size; i++)
	{
		for (j = 0; j<sizeX; j++)
		{
			printf("%0.3e ", matrix[i*sizeX+j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("\n");

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


void calAllForcesGravitational_GPU(float *h_Forces,float *h_dist,float *h_mass,float *h_spacePoints,int M,int N)
{
	size_t sizeN = N * sizeof(float); 			 
    size_t sizeM = M * sizeof(float); 

    /*CONFIGURATION FOR KERNEL FOR TO WORK WHIT MATRIX 5632 x 5632 */
    dim3 grid(256,256);
    dim3 block(22,22);

    /*CREATE MATRIX WHERE SAVE ALL VALUES OF ALL FORCES CALCULATED*/
	float *d_Forces;
    cudaMalloc((void **)&d_Forces,sizeM*sizeM);
    cudaMemset(d_Forces,0,sizeM*sizeM); 	

    /*CREATE MATRIX WHERE SAVE ALL VALUES OF THE DISTANCES*/
	float *d_dist;
    cudaMalloc((void **)&d_dist,sizeM*sizeM);
    cudaMemset(d_dist,0,sizeM*sizeM); 					

	/*COPY VALUES OF SPACEPOINTS OF HOST TO DEVICE*/
    float *d_spacePoints;
    cudaMalloc((void **)&d_spacePoints,sizeM*sizeN);
    cudaMemset(d_spacePoints,0,sizeM*sizeN); 
    cudaMemcpy(d_spacePoints,h_spacePoints,sizeM*sizeN,cudaMemcpyHostToDevice);

    /*COPY VALUES OF MASS OF HOST TO DEVICE*/
    float *d_mass;
    cudaMalloc((void **)&d_mass,sizeM);
    cudaMemcpy(d_mass,h_mass,sizeM,cudaMemcpyHostToDevice);
    printf("desde gpu where N=%d y M=%d\n",M,N);
	printMatrix(h_spacePoints,M,N); 	
	getchar();
    kernel<<< grid,block>>>(d_Forces,d_spacePoints,d_dist,d_mass,M,N);
	cudaThreadSynchronize();
    /*RECOVERY FORCES AND DISTANCE TO HOST*/
    cudaMemcpy(h_Forces,d_Forces,sizeM*sizeM,cudaMemcpyDeviceToHost);
    checkCUDAError("cudaMemcpy");
    cudaMemcpy(h_dist,d_dist,sizeM*sizeM,cudaMemcpyDeviceToHost);
    checkCUDAError("cudaMemcpy");

    cudaFree(d_Forces);
    cudaFree(d_mass);
    cudaFree(d_dist);
    cudaFree(d_spacePoints);
}


__global__ void kernel(float *F,float *points,float *dist,float *mass,int M,int numVar)
{
	int i;
	float rest,dSum;
	int tidx = (blockDim.x * blockIdx.x) + threadIdx.x;
	int tidy = (blockDim.y * blockIdx.y) + threadIdx.y;
    if(tidx <M && tidy<M)
    {
    	if( tidx != tidy )
     	{
     		dSum=0;
     		for( i=0; i < numVar; i++ )
     		{
    			rest = (points[(tidx * numVar) + i]-points[(tidy * numVar) + i]);
    			dSum= dSum + (rest*rest);//al cuadrado
     		}
     		dist[tidx*M+tidy] = sqrt(dSum);
     		F[tidx*M+tidy] = (G * mass[tidx] * mass[tidy]) / dSum;
     	} 
     	else
     	{
     		dist[tidx*M+tidy]=3;
     	}
    }
}

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
    }                         
}