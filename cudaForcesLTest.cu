#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cuda.h>
#include "utilityarray.h"
#include "operations.h" 
#define R 5636
#define C 1000

typedef struct
{
	float minD;				//GUARDA EL MINIMO DEL VECTOR.
	int md;					//GUARDA EL INDICE DONDE ENCONTRO EL MINIMO.
}Particules;

void calAllForcesGravitational_GPU(float *h_Forces,float *h_dist,float *h_mass,float *h_spacePoints,int M,int N);
void checkCUDAError(const char *msg);
__global__ void kernel(float *F,float *points,float *dist,float *mass,int M,int numVar);
Particules min(float *vector,int SIZE);


float x1[R][C]={ 	
	
int main()
{  
    int i,k,j;
    int J=R;
    size_t sizeM = J * sizeof(float);     
    size_t sizeN= C* sizeof(float); 
   	int SIZE=J;
	int INTERACTIONS_EXIST = 1;
	int iters=0;
	float radius = 0.5;
	float modifier = 0.9;
    float d, distp, t;
	int op1, op2, a, b;
 	
 	// DECLARACION E INICIALIZACION DE VARIABLES DEL TEMPORIZADOR EN EL GPU
  	float time; 
  	cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    //////////////////////////////////////

    float *x;
    x=(float*)malloc(sizeM*sizeN);

    float *F;
    F = (float*)malloc(sizeM*sizeM);  

    float *dist;
    dist= (float*)malloc(sizeM*sizeM);  

    float *m;
    m= (float*)malloc(sizeM); 

    int *index;
    initializateArray(m,J,100);
    //COPY ARRAY GLOBAL TO ARRAY LOCAL
    for(i=0;i<R;i++)
    {
        for(j=0;j<C;j++)
        {
            x[i*C+j]=x1[i][j];
        }
    }
  
	while(INTERACTIONS_EXIST ==1)
	{
		 
		INTERACTIONS_EXIST = 0;
		iters = iters + 1;
		initializateArray(F, SIZE*SIZE, 0);
		initializateArray(dist, SIZE*SIZE, 0);

		/***Calculate all interacting forces in the system***/
		calAllForcesGravitational_GPU(F,dist,m,x,SIZE,C);
		/*cudaEventRecord(start, 0);
		calAllForcesGravitational_CPU(F,dist,m,x,SIZE,C);
		cudaEventRecord(stop, 0);    
    	cudaEventSynchronize(stop); 
    	cudaEventElapsedTime( &time, start, stop );
    	printf("TIEMPO DE EJECUCION EN CPU: %f ms\n",time);
    	getchar();
		/******************Reorder particule Data*******************/
		//printMatrix(F,SIZE,SIZE,"FUERZAS GRAVITACIONALES CALCULADAS");
		index=SortSumAllColumns(F,SIZE);
		//printArray(index,SIZE,"INDICE DE PARTICULAS ORDENADAS");
		/******************Unificate Particules*********************/
		int newSize=SIZE;
		for(i=0;i<SIZE;i++)
		{
			int idx=index[i];
			Particules p=min(dist+(idx*SIZE),SIZE);
			op1=idx;
			op2=p.md;
			if(m[op1]!=-1 && m[op2]!=-1 && op1!=op2 && p.minD<radius)
			{
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
				/******************move particulas according to masses************/
				float dSum;
				for(k=0;k<C;k++)
				{
					dSum= dSum + pow((x[a*C+k]-x[b*C + k]),2);
				}
				d=sqrt(dSum);
				distp=(m[b] / (m[a]+m[b]))*d;
				t=distp/d;
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

		/******************CLEAN UP ELEMENT FROM MASS AND SPACE POINTS*************/
		clearMass(m,SIZE);
		clearSpacePoints(x,SIZE,C);

		//printArrayfloat(m,newSize,"VALORES DE LAS MASAS");
		//printMatrix(x,newSize,C,"PUNTOS EN EL ESPACIOS");
		radius=radius*modifier;
		SIZE=newSize;
	}
	//printMatrix(x,SIZE,C,"RESULTADO FINAL");
	cudaEventDestroy( start ); // GC DEL TEMPORIZADOR
    cudaEventDestroy( stop );  // GC DEL TEMPORIZADOR
	free(m);
	free(x);
	free(F);
	free(dist);
	 
	return 0;
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




void calAllForcesGravitational_GPU(float *h_Forces,float *h_dist,float *h_mass,float *h_spacePoints,int M,int N)
{
	size_t sizeN = N * sizeof(float); 			 
    size_t sizeM = M * sizeof(float); 
    // GPU'S TIMER VARIABLES
  	float time; 
  	cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    //////////////////////////////////////

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
    /*********************APPLY TIMER***********************************/
    cudaEventRecord(start, 0); 
    kernel<<< grid,block>>>(d_Forces,d_spacePoints,d_dist,d_mass,M,N);
	cudaEventRecord(stop, 0);    
    cudaEventSynchronize(stop); 
    cudaEventElapsedTime( &time, start, stop );
	cudaThreadSynchronize();
	printf("TIEMPO DE EJECUCION EN GPU: %f ms\n",time);
    getchar();

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