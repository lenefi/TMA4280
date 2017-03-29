/**
 * C program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>


#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

// Function prototypes

void transpose(real **bt, real **b, int *nrows, int *displ, int size, int rank);
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);


int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage:\n");
        printf("  poisson n\n\n");
        printf("Arguments:\n");
        printf("  n: the problem size (must be a power of 2)\n");
    }

    // The number of grid points in each direction is n+1
    // The number of degrees of freedom in each direction is n-1
    int n = atoi(argv[1]);

//Initalise MPI
int rank, size;

MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&size);

//global size
size_t M = n;
// local size
size_t m = M/size;
int part = m;
// local row-wise block
real **b = mk_2D_array(m, M, false);
real **bt = mk_2D_array(m, M, false);
int i,j;


//distribute rows
    int *nrows = calloc(size, sizeof(int));
    for (size_t i = 0; i < size; i++){
        nrows[i] = part;
}



//'leftover' columns
    int leftover = M % size;
   printf("leftover %u\n", leftover);
   for (size_t i = 1;  i <= leftover; i++)
   {
      nrows[size - i]++;
   }

//Calculate displacment vector for call to MPI_Alltoallv
    int *displ = calloc(size+1, sizeof(int));
    displ[0] = 0;
    for (size_t i = 1; i < size+1; i++){
      printf("displ= %d\n",displ[i-1]);
      printf("nrows= %d\n",nrows[i-1]);
      displ[i] = displ[i-1] + nrows[i-1];  
 	}

//Print matrix b
printf("b=\n");

for(i=0; i<m; i++){
    for(j=0; j<M; j++){
        b[i][j]= (displ[rank]+i)*M + (j+1);
        printf("%f ",b[i][j]); }
    printf("\n");	
}

transpose(bt, b, nrows, displ, size, rank);

MPI_Barrier(MPI_COMM_WORLD);

//Print matrix bt
printf("%u : bt@%0x= \n", rank, bt);	
for(i=0;i<m;i++){
    for(j=0;j<M;j++){
	printf("%f ",bt[i][j]); 
	}
	printf("\n");	
}

MPI_Finalize();
    return 0;
}

//Transpose function
void transpose(real **bt, real **b, int *nrows, int *displ, int size, int rank)
{
printf("on rank %u : nrows = %u; displ = %u \n", rank, nrows[rank], displ[rank]);
MPI_Alltoallv(b[0],nrows,displ,MPI_DOUBLE,bt[0],nrows,displ,MPI_DOUBLE,MPI_COMM_WORLD);				
}

//Generate array-function
real *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

real **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    real **ret = (real **)malloc(n1 * sizeof(real *));

    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }

    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}


