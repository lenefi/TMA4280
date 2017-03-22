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

void transpose(real **bt, real **b, int *nrows, int *displ, int size, int rank,size_t m);
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

size_t m=4;
real **b = mk_2D_array(m, m, false);
real **bt = mk_2D_array(m, m, false);
int i,j;

//Initalise MPI
int rank, size;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&size);

//distribute rows
    int part = m/size;
    int *nrows = calloc(size, sizeof(int));
    for (size_t i = 0; i < size; i++){
        nrows[i] = part;
printf("nrows= %d ",nrows[i]);
}printf("\n");



//'leftover' columns
    int leftover = m % size;
   for (size_t i = 1;  i <= leftover; i++)
      nrows[size - i]++;


//Calculate displacment vector for call to MPI_Alltoallv
    int *displ = calloc(size+1, sizeof(int));
    displ[0] = 0;
    for (size_t i = 1; i < size+1; i++){
    displ[i] = displ[i-1] + nrows[i-1];
 	printf("disp= %d ",displ[i]);}
printf("\n");

//Print matrix b
printf("b=\n");

for(i=0; i<m; i++){
    for(j=0; j<m; j++){
        b[i][j]=(i)*m + (j+1);
        printf("%f ",b[i][j]); }
    printf("\n");	
}


//transpose(bt, b, nrows, displ, size, rank);

//Print matrix bt
printf("bt= \n");	
for(i=0;i<m;i++){
    for(j=0;j<m;j++){
	printf("%f ",bt[i][j]); 
	}
	printf("\n");	
}

MPI_Finalize();
    return 0;
}

//Transpose function
void transpose(real **bt, real **b, int *nrows, int *displ, int size, int rank,size_t m)
{
int j;
size_t i;
	while ( i< m*m){
    		for (j < 0; j < size-1; j++){
    	
 	MPI_Alltoallv(b[i],nrows,displ,MPI_DOUBLE,bt[i],nrows,displ,MPI_DOUBLE,MPI_COMM_WORLD);
					}
		i+=displ[j];
				}
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


