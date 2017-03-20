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

size_t m=6;
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
}



//'leftover' columns
    int leftover = m % size;
   for (size_t i = 1;  i <= leftover; i++)
      nrows[size - i]++;


//Calculate displacment vector for call to MPI_Alltoallv
    int *displ = calloc(size+1, sizeof(int));
    displ[0] = 0;
    for (size_t i = 1; i < size+1; i++)
    displ[i] = displ[i-1] + nrows[i-1];
 

//Print matrix b
printf("b=\n");

for(i=0; i<m; i++){
    for(j=0; j<m; j++){
        b[i][j]=(i)*m + (j+1);
        printf("%f ",b[i][j]); }
    printf("\n");	
}


transpose(bt, b, nrows, displ, size, rank);

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
//LARS, HAR LAGT INN FUNKSJON HER!
//Transpose function
void transpose(real **bt, real **b, int *nrows, int *displ, int size, int rank)
{
    int * temp_displ = calloc(size+1, sizeof(int));
    int * temp_nrows = calloc(size, sizeof(int));
    double *temp = calloc(nrows[rank]*size, sizeof(double));

temp_displ[0] = 0;
    for (int i = 1; i < size+1; i++){
    temp_displ[i] = temp_displ[i-1] + nrows[rank];
    temp_nrows[i-1] = nrows[rank];
}
    
for (size_t i; i < nrows[size-1]; i++){
//If number of row is smaller than rank of process, store into temp adress    
    if (i <nrows[rank]){
 	MPI_Alltoallv(b[i],nrows,displ,MPI_DOUBLE,temp,temp_nrows,temp_displ,MPI_DOUBLE,MPI_COMM_WORLD);

}
//If not smaller than process, store into temp adress    
    else {
 	MPI_Alltoallv(b[i-1], nrows, displ, MPI_DOUBLE, temp,temp_nrows,temp_displ,MPI_DOUBLE,MPI_COMM_WORLD);

}
//Insert to adress in bt, transposed column,row.
   for (size_t r=0; r<size; r++){
	for (size_t c=0; c<nrows[rank]; c++){
	    if (displ[r]+1 < displ[r+1]){
		bt[c][displ[r]+1] = temp[temp_displ[r]+c];
			}	
		}
	}
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


