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
#include "mpi.h"

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
//void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

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
	int m = n - 1;
	int nn = 4 * n;
	real h = 1.0 / n;
	//Used for MPI partitioning
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	int rows=m/size;
	int index=0;
	int tag=100;
	

	
	//Distribute number of rows to each process	
	double *nrows = calloc(size, sizeof(int));
	for(size_t i=0; i<size; i++){
		nrows[i]=rows;
	}
	//Distribute restrows
	int rest=m%size;
	for(size_t i =0; i<rest; i++){
		nrows[size-i]++;
	}
	
	// Grid points
	real *grid = mk_1D_array(n+1, false);
	for (size_t i = 0; i < n+1; i++) {
		grid[i] = i * h;
	}
	
	//Make a displacement vector to keep track for each rank
	size_t *displacement=calloc(size+1, sizeof(size_t));
	displacement[0]=0;	
	for(size_t i = 1; i<size; i++){
		displacement[i]=displacement[i-1]+nrows[i-1];
	}

	 // The diagonal of the eigenvalue matrix of T
	real *diag = mk_1D_array(m, false);
	for (size_t i = 0; i < m; i++) {
		diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
	}

	// Initialize the right hand side data
	real **b = mk_2D_array(nrows[rank], m, false);
	real **bt = mk_2D_array(nrows[rank], m, false);
	real *z = mk_1D_array(nn, false);               // Hva er denne til?
	for (size_t i = 0; i < nrows[rank]; i++) {
		for (size_t j = 0; j < m; j++) {
			b[i][j]=h*h*rhs(grid[displacement[rank]+i], grid[j]);	   		 
			//b[i][j] = h * h * rhs(grid[i], grid[j]);
		}
	}

	// Calculate Btilde^T = S^-1 * (S * B)^T       Ikke ferdig
	for (size_t i = 0; i < nrows[rank]; i++) {
		fst_(b[i], &n, z, &nn);
	}
	//transpose(bt, b, m);
	for (size_t i = 0; i < nrows[rank]; i++) {               // Ikke ferdig
		fstinv_(bt[i], &n, z, &nn);
	}

	// Solve Lambda * Xtilde = Btilde
	for (size_t i = 0; i < nrows[rank]; i++) {
		for (size_t j = 0; j < m; j++) {
			bt[i][j] = bt[i][j] / (diag[i] + diag[displacement[rank]+j]); // Hvor er algoritmen for dette
		}
	}

	// Calculate X = S^-1 * (S * Xtilde^T)
	for (size_t i = 0; i < nrows[rank]; i++) {
		fst_(bt[i], &n, z, &nn);
	}
	
	

//	transpose(b, bt, m);

	for (size_t i = 0; i < nrows[rank]; i++) {
		fstinv_(b[i], &n, z, &nn);
	}

	// Calculate maximal value of solution
	double u_max = 0.0;
	for (size_t i = 0; i < nrows[rank]; i++) {
		for (size_t j = 0; j < m; j++) {
			u_max = u_max > b[i][j] ? u_max : b[i][j];
		}
	}

	// All to find reduce to find maximum global sum
	double global_u_max=0.0;
	MPI_Reduce(&u_max, &global_u_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	if(rank==0){	
		printf("global_u_max = %e\n", global_u_max);
	}
	MPI_Finalize();
	return 0;
}

real rhs(real x, real y) {
    return 2 * (y - y*y + x - x*x);
}


/*
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

*/

/*
void transpose(real **bt, real **b, size_t m)
{
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = b[j][i];
        }
    }
}
*/

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
