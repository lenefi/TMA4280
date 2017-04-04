
#include "mpi.h"

#include <cstdio>
#include <cstring>
#include <iomanip>
#include <sstream>

typedef double real;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
	int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int m = 2*size;
  int n = m / size; 
  int r = m % size; 
  int np = n + (rank < r);
  int * nrows =  new int[size];
  int * bsize =  new int[size];
  int * displ =  new int[size];
	for (int i=0; i < size; ++i)
  {
		nrows[i] = n + (i<r);
		bsize[i] = np * nrows[i];
		//Print block info
    //printf("%2u: nrows = %u; bsize = %u\n",i, nrows[i], bsize[i]);
  }
  displ[0] = 0;
	for (int i=1; i < size; ++i)
  {
		displ[i] = displ[i-1]+bsize[i-1];
  }

  /* Initialize matrices */
	real * Avalues = new real[np*m];
	real ** A = new real*[np];
	real * Atvalues = new real[np*m];
	real ** At = new real*[np];
	real Aij = 1.0; 
	for (int i=0; i < rank; ++i)
  {
		Aij+= nrows[i]*m;
	}
	for (int i=0; i < np; ++i)
  {
		A[i] = Avalues + i * m;
		At[i] = Atvalues + i * m;
		for (int j=0; j < m; ++j, ++Aij)
    {
			A[i][j] = Aij;
    }
  }

  /* Copy sub-blocks with packing */
  real * Apck = Atvalues;
  for (int p = 0, off_rp = 0; p < size; ++p, off_rp+=nrows[p])
	{
    //Print row offset
    //printf("offrp = %u\n", off_rp);
	  for (int i = 0; i < np; ++i, Apck+=nrows[p])
	  {
			memcpy(Apck, A[i] + off_rp, nrows[p]*sizeof(real)); 
	  }
	}

  /* Print sub-blocks with packing */
	std::stringstream st;
  Apck = Atvalues;
  for (int p = 0, off_bp = 0; p < size; ++p, off_bp += bsize[p])
  {
		for (int ij = 0; ij < bsize[p]; ++ij)
 		{ 
			st << std::setw(8) << Apck[off_bp + ij];
		}
		st << "\n";
  }
	printf("S%u\n%s\n", rank, st.str().c_str());

  /* Exchange blocks */
  MPI_Alltoallv(Atvalues, bsize, displ, MPI_DOUBLE, Avalues, bsize, displ, MPI_DOUBLE, MPI_COMM_WORLD);

  /* Transpose blocks */
  Apck = Avalues;
  for (int p = 0, off_rp = 0; p < size; ++p, off_rp+=nrows[p], Apck +=bsize[p])
	{
	  for (int i = 0; i < np; ++i)
	  {
			for (int j = 0; j < nrows[p]; ++j)
      {
				At[i][off_rp + j] = Apck[j*nrows[p]+i];
      }
	  }
	}

	MPI_Barrier(MPI_COMM_WORLD);
  /* Display matrix */
  int ij = 0;
	std::stringstream ss;
	for (int i=0; i < np; ++i)
  {
		for (int j=0; j < m; ++j, ++ij)
    {
			ss << std::setw(8) << Atvalues[ij];
    }
		ss << "\n";
  }
	printf("A^t%u\n%s", rank, ss.str().c_str());

	delete [] bsize;
	delete [] nrows;
	delete [] A;
	delete [] Avalues;
	MPI_Finalize();	
	return 0;
}
