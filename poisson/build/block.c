




int main(int n1, int n2){

int rank, size;
MPI_IInit(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&size);

int N=atoi(argv[1]);
int blocked = 0;
if(argc > 2)
	blocked = atoi(argv[2]);
int sizes[2];
MPI_Datatype filetype;
createFileView(&filetype, sizes,N,rank,sizes,blocked);

int periodic[2]; periodic[0]=periodic[1]=0;
MPI_Comm comm;
MPI_Cart_create(MPI_COMM_WORLD,2,sizes,periodic,0,&comm);

int coord[2];
MPI_Cart_coords(comm,rank,2,coords);

int rows = N/sizes[0];
int cols = N/sizes[1];
int bufsize = rows*cols;
int K;
if (argc > 3)
K=atoi(argv[3]);
K=1;

double**A = createMatrix(rows,cols);
setupMatrix(A,rows,cols,coord[0]*rows,coord[1]*cols,N);

MPI_File f;
MPI_File_open(MPI_COMM_WORLD, "combined.mat",MPI_MODE_CREATE|MPI_MODE_RDWR,MPI_INFO_NULL,&f);
MPI_File_set_view(f,0,MPI_COUDLE,filetype,"native",MPI_INFO_NULL);

for(int n=0;n<K;++n)
    MPI_File_write(f,A[0],bufsize,MPI_DOUBLE,MPI_STATUS_IGNORE);

MPI_File_close(&f);

MPI_Finalize();

return 0;
}
