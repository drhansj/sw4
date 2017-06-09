#include <mpi.h>
#include <string>
#include <iostream>
#include "EW.h"

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sched.h>
#include <omp.h>

using namespace std;

/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c */
static char *cpuset_to_cstr(cpu_set_t *mask, char *str)
{
  char *ptr = str;
  int i, j, entry_made = 0;
  for (i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, mask)) {
      int run = 0;
      entry_made = 1;
      for (j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, mask)) run++;
        else break;
      }
      if (!run)
        sprintf(ptr, "%d,", i);
      else if (run == 1) {
        sprintf(ptr, "%d,%d,", i, i + 1);
        i++;
      } else {
        sprintf(ptr, "%d-%d,", i, i + run);
        i += run;
      }
      while (*ptr != 0) ptr++;
    }
  }
  ptr -= entry_made;
  *ptr = 0;
  return(str);
}

int main( int argc, char** argv )
{
   //MPI_Init(&argc, &argv);
   int myRank, thread;
   double  time_start, time_end;
   cpu_set_t coremask;
   char clbuf[7 * CPU_SETSIZE], hnbuf[64];

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

/*Rahul - adding stuff to print the thread affinity */
   memset(clbuf, 0, sizeof(clbuf));
   memset(hnbuf, 0, sizeof(hnbuf));
   (void)gethostname(hnbuf, sizeof(hnbuf));
   if(myRank == 0)
       printf("Printing thread affinity in the form of MPI-Task\t thread-in-MPI_Task\t core-id mapped to thread\n");
#pragma omp parallel private(thread, coremask, clbuf)
   {
     thread = omp_get_thread_num();
     (void)sched_getaffinity(0, sizeof(coremask), &coremask);
     cpuset_to_cstr(&coremask, clbuf);
     #pragma omp barrier
     printf("MPI-Task %d, thread %d, on %s. (core id = %s)\n",
 	   myRank, thread, hnbuf, clbuf);
   }

   MPI_Barrier(MPI_COMM_WORLD);
   time_start = MPI_Wtime();

   string filename;
   if( argc <= 1 )
   {
      if( myRank == 0 )
      {
	 cout  << "ERROR: ****No input file specified!" << endl;
	 for (int i = 0; i < argc; ++i)
	    cout << "Argv[" << i << "] = " << argv[i] << endl;
      }
      //MPI_Finalize();
      //return 1;
   }
   else
   {
      filename = argv[1];
      EW simulation(filename);
      //MPI_Finalize();
      //return 0;
   }

   MPI_Barrier(MPI_COMM_WORLD);
   time_end = MPI_Wtime();
   if(myRank == 0) cout <<  " Total running time: " << time_end - time_start << endl;

   MPI_Finalize();

   return 0;
}
