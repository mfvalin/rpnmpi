/* RPN_MPI - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 2021  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 */
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <fcntl.h>
#include <mpi.h>

// given a communicator (possibly spanning more than 1 SMP node)
// on each node, get a shared memory area allowing each process
// to have the amount of memory it needs
//
// if a node communicator is available, it may be passed in [INOUT] argument 2 (n_comm)
//
// if no node communicator is available, a communicator variable with the value MPI_COMM_NULL
// may be passed, and upon exit, its value will be replaced with a valid node communicator
//
// the process with rank 0 on the SMP node will try to allocate enough
// shared memory to satisfy the requests of all Processes on the node.
// a pointer to the total memory is returned to all participants
// it is up to the processes 
//
// g_comm   [IN]    : global communicator
// n_comm   [INOUT] : node communicator (may be MPI_COMM_NULL upon entry)
//                    the node communicator is returned if it was MPI_COMM_NULL
// shm_size [IN]    : size in 32 bit words of desired shared memory area for this process
// returns a pointer to the shared memory region (NULL if an error occurs)
void *MPI_Shmget_memory(MPI_Comm g_comm, MPI_Comm *n_comm, unsigned int shm_size)
{
  size_t size ;
  uint64_t size_one, size_all ;
  int id, myrank, error, errors ;
  struct shmid_ds shm_buf ;
  void *ptr = NULL ;

  if(*n_comm == MPI_COMM_NULL){            // node communicator not supplied yet
    MPI_Comm_rank(g_comm,&myrank) ;        // global rank is used as key for global communicator split
    MPI_Comm_split_type(g_comm, MPI_COMM_TYPE_SHARED, myrank, MPI_INFO_NULL, n_comm) ; // split g_comm by node
  }

  size_one = shm_size * sizeof(uint32_t) ;
  MPI_Allreduce(&size_one, &size_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, *n_comm) ; // add all sizes
  size = size_all ;
  
  MPI_Comm_rank(*n_comm, &myrank) ;                /* rank on node */
  if(myrank == 0) {                                /* rank 0 allocates shared memory and broadcasts the shared memory id */
    id = shmget(IPC_PRIVATE, size, IPC_CREAT|S_IRUSR|S_IWUSR) ;
    if(id >= 0) ptr = shmat(id,NULL,0) ;           /* attach memory segment */
    if(id >= 0) shmctl(id,IPC_RMID,&shm_buf) ;     /* mark segment for deletion to make sure it is released when all processes terminate */
    MPI_Bcast(&id,1,MPI_INTEGER,0,*n_comm) ;       /* all processes get segment id */
  }else{
    MPI_Bcast(&id,1,MPI_INTEGER,0,*n_comm) ;       /* all processes get segment id */
    if(id >= 0) ptr = shmat(id,NULL,0) ;           /* attach memory segment */
  }
  if(ptr == (void *) -1) ptr = NULL ;

  error = (ptr == NULL) ? 1 : 0 ;
  MPI_Allreduce(&error, &errors, 1, MPI_INTEGER, MPI_SUM, *n_comm) ; // add all error flags
  if(errors > 0) {                                                   // there were errors
    if(ptr != NULL) shmdt(ptr) ;                                     // detach memory segment if it was attached
    if(myrank == 0) printf("ERROR: could not get shared memory segment\n");
  }

  MPI_Barrier(g_comm);                                               /* all processes hgave attached the segment */

  return ptr;                                        /* return pointer tio shared memory area */
}

// Fortran wrapper for above
//F_StArT
// interface
//   function MPI_Shmget_memory(g_comm, n_comm, shm_size) result(ptr) BIND(C,name='MPI_Shmget_memory_f08')
//     import :: C_PTR, C_INT
//     implicit none
//     integer(C_INT), intent(IN), value :: g_comm
//     integer(C_INT), intent(INOUT)     :: n_comm
//     integer(C_INT), intent(IN), value :: shm_size
//     type(C_PTR) :: ptr
//   end function MPI_Shmget_memory
// end interface
//F_EnD
void *MPI_Shmget_memory_f08(int g_comm, int *n_comm, unsigned int shm_size)
{
  void *ptr ;
  MPI_Fint f_g_comm = g_comm ;
  MPI_Fint f_n_comm = *n_comm ;
  MPI_Comm c_g_comm = MPI_Comm_f2c(f_g_comm) ;   // translate Fortran communicators
  MPI_Comm c_n_comm = MPI_Comm_f2c(f_n_comm) ;

  ptr = MPI_Shmget_memory(c_g_comm, &c_n_comm, shm_size) ;
  if(MPI_Comm_f2c(f_n_comm) == MPI_COMM_NULL) {  // inbound value was MPI_COMM_NULL, replace with proper value
    *n_comm = MPI_Comm_c2f(c_n_comm) ;
  }
  return ptr ;
}
