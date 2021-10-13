//     functions for C and FORTRAN MPI
//     Copyright (C) 2021  Recherche en Prevision Numerique
// 
//     This software is free software; you can redistribute it and/or
//     modify it under the terms of the GNU Lesser General Public
//     License as published by the Free Software Foundation,
//     version 2.1 of the License.
//
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//     Lesser General Public License for more details.
//
// Author : 
//     M.Valin, Recherche en Prevision Numerique, Oct 2021
//
#include <stdio.h>
#include <mpi.h>

#if ! defined(DEBUG)
#define DEBUG 0
#endif

void MPI_Track_Window(MPI_Win win, int insert) ;

// intercept user call to some MPI one sided windows functions

int MPI_Win_free(MPI_Win *win){
printf("MPI_Win_free C %p\n", win) ;
  MPI_Track_Window(*win, 0) ;     // untrack window
  return PMPI_Win_free(win) ;
}

int MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win *win){
  int ierr = PMPI_Win_create(base, size, disp_unit, info, comm, win) ;
  if(DEBUG) printf("MPI_Win_create C %p\n", win) ;
  MPI_Track_Window(*win, 1) ;     // track new window
  return ierr ;
}

int MPI_Win_allocate (MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win){
  int ierr = PMPI_Win_allocate (size, disp_unit, info, comm, baseptr, win) ;
  if(DEBUG) printf("MPI_Win_allocate C %p\n", win) ;
  MPI_Track_Window(*win, 1) ;     // track new window
  return ierr ;
}

int MPI_Win_allocate_shared (MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win){
  int ierr = PMPI_Win_allocate_shared (size, disp_unit, info, comm, baseptr, win) ;
  if(DEBUG) printf("MPI_Win_allocate_shared C %p\n", win) ;
  MPI_Track_Window(*win, 1) ;     // track new window
  return ierr ;
}

int MPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win){
  int ierr = PMPI_Win_create_dynamic(info, comm, win) ;
  if(DEBUG) printf("MPI_Win_create_dynamic C %p\n", win) ;
  MPI_Track_Window(*win, 1) ;     // track new window
  return ierr ;
}
