//     functions for C and FORTRAN programming
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

#define MAXWINS 2048
#define MAXFNS 32

static MPI_Win winlist[MAXWINS] ;      // windows list
static int init = 1 ;                  // initialize flag
static int used = 0 ;                  // number of registered windows in list

typedef void (*fnptr)(void) ;
static fnptr fn_list[MAXFNS] ;
static int used_fns = 0 ;

void RPN_MPI_free_tracked_windows() ;

int MPI_At_Finalize(void (*fn)(void)){   // register a function to be executed before MPI_Finalize
  if(used_fns < MAXFNS){
    fn_list[used_fns] = fn ;             // store in function list if not full
  }
  used_fns++ ;
  return 0 ;
}

void MPI_At_Finalize_exec(){             // execute registered functions (in reverse order of registration)
  while(used_fns > 0){
    used_fns-- ;
    fn_list[used_fns]() ;
  }
}

static void dummy1(){             // demo at finalize function
  printf("dummy1 function\n");
}

static void dummy2(){             // demo at finalize function
  printf("dummy2 function\n");
}

static void dummy3(){             // demo at finalize function
  printf("dummy3 function\n");
}

int MPI_Init(int *argc, char ***argv){             // intercept MPI_Init call from user
  if(DEBUG) printf("MPI_Init C\n");
  MPI_At_Finalize(RPN_MPI_free_tracked_windows) ;
  if(DEBUG) MPI_At_Finalize(dummy1) ;
  if(DEBUG) MPI_At_Finalize(dummy2) ;
  if(DEBUG) MPI_At_Finalize(dummy3) ;
  return PMPI_Init(argc, argv) ;
}

int MPI_Finalize(){                                // intercept MPI_Finalize call from user
  if(DEBUG) printf("MPI_finalize C\n");
  MPI_At_Finalize_exec() ;                         // execute registered functions
  return PMPI_Finalize() ;
}

void RPN_MPI_track_window(  // insert/remove C MPI window into/from tracked list
  MPI_Win win,              // MPI C one sided window
  int insert)               // 1 = insert, 0 = remove
{
  int i ;
  if(init){                 // initialize list entries to null windows
    for(i = 0 ; i < MAXWINS ; i++) winlist[i] = MPI_WIN_NULL ;
    init = 0 ;
  }
  if(DEBUG) printf("RPN_MPI_track_window, insert= %d, used = %d \n", insert, used);
  if(insert){                                          // insert

    for(i = 0 ; i < used ; i++){
      if(winlist[i] == win){                           // already in list
      if(DEBUG) printf("entry %d is a duplicate\n",i);
        return ;
      }
    }
    for(i = 0 ; i < used ; i++){
      if(winlist[i] == MPI_WIN_NULL){                  // resuse a free entry in list
        winlist[i] = win ;
      if(DEBUG) printf("reusing entry %d\n", i);
        return ;
      }
    }
    if(used >= MAXWINS) return ;                       // table full, OUCH
    winlist[used] = win ;                              // add a new entry
    if(DEBUG) printf("adding entry %d\n", used);
    used++ ;
    return ;

  }else{                                               // remove

    for(i = 0 ; i < used ; i++){
      if(winlist[i] == win){                           // entry found in table, replace with null window
        winlist[i] = MPI_WIN_NULL ;
        if(DEBUG) printf("removing entry %d\n", i);
        return ;
      }
    }
  }

}

void RPN_MPI_track_window_f(MPI_Fint win_f, int insert){  // insert/remove Fortran MPI window into/from tracked list
  MPI_Win win = MPI_Win_f2c(win_f) ;                      // translate Fortran MPI window into C window
  RPN_MPI_track_window(win, insert) ;
  if(DEBUG) printf("RPN_MPI_track_window_f used = %d, insert = %d \n", used, insert);
}

void RPN_MPI_free_tracked_windows(){                      // free all registered windows not freed yet
  int i, count ;
  count = 0 ;
  if(DEBUG) printf("RPN_MPI_free_tracked_windows : start\n");
  if(used == 0) return ;
  for(i = 0 ; i < used ; i++){
    if(winlist[i] != MPI_WIN_NULL){                       // entry still valid
      if(DEBUG) printf("RPN_MPI_track_window_f, freeing %d\n",i);
      PMPI_Win_free( &(winlist[i]) ) ;
      winlist[i] = MPI_WIN_NULL ;
      count++ ;
    }
  }
  used = 0 ;
  printf("RPN_MPI_free_tracked_windows : freed %d windows\n", count) ;
}

// intercept user call to some MPI functions

int MPI_Win_free(MPI_Win *win){
printf("MPI_Win_free C %p\n", win) ;
  RPN_MPI_track_window(*win, 0) ;     // untrack window
  return PMPI_Win_free(win) ;
}

int MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win *win){
  int ierr = PMPI_Win_create(base, size, disp_unit, info, comm, win) ;
  if(DEBUG) printf("MPI_Win_create C %p\n", win) ;
  RPN_MPI_track_window(*win, 1) ;     // track new window
  return ierr ;
}

int MPI_Win_allocate (MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win){
  int ierr = PMPI_Win_allocate (size, disp_unit, info, comm, baseptr, win) ;
  if(DEBUG) printf("MPI_Win_allocate C %p\n", win) ;
  RPN_MPI_track_window(*win, 1) ;     // track new window
  return ierr ;
}

int MPI_Win_allocate_shared (MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win){
  int ierr = PMPI_Win_allocate_shared (size, disp_unit, info, comm, baseptr, win) ;
  if(DEBUG) printf("MPI_Win_allocate_shared C %p\n", win) ;
  RPN_MPI_track_window(*win, 1) ;     // track new window
  return ierr ;
}

int MPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win){
  int ierr = PMPI_Win_create_dynamic(info, comm, win) ;
  if(DEBUG) printf("MPI_Win_create_dynamic C %p\n", win) ;
  RPN_MPI_track_window(*win, 1) ;     // track new window
  return ierr ;
}

