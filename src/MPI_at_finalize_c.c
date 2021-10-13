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

#define MAXFNS 32

typedef void (*fnptr)(void) ;
static fnptr fn_list[MAXFNS] ;
static int used_fns = 0 ;

#define MAXWINS 2048

static MPI_Win winlist[MAXWINS] ;      // windows list
static int init = 1 ;                  // initialize flag
static int used = 0 ;                  // number of registered windows in list

void MPI_Free_Tracked_Windows() ;
void MPI_Track_Window_f(MPI_Fint win_f, int insert) ;
void MPI_Track_Window(MPI_Win win, int insert) ;

// register a function to be executed before MPI_Finalize
int MPI_At_Finalize( void (*fn)(void) ){
  if(used_fns < MAXFNS){
    fn_list[used_fns] = fn ;           // store in function list if not full
  }else{
    return 1 ;                         // table is full
  }
  used_fns++ ;
  return 0 ;
}

// execute functions registered with MPI_At_Finalize (in reverse order of registration)
void MPI_At_Finalize_exec(){
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
  MPI_At_Finalize(MPI_Free_Tracked_Windows) ;  // register function to auto free one sided windows not freed yet
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

void MPI_Track_Window(  // insert/remove C MPI window into/from tracked list
  MPI_Win win,              // MPI C one sided window
  int insert)               // 1 = insert, 0 = remove
{
  int i ;
  if(init){                 // initialize list entries to null windows
    for(i = 0 ; i < MAXWINS ; i++) winlist[i] = MPI_WIN_NULL ;
    init = 0 ;
  }
  if(DEBUG) printf("MPI_Track_Window, insert= %d, used = %d \n", insert, used);
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

void MPI_Track_Window_f(MPI_Fint win_f, int insert){  // insert/remove Fortran MPI window into/from tracked list
  MPI_Win win = MPI_Win_f2c(win_f) ;                      // translate Fortran MPI window into C window
  MPI_Track_Window(win, insert) ;
  if(DEBUG) printf("MPI_Track_Window_f used = %d, insert = %d \n", used, insert);
}

// free all registered windows not freed yet
void MPI_Free_Tracked_Windows(){
  int i, count ;
  count = 0 ;
  if(DEBUG) printf("MPI_Free_Tracked_Windows : start\n");
  if(used == 0) return ;
  for(i = 0 ; i < used ; i++){
    if(winlist[i] != MPI_WIN_NULL){                       // entry still valid
      if(DEBUG) printf("MPI_Track_Window_f, freeing %d\n",i);
      PMPI_Win_free( &(winlist[i]) ) ;
      winlist[i] = MPI_WIN_NULL ;
      count++ ;
    }
  }
  used = 0 ;
  printf("MPI_Free_Tracked_Windows : freed %d windows\n", count) ;
}
