#!/bin/bash
typeset -a array=('' '(1)' '(1,1)' '(1,1,1)' '(1,1,1,1)')
typeset -a arrayp=('' '(:)' '(:,:)' '(:,:,:)' '(:,:,:,:)')
typeset -a arrayy=('' '(lbound(what,1))' '(lbound(what,1),lbound(what,2))' '(lbound(what,1),lbound(what,2),lbound(what,3))' '(lbound(what,1),lbound(what,2),lbound(what,3),lbound(what,4))')
type=integer
kind=4
cat <<EOT
!!end interface     !InTf!
!! interface rpn_mpi_ptr     !InTf!
#define IN_RPN_MPI_ptr
EOT
for typex in integer real
do
for kind in 4 8
do

for dim in 1 2 3 4
do
  [[ $typex == i* ]] && type='i'
  [[ $typex == r* ]] && type='r'
  cat <<EOT
  subroutine RPN_MPI_ptr_${type%??????}${kind}_${dim}d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpnmpi_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_MPI_types.inc'
  $typex(KIND=$kind), dimension${array[$dim]},intent(IN), target :: what    !InTf!
  type(rpnmpi_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what${array[$dim]})
  return
  end subroutine RPN_MPI_ptr_${type%??????}${kind}_${dim}d    !InTf!

EOT
done


done
done
cat <<EOT
!!end interface rpn_mpi_ptr      !InTf!
!!interface     !InTf!
EOT
