#!/bin/bash
if [[ -z $1 ]] ; then
  echo "interface"
  grep '!InTf!' | sed -e 's/^!!//' -e 's/^ [ ]*/      /' -e 's/!InTf!//' -e 's/ [ ]*!.*//' -e 's/^[ ]*(/         (/'
  echo "end interface"
else
  [[ -f $1 ]] || exit 0
  if grep -q '!InTf!' $1 ; then
    echo "#if ! defined(IN_${1%.*})"
#    echo "interface"
    grep '!InTf!' $1 | sed -e 's/^!!//' -e 's/^ [ ]*/      /' -e 's/!InTf!//' -e 's/ [ ]*!.*//' -e 's/^[ ]*(/         (/'
#    echo "end interface"
    echo "#endif"
  fi
fi
