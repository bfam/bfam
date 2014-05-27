#!/bin/bash

if [ $# != 7 ]
then
  echo "convergence.sh BEARD LUA_FILE NUM_PROCS DG_N1 DG_N2 DG_N3 LP"
  exit
fi


error_token="error:"
energy_token="d_energy:"
err1=0
err2=0

SED=sed
type gsed >/dev/null 2>&1 && SED=gsed


for num in {1..6}; do
  cat $2 |                                                         \
    $SED "\$amin_level=0\nmax_level=min_level+$7\nN1=$4\nN2=$5\nN3=$6\nstatic_refinement=$num-1\nlsrk_method = \"KC54\"\ndt_fudge = 0.5"\
    > $2_$3_$4_$5_$6_$7_lsrk_tmp.lua

  OUT=$(mpirun -n $3 $1 $2_$3_$4_$5_$6_$7_lsrk_tmp.lua | grep error:)

  RES="0"

  for token in $OUT
  do
    if [ "$RES" == "1" ]; then
      err1=$(echo ${token} | sed 's/e/\*10\^/' | sed 's/+//')
    elif [ "$RES" == "2" ]; then
      energy=$(echo ${token} | sed 's/e/\*10\^/' | sed 's/+//')
    fi

    if [ "$token" == "$error_token" ]; then
      RES="1"
    elif [ "$token" == "$energy_token" ]; then
      RES="2"
    else
      RES="0"
    fi
  done


  if [[ $err2 > 0 ]]; then
    rate=$(echo "l(($err1)/($err2))/l(0.5)" | bc -l)

    echo "level: $num denergy: $energy error: $err1 rate: $rate"
  else
    echo "level: $num denergy: $energy error: $err1"
  fi

  err2=$err1
done
