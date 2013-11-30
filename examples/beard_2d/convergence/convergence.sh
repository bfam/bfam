#!/bin/bash

if [ $# != 4 ]
then
  echo "convergence.sh BEARD NUM_PROCS DG_N LP"
  exit
fi


error_token="error:"
energy_token="d_energy:"
err1=0
err2=0

for num in {1..10}; do
  cat convergence.lua |                                                         \
    sed "\$amin_level=0\nmax_level=min_level+$4\nN=$3\nstatic_refinement=$num-1"\
    > tmp_$2_$3_$4.lua

  OUT=$(mpirun -n $2 $1 tmp_$2_$3_$4.lua | grep error:)

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
