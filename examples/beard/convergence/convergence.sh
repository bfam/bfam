#!/bin/bash

if [ $# != 3 ]
then
  echo "convergence.sh BEARD NUM_PROCS DG_N"
  exit
fi


error_token="error:"
err1=0
err2=0

for num in {1..5}; do
  cat convergence.lua | sed "\$amin_level=$num\nmax_level=$num\nN = $3" > tmp.lua

  OUT=$(mpirun -n $2 $1 tmp.lua | grep error:)

  RES="0"

  for token in $OUT
  do
    if [ "$RES" == "1" ]; then
      err1=$(echo ${token} | sed 's/e/\*10\^/' | sed 's/+//')
    fi

    if [ "$token" == "$error_token" ]; then
      RES="1"
    fi
  done


  if [[ $err2 > 0 ]]; then
    rate=$(echo "l(($err1)/($err2))/l(0.5)" | bc -l)

    echo "level: $num error: $err1 rate: $rate"
  else
    echo "level: $num error: $err1"
  fi

  err2=$err1
done
