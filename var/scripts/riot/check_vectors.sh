#!/bin/bash

#----------------------------------------------
#Check for presence of necessary vectors
#----------------------------------------------
cd $CWD

it="$1"; NSAMP="$2"; rand_stage="$3"; innerit="$4"

it0=$it
if [ $it -lt 10 ]; then it0="0"$it0; fi

innerit0=$innerit
if [ $innerit -lt 10 ]; then innerit0="0"$innerit0; fi
if [ $innerit -lt 100 ]; then innerit0="0"$innerit0; fi
if [ $innerit -lt 1000 ]; then innerit0="0"$innerit0; fi

vectors2=("" "" "")

if [ $rand_type -eq 3 ]; then
   if [ $rand_stage -eq 2 ]; then
      vectors=("ghat.e*.p" "dummy")
      nvec=("$NSAMP" "0")
   fi
   if [ $rand_stage -eq 4 ]; then
      vectors=("yhat.e*.p" "qhat.e*.p")
      vectors2=(".iter$innerit0" "" "")
      nvec=("$NSAMP" "$((NSAMP*innerit))")
   fi
fi

if [ $rand_type -eq 6 ]; then
   if [ $RIOT_RESTART -eq 1 ]; then
     vectors=("omega.e*.p" "yhat.e*.p")
     nvec=("$NSAMP" "$NSAMP")
   else
     vectors=("omega.e*.p" "yhat.e*.p" "ghat.p")
     nvec=("$NSAMP" "$NSAMP" "1")
   fi
fi

if [ $rand_type -eq 2 ]; then
   if [ $RIOT_RESTART -eq 1 ]; then
      vectors=("ahat.e*.p")
      nvec=("$NSAMP")
   else
      vectors=("ahat.e*.p" "ghat.p")
      nvec=("$NSAMP" "1")
   fi
fi

vcount=0
for var in ${vectors[@]}
do
   if [ ${nvec[$vcount]} -gt 0 ]; then 
      echo "Working on $var files"
      suffix=${vectors2[$vcount]}
      #Test for the presence of each vector type
      ls ../vectors_$it0/$var"0000"*$suffix
      dummy=`ls ../vectors_$it0/$var"0000"*$suffix | wc -l`

      echo "$dummy present of ${nvec[$vcount]} $var"0000"*$suffix files"
      if [ $dummy -ne ${nvec[$vcount]} ]; then 
         echo "ERROR: Missing or extra $var""0000"
         echo $((rand_stage*100+vcount*10+4)); exit $((rand_stage*100+vcount*10+4));
      fi
   fi
   vcount=$((vcount+1))
done
