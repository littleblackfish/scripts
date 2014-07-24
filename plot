#!/bin/bash


command="gnuplot -p -e \"plot "

for y in $@
do
#  echo $y
  for x in $y*
  do
    command=$command\'$x\'\ w\ l\ ,  
  done

done
command=$command\'asdvbu\'\"

echo $command
#echo $@

eval $command

