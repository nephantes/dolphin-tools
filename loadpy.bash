#!/bin/bash
a=$( module list 2>&1)
echo $a

if [[ $a != *"python/2.7.5"* ]];
then
   module load python/2.7.5
   b=($( python --version))
   echo $b   
fi
  
