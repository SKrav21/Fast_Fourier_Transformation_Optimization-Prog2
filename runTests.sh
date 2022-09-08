#!/bin/bash
  echo -e "Size\t\tNo-opt\t\tCode-Motion\tMem ref\t\t2x1\t\t2x2"
  for n in {10..50..10}; do
     gcc -O0 -DN=$n -oprog2 prog2.c -lm
     ./prog2
  done
