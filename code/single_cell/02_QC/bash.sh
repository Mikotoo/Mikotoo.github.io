#!/bin/bash
#CSUB -J test
#CSUB -q fat2
#CSUB -o test.out
#CSUB -e test.error
#CSUB -n 64
#CSUB -R span[hosts=1]

#python integrate.py
python dimenRedu.py