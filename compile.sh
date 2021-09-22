#!/bin/zsh

gcc $(gsl-config --cflags) $(gsl-config --libs) -Wall -o kinfer_unr kinfer_unr.c 
gcc $(gsl-config --cflags) $(gsl-config --libs) -Wall -o kinfer_po kinfer_po.c 
gcc $(gsl-config --cflags) $(gsl-config --libs) -Wall -o kinfer_hsp kinfer_hsp.c 
gcc $(gsl-config --cflags) $(gsl-config --libs) -Wall -o kinfer_fsp kinfer_fsp.c 
gcc $(gsl-config --cflags) $(gsl-config --libs) -Wall -o kinfer_hun kinfer_hun.c
gcc $(gsl-config --cflags) $(gsl-config --libs) -Wall -o kinfer_hsc kinfer_hsc.c 