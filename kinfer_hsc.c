/*--------------------------------------------------------
 Mar.3, 2020
 Simulate genotype data of HSC for "Relate"

 NOTE:

 [MISC]
 210201 Including SFS
 203010 Start coding

 [Operation]
 gcc -DHAVE_INLINE -lgsl -lm -lgslcblas  -Wall -o kinfer_hsc kinfer_hsc.c
./kinfer_hsc [maker_number] [error_rate]
 --------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <sys/stat.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <unistd.h> // for getpid()

int usage () ;
int debug = 0 ;

int main( int argc, char *argv[] ) {
  if ( argc != 3) {
    fprintf(stderr,"Command line arguments are incorrect\n");
    usage();
  }
  time_t start_time = time ( NULL ) ;
  double progress_time ;
  gsl_rng *r = gsl_rng_alloc ( gsl_rng_mt19937 ) ;
  gsl_rng_set ( r, start_time * getpid()  ) ;

  const int maker_number =      atoi(argv[1]) ; //1
  const double error_rate =     atof(argv[2]) ; //2

  printf ( "maker_number = %d, error_rate = %.5f\n", maker_number, error_rate) ;
  unsigned int G1 = 0 ;
  unsigned int G2 = 0 ;

  for ( int i = 0; i < maker_number; i++ ) {

    // p1 -- p2 -> c1
    // p1 -- p3 -> c2
    // c2 -- p4 -> c3
    // c1 -- p5 -> c4
    // compare c3 vs c4

    unsigned int allele_p1[ 2 ] ; //A1A2
    unsigned int allele_p2[ 2 ] ; //A3A4
    unsigned int allele_p3[ 2 ] ; //A5A6
    unsigned int allele_p4[ 2 ] ; //A7A8
    unsigned int allele_p5[ 2 ] ; //A9A10

    double AF ;
    double SFS[ 4 ] ;
    SFS[0] = 0.675 ;
    SFS[1] = 0.325 ;
    SFS[2] = 0.0 ;
    SFS[3] = 0.0 ;
    double freq_threshold[ 4 ] ;
    freq_threshold[0] = 0.125 ;
    freq_threshold[1] = 0.25 ;
    freq_threshold[2] = 0.375 ;
    freq_threshold[3] = 0.5 ;

    double ran_AF = gsl_rng_uniform ( r ) ;
    double sum_tmp = 0.0 ;
    for ( int j = 3; j >= 0; j-- ) {
      sum_tmp += SFS[ j ] ;
      if ( ran_AF < sum_tmp ) {
        AF = freq_threshold[ j ] ;
        break ;
      }
    }

    if (debug) printf ("AF = %.3f\n", AF) ;

    if (debug) printf ("== p1 == \n") ;
    for ( int a = 0; a < 2; a++ ) {
      if ( gsl_rng_uniform ( r ) < AF ) { // allele type is 1 or 2
        allele_p1[ a ] = 1 ;
      } else {
        allele_p1[ a ] = 0 ;
      }
      if (debug) printf ( "%d\n", allele_p1[ a ] ) ;
    }
    if (debug) printf ("== p2 == \n") ;
    for ( int a = 0; a < 2; a++ ) {
      if ( gsl_rng_uniform ( r ) < AF ) { // allele type is 1 or 2
        allele_p2[ a ] = 1 ;
      } else {
        allele_p2[ a ] = 0 ;
      }
      if (debug) printf ( "%d\n", allele_p2[ a ] ) ;
    }
    if (debug) printf ("== p3 == \n") ;
    for ( int a = 0; a < 2; a++ ) {
      if ( gsl_rng_uniform ( r ) < AF ) { // allele type is 1 or 2
        allele_p3[ a ] = 1 ;
      } else {
        allele_p3[ a ] = 0 ;
      }
      if (debug) printf ( "%d\n", allele_p3[ a ] ) ;
    }
    if (debug) printf ("== p4 == \n") ;
    for ( int a = 0; a < 2; a++ ) {
      if ( gsl_rng_uniform ( r ) < AF ) { // allele type is 1 or 2
        allele_p4[ a ] = 1 ;
      } else {
        allele_p4[ a ] = 0 ;
      }
      if (debug) printf ( "%d\n", allele_p4[ a ] ) ;
    }
    if (debug) printf ("== p5 == \n") ;
    for ( int a = 0; a < 2; a++ ) {
      if ( gsl_rng_uniform ( r ) < AF ) { // allele type is 1 or 2
        allele_p5[ a ] = 1 ;
      } else {
        allele_p5[ a ] = 0 ;
      }
      if (debug) printf ( "%d\n", allele_p4[ a ] ) ;
    }

    unsigned int allele_c1[ 2 ] ;
    double rand = gsl_rng_uniform ( r ) ;
    if ( rand < 0.25 ) {
      allele_c1[ 0 ] = allele_p1[ 0 ] ;
      allele_c1[ 1 ] = allele_p2[ 0 ] ;
    } else if ( rand < 0.5 ) {
      allele_c1[ 0 ] = allele_p1[ 0 ] ;
      allele_c1[ 1 ] = allele_p2[ 1 ] ;
    } else if ( rand < 0.75 ) {
      allele_c1[ 0 ] = allele_p1[ 1 ] ;
      allele_c1[ 1 ] = allele_p2[ 0 ] ;
    } else {
      allele_c1[ 0 ] = allele_p1[ 1 ] ;
      allele_c1[ 1 ] = allele_p2[ 1 ] ;
    }

    for ( int a = 0; a < 2; a++ ) {
      if ( gsl_rng_uniform ( r ) < (error_rate) ) {
        if ( allele_c1[ a ] == 0 )  {
          allele_c1[ a ] = 1 ;
        } else {
          allele_c1[ a ] = 0 ;
        }
      }
    }

    unsigned int allele_c2[ 2 ] ;
    rand = gsl_rng_uniform ( r ) ;
    if ( rand < 0.25 ) {
      allele_c2[ 0 ] = allele_p1[ 0 ] ;
      allele_c2[ 1 ] = allele_p3[ 0 ] ;
    } else if ( rand < 0.5 ) {
      allele_c2[ 0 ] = allele_p1[ 0 ] ;
      allele_c2[ 1 ] = allele_p3[ 1 ] ;
    } else if ( rand < 0.75 ) {
      allele_c2[ 0 ] = allele_p1[ 1 ] ;
      allele_c2[ 1 ] = allele_p3[ 0 ] ;
    } else {
      allele_c2[ 0 ] = allele_p1[ 1 ] ;
      allele_c2[ 1 ] = allele_p3[ 1 ] ;
    }

    for ( int a = 0; a < 2; a++ ) {
      if ( gsl_rng_uniform ( r ) < (error_rate) ) {
        if ( allele_c2[ a ] == 0 )  {
          allele_c2[ a ] = 1 ;
        } else {
          allele_c2[ a ] = 0 ;
        }
      }
    }

    unsigned int allele_c3[ 2 ] ;
    rand = gsl_rng_uniform ( r ) ;
    if ( rand < 0.25 ) {
      allele_c3[ 0 ] = allele_c2[ 0 ] ;
      allele_c3[ 1 ] = allele_p4[ 0 ] ;
    } else if ( rand < 0.5 ) {
      allele_c3[ 0 ] = allele_c2[ 0 ] ;
      allele_c3[ 1 ] = allele_p4[ 1 ] ;
    } else if ( rand < 0.75 ) {
      allele_c3[ 0 ] = allele_c2[ 1 ] ;
      allele_c3[ 1 ] = allele_p4[ 0 ] ;
    } else {
      allele_c3[ 0 ] = allele_c2[ 1 ] ;
      allele_c3[ 1 ] = allele_p4[ 1 ] ;
    }

    for ( int a = 0; a < 2; a++ ) {
      if ( gsl_rng_uniform ( r ) < (error_rate) ) {
        if ( allele_c3[ a ] == 0 )  {
          allele_c3[ a ] = 1 ;
        } else {
          allele_c3[ a ] = 0 ;
        }
      }
    }

    unsigned int allele_c4[ 2 ] ;
    rand = gsl_rng_uniform ( r ) ;
    if ( rand < 0.25 ) {
      allele_c4[ 0 ] = allele_c1[ 0 ] ;
      allele_c4[ 1 ] = allele_p5[ 0 ] ;
    } else if ( rand < 0.5 ) {
      allele_c4[ 0 ] = allele_c1[ 0 ] ;
      allele_c4[ 1 ] = allele_p5[ 1 ] ;
    } else if ( rand < 0.75 ) {
      allele_c4[ 0 ] = allele_c1[ 1 ] ;
      allele_c4[ 1 ] = allele_p5[ 0 ] ;
    } else {
      allele_c4[ 0 ] = allele_c1[ 1 ] ;
      allele_c4[ 1 ] = allele_p5[ 1 ] ;
    }

    for ( int a = 0; a < 2; a++ ) {
      if ( gsl_rng_uniform ( r ) < (error_rate) ) {
        if ( allele_c4[ a ] == 0 )  {
          allele_c4[ a ] = 1 ;
        } else {
          allele_c4[ a ] = 0 ;
        }
      }
    }


    if (debug) printf("c1: %d%d\n", allele_c1[ 0 ],allele_c1[ 1 ]) ;
    if (debug) printf("c2: %d%d\n", allele_c2[ 0 ],allele_c2[ 1 ]) ;
    if (debug) printf("c3: %d%d\n", allele_c3[ 0 ],allele_c3[ 1 ]) ;
    if (debug) printf("c4: %d%d\n", allele_c3[ 0 ],allele_c4[ 1 ]) ;

    G1 = 0 ;
    G2 = 0 ;
    // 00 -> 0
    // 01 -> 1
    // 11 -> 2

    if ( (allele_c4[ 0 ] == 0)&&(allele_c4[ 1 ] == 0) ) {
      G1 = 0 ;
    } else if ( (allele_c4[ 0 ] == 1)&&(allele_c4[ 1 ] == 1) ) {
      G1 = 2 ;
    } else {
      G1 = 1 ;
    }

    if ( (allele_c3[ 0 ] == 0)&&(allele_c3[ 1 ] == 0) ) {
      G2 = 0 ;
    } else if ( (allele_c3[ 0 ] == 1)&&(allele_c3[ 1 ] == 1) ) {
      G2 = 2 ;
    } else {
      G2 = 1 ;
    }

    if (debug) printf("%d\t%d\n",G1,G2) ;
    //--- save ---
    FILE *fp ;
    if ( i == 0 ) {
      fp = fopen ( "hsc.txt", "w" ) ;
    } else {
      fp = fopen ( "hsc.txt", "a" ) ;
    }
    fprintf ( fp, "%d\t%d\n", G1, G2) ;
    fclose ( fp ) ;
  }

  progress_time = difftime ( time(NULL), start_time ) ;
  printf ( "progress time [%.0f]sec\n", progress_time ) ;
  return ( 0 ) ;
}
//------------------------------------------------------------------------------------------
int usage ( ) {
  fprintf ( stderr, "usage: kinfer_hsc maker_number error_rate\n" );
  exit ( 1 ) ;
}
