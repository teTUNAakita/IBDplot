/*--------------------------------------------------------
 Mar.3, 2020
 Simulate genotype data of HSP for "Relate"

 NOTE:

 [MISC]
 200428 Including SFS
 200303 Start coding

 [Operation]
 gcc -DHAVE_INLINE -lgsl -lm -lgslcblas  -Wall -o kinfer_hsp kinfer_hsp.c
./kinfer_hsp [maker_number] [error_rate]
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
  //const int ind_number =        atoi(argv[1]) ; //3

  printf ( "maker_number = %d, error_rate = %.5f\n", maker_number, error_rate) ;
  unsigned int G1 = 0 ;
  unsigned int G2 = 0 ;

  for ( int i = 0; i < maker_number; i++ ) {

    // p1 -- p2 -> c1
    // p1 -- p2 -> c2
    // compare c1 vs c2

    unsigned int allele_p1[ 2 ] ; //A1A2
    unsigned int allele_p2[ 2 ] ; //A3A4

    double AF ;
    double SFS[ 8 ] ;
    SFS[0] = 0.036 ;
    SFS[1] = 0.100 ;
    SFS[2] = 0.064 ;
    SFS[3] = 0.118 ;
    SFS[4] = 0.109 ;
    SFS[5] = 0.200 ;
    SFS[6] = 0.100 ;
    SFS[7] = 0.273 ;
    double freq_threshold[ 8 ] ;
    freq_threshold[0] = 0.125 ;
    freq_threshold[1] = 0.175 ;
    freq_threshold[2] = 0.225 ;
    freq_threshold[3] = 0.275 ;
    freq_threshold[4] = 0.325 ;
    freq_threshold[5] = 0.375 ;
    freq_threshold[6] = 0.425 ;
    freq_threshold[7] = 0.475 ;
    double ran_AF = gsl_rng_uniform ( r ) ;
    double sum_tmp = 0.0 ;
    for ( int j = 7; j >= 0; j-- ) {
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
      allele_c2[ 1 ] = allele_p2[ 0 ] ;
    } else if ( rand < 0.5 ) {
      allele_c2[ 0 ] = allele_p1[ 0 ] ;
      allele_c2[ 1 ] = allele_p2[ 1 ] ;
    } else if ( rand < 0.75 ) {
      allele_c2[ 0 ] = allele_p1[ 1 ] ;
      allele_c2[ 1 ] = allele_p2[ 0 ] ;
    } else {
      allele_c2[ 0 ] = allele_p1[ 1 ] ;
      allele_c2[ 1 ] = allele_p2[ 1 ] ;
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
    if (debug) printf("c1: %d%d\n", allele_c1[ 0 ],allele_c1[ 1 ]) ;
    if (debug) printf("c2: %d%d\n", allele_c2[ 0 ],allele_c2[ 1 ]) ;

    G1 = 0 ;
    G2 = 0 ;
    // 00 -> 0
    // 01 -> 1
    // 11 -> 2

    if ( (allele_c1[ 0 ] == 0)&&(allele_c1[ 1 ] == 0) ) {
      G1 = 0 ;
    } else if ( (allele_c1[ 0 ] == 1)&&(allele_c1[ 1 ] == 1) ) {
      G1 = 2 ;
    } else {
      G1 = 1 ;
    }

    if ( (allele_c2[ 0 ] == 0)&&(allele_c2[ 1 ] == 0) ) {
      G2 = 0 ;
    } else if ( (allele_c2[ 0 ] == 1)&&(allele_c2[ 1 ] == 1) ) {
      G2 = 2 ;
    } else {
      G2 = 1 ;
    }

    if (debug) printf("%d\t%d\n",G1,G2) ;
    //--- save ---
    FILE *fp ;
    if ( i == 0 ) {
      fp = fopen ( "fsp.txt", "w" ) ;
    } else {
      fp = fopen ( "fsp.txt", "a" ) ;
    }
    fprintf ( fp, "%d\t%d\n", G1, G2) ;
    fclose ( fp ) ;

    FILE *fp1 ;
    if ( i == 0 ) {
      fp1 = fopen ( "AF.txt", "w" ) ;
    } else {
      fp1 = fopen ( "AF.txt", "a" ) ;
    }
    fprintf ( fp1, "%.3f\n", AF) ;
    fclose ( fp1 ) ;
  }

  progress_time = difftime ( time(NULL), start_time ) ;
  printf ( "progress time [%.0f]sec\n", progress_time ) ;
  return ( 0 ) ;
}
//------------------------------------------------------------------------------------------
int usage ( ) {
  fprintf ( stderr, "usage: kinfer_hsp maker_number error_rate\n" );
  exit ( 1 ) ;
}
