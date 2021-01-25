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
    // p1 -- p3 -> c2

    unsigned int allele_p1[ 2 ] ; //A1A2
    unsigned int allele_p2[ 2 ] ; //A3A4
    unsigned int allele_p3[ 2 ] ; //A5A6

    double AF ;
    // E3
    // 33% -> 0.075
    // 43% -> 0.15
    // 12% -> 0.25
    // 7%  -> 0.35
    // 5%  -> 0.45

    // L1
    // 0.46	0.39	0.09	0.03	0.03
    // 46% -> 0.075
    // 39% -> 0.15
    // 9%  -> 0.25
    // 3%  -> 0.35
    // 3%  -> 0.45
    double SFS[ 4 ] ;
    SFS[0] = 0.675 ;
    SFS[1] = 0.325 ;
    SFS[2] = 0.0 ;
    SFS[3] = 0.0 ;
    double freq_threshold[ 5 ] ;
    freq_threshold[0] = 0.125 ;
    freq_threshold[1] = 0.25 ;
    freq_threshold[2] = 0.37.5 ;
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
      fp = fopen ( "hsp.txt", "w" ) ;
    } else {
      fp = fopen ( "hsp.txt", "a" ) ;
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


/*
 int save = 1 ;
 int print = 0 ;
 int TM = 20 ;

 int HAP_ID = 0 ; // all hapltype ID appeared so far

 struct haplotype{
 int ID ;
 int allele ;
 int num_in_pop ;
 int num_child ;
 double freq_in_pop ;
 double next_freq ;
 } ;

 struct haplotype *hap, newHap ;

 int main( void ) {
 time_t start_time = time( NULL ) ;
 double progress_time ;
 gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937 ) ;
 gsl_rng_set( r, start_time * getpid() ) ;

 // Basic params
 int N = 100 ; // corresponding to 2N in the case of diploid
 double u = 0.05 ; // mutation rate per locus

 int T = 5 * N ;
 int printT = 100 ;
 double lambda = 5.0 ;
 double phi = 10.1 ;
 double NB_para = phi / (lambda + phi) ;

 // Variables
 int iM, iH, iN ;
 int H0 ;
 double x ;
 int M ;

 // Functions
 int sel_hap() ;
 int check_newhap() ;
 int modify_hap() ;
 //int modify_polysites() ;
 //void freq_spec() ;
 //void mutation_record() ;

 // Initialize
 int t = 0 ;
 int nhap = 1 ;
 int num_child_tmp ;
 int total_child ;
 double *Prob ;
 unsigned int *Pop ;
 Prob = (double *) calloc ( nhap, sizeof(double) ) ;
 Pop = (unsigned int *) calloc ( nhap, sizeof(unsigned int) ) ;
 hap = (struct haplotype *) malloc ( nhap*sizeof(struct haplotype) );
 hap[0].ID = 0 ;
 hap[0].num_in_pop = N ;
 hap[0].num_child = 0 ;
 hap[0].freq_in_pop = 1.0 ;
 hap[0].next_freq = 0.0 ;
 //hap[0].allele = (int *) malloc ( sizeof(int) ) ;
 hap[0].allele = HAP_ID ;
 newHap.allele = 0 ;

 // Main iteration
 while ( t < T ) {
 // 1. Multiple mutations
 M = gsl_ran_binomial ( r, u, N ) ;
 //printf( "M = %d\n", M ) ;

 for ( iM = 0; iM < M; iM++ ) {

 // Select haplotype (H0)
 x = gsl_rng_uniform ( r ) ;
 H0 = sel_hap( x, nhap ) ;

 // Make new haplotype
 HAP_ID ++ ;
 if ( HAP_ID == 1000 ) HAP_ID = 0 ;
 //newHap.allele = hap[H0].allele ; // Copy H0 allele
 //printf( "new_allele = %d\n", newHap.allele ) ;

 // Check & create new haplotype
 nhap = check_newhap( nhap, H0, N ) ;
 } // for(iM = 0; iM < M; iM++)
 if ( (print == 1) && (t % printT == 0) ) {
 printf( "BEFORE sampling: nHap = %d, t = %d\n", nhap, t ) ;
 for ( iH = 0; iH < nhap; iH++ ) {
 printf( "%d", hap[iH].allele ) ;
 printf( "\tID:%d\tFreq:%.3f\tn:%d", hap[iH].ID, hap[iH].freq_in_pop, hap[iH].num_in_pop ) ;
 printf( "\n" ) ;
 }
 printf( "\n" ) ;
 }
 //break ;
 // 2. Reproduction
 total_child = 0 ;
 for ( iH = 0; iH < nhap; iH++ ) {
 num_child_tmp = 0 ;
 for ( iN = 0; iN < hap[iH].num_in_pop; iN++ ) {
 num_child_tmp = num_child_tmp + gsl_ran_negative_binomial ( r, NB_para, phi ) ;
 }
 hap[iH].num_child = num_child_tmp ;
 total_child = total_child + num_child_tmp ;
 }
 if ( (print == 1) && (t % printT == 0) ) {
 printf( "reproduction: nHap = %d, t = %d\n", nhap, t ) ;
 for ( iH = 0; iH < nhap; iH++ ) {
 printf( "%d", hap[iH].allele ) ;
 printf( "\tID:%d\tchild:%d", hap[iH].ID, hap[iH].num_child ) ;
 printf( "\n" ) ;
 }
 printf( "\n" ) ;
 }

 // Save
 if ( (save == 1) && (TM > 0) && (t == ( T - TM - 1) ) ) {
 FILE *fp_hap ;
 fp_hap = fopen( "hap_number_TM.txt", "w" ) ;
 for ( iH = 0; iH < nhap; iH++ ) {
 fprintf( fp_hap, "%d\t%d\n", hap[iH].allele, hap[iH].num_child ) ;
 }
 fclose( fp_hap );
 }
 if ( (save == 1) && (t == (T-1)) ) {
 FILE *fp_hap ;
 fp_hap = fopen( "hap_number.txt", "w" ) ;
 for ( iH = 0; iH < nhap; iH++ ) {
 fprintf( fp_hap, "%d\t%d\n", hap[iH].allele, hap[iH].num_child ) ;
 }
 fclose( fp_hap );
 }

 // Allocation
 Prob = (double *)realloc( Prob, nhap*sizeof(double) ) ;
 Pop = (unsigned int *)realloc( Pop, nhap*sizeof(unsigned int) ) ;

 for ( iH = 0; iH < nhap; iH++ ) {
 Prob[iH] = (double)(hap[iH].num_child) / total_child ;
 //printf( "%.3f\n", Prob[iH] ) ;
 Pop[iH] = (unsigned int)(hap[iH].num_child) ;
 }

 // 3. Binomial sampling
 gsl_ran_multinomial(r, nhap, N, Prob, Pop) ;
 for ( iH = 0; iH < nhap; iH++ ) {
 hap[iH].num_in_pop = (int)Pop[iH] ;
 hap[iH].freq_in_pop = (double)(hap[iH].num_in_pop)/(double)N ;
 }
 if ( (print == 1) && (t % printT == 0) ) {
 printf( "AFTER sampling: nHap = %d, t = %d\n", nhap, t ) ;
 for ( iH = 0; iH < nhap; iH++ ) {
 printf( "%d", hap[iH].allele ) ;
 printf( "\tID:%d\tFreq:%.3f\tn:%d", hap[iH].ID, hap[iH].freq_in_pop, hap[iH].num_in_pop ) ;
 printf( "\n" ) ;
 }
 printf( "\n" ) ;
 }
 // 4. Modify haplotype
 nhap = modify_hap( nhap ) ;

 if ( (print == 1) && (t % printT == 0) ) {
 printf( "AFTER modifying: nHap = %d, t = %d\n", nhap, t ) ;
 for ( iH = 0; iH < nhap; iH++ ) {
 printf( "%d", hap[iH].allele ) ;
 printf( "\tID:%d\tFreq:%.3f\tn:%d", hap[iH].ID, hap[iH].freq_in_pop, hap[iH].num_in_pop ) ;
 printf( "\n" ) ;
 }
 printf( "\n" ) ;
 }

 t++ ;

 } // while t


 progress_time = difftime ( time(NULL), start_time ) ;
 printf ( "progress time [%.0f]sec\n", progress_time ) ;
 return ( 0 );
 }
 //-------------------------------------------------------------------------------------
 int sel_hap( double x, int nhap ) {
 double cumFreq ;
 int iH ;
 int H0 = 0 ;
 cumFreq = 0.0;
 for ( iH = 0; iH < nhap; iH++ ) {
 cumFreq += (hap+iH) -> freq_in_pop ;
 if ( cumFreq > x ) {
 H0 = iH ; // H0 is a haptype that will occurs mutation
 break ;
 }
 }
 return ( H0 ) ;
 }

 int check_newhap ( int nhap, int H0, int N, int new_allele ) {
 if ( hap[H0].num_in_pop == 1) {
 //extinction of a hap being mutated -> no change!
 } else {
 // New haplotype
 hap[H0].num_in_pop -- ;
 hap[H0].freq_in_pop = (double) hap[H0].num_in_pop/N ;
 nhap++ ;
 hap = (struct haplotype *) realloc( hap, nhap*sizeof(struct haplotype) ) ;
 hap[nhap-1].ID = nhap-1 ;
 hap[nhap-1].num_in_pop = 1 ;
 hap[nhap-1].freq_in_pop = (double) hap[nhap-1].num_in_pop/N ;
 //hap[nhap-1].allele = (int *) malloc( sizeof(int) ) ;
 hap[nhap-1].allele = HAP_ID ;
 }
 return ( nhap ) ;
 }

 int modify_hap ( int nhap ) {
 int iH ;
 int check_zero1, check_zero2 ;

 check_zero1 = 0 ;
 while ( check_zero1 == 0 ) { // Roop unless hap with pop_num=0 (hap0) disappear
 while ( hap[nhap-1].num_in_pop == 0 ){ // Remove last hap as long as it is hap0
 //free( hap[nhap-1].allele ) ;
 nhap-- ;
 hap = (struct haplotype *) realloc( hap, nhap*sizeof(struct haplotype) ) ;
 }
 for ( iH = 0; iH < nhap; iH++ ) { // Substitution (hap0 -> nhap-1)
 if ( hap[iH].num_in_pop == 0 ) {
 hap[iH].num_in_pop = hap[nhap-1].num_in_pop ;
 hap[iH].freq_in_pop = hap[nhap-1].freq_in_pop ;
 hap[iH].allele = hap[nhap-1].allele ;
 //free( hap[nhap-1].allele ) ;
 nhap-- ;
 hap = (struct haplotype *) realloc( hap, nhap*sizeof(struct haplotype) );
 break;
 }
 }

 check_zero2 = 0 ;
 for ( iH = 0; iH < nhap; iH++ ) { // Check hap0
 if ( hap[iH].num_in_pop == 0 ) {
 check_zero2 = 1 ;
 }
 }
 if ( check_zero2 == 0 ) check_zero1 = 1 ;
 }
 return ( nhap ) ;
 }
 */
