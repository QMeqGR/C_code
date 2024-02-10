/*
  version changes:
  5.2.2
  22 Jan 2016
  Bug fix for -o switch. It set all the other occupancies to
  the same value if you only set one occupancy. (when going through
  the atom loop)

  5.2.1
  Fri Mar 21 15:54:41 CDT 2014
  Important bugfix. Previous version could overestimate
  intensity because the peak would fit in overlapping bins.
  Now only the first matching bin is found. This should only be
  a problem in structures with lots of peaks and wide tolerance.

  5.2.0
  Tue Mar 18 16:40:50 CDT 2014
  Allow modification of occupancies of multiple atoms using
  input string such as "12.050:1.300:7.099" meaning
  Z=12 has occupancy  50%
  Z=1  has occupancy 300%
  Z=7  has occupancy  99%

  5.1.0
  Thu Jan  3 11:02:26 CST 2013
  Modify the occupancy of Z atoms using command line switch.
  This will alter the occupancy of ALL atoms of type Z.

  5.0.0
  Fri Mar 16 16:03:35 CDT 2012
  Add scattering cross sections for all atoms for neutron
  use -1 for deuterium
  
  4.5.1
  Thu Mar 25 15:21:37 CDT 2010
  Reduce q-min default value to 0.1 as I was missing peaks
  in important cases.
  
  4.5
  Thu Aug 28 15:19:29 CDT 2008
  Add d-spacing to q output.
  Add option for 'quiet' (no header)

  4.4
  Tue Oct 31 15:10:20 PST 2006
  Added code for biso values for individual atoms

  4.3
  Tue Oct 25 13:50:56 PDT 2005
  Added switch for scaling the lattice vectors.

*/

/**************************************************************/
/*                                                            */
/*  powder.c                                                  */
/*                                                            */
/*  version 0.5                                               */
/*  Generates the q values and the intensity for diffraction  */
/*  from  cubic or tetragonal lattices.                       */
/*                                                            */
/*  General procedure:                                        */
/*  I calculate the q values first, since I know the          */
/*  d-spacing formula for the lattices, then run through      */
/*  the h,k,l and bin the intensities.                        */
/*                                                            */
/*  Note: must use x--a   y--b   z--c                         */
/*        in the data files                                   */
/*                                                            */
/**************************************************************/
/*                                                            */
/*  file format for xray or neutron (lookup table for neutron */
/*  scattering factors is included in this source file)       */
/*                                                            */
/*  n                                                         */
/*  x1 y1 z1 Z                                                */
/*  x2 y2 z2 Z                                                */
/*                                                            */
/*  xn yn zn Z                                                */
/*                                                            */
/**************************************************************/
/*                                                            */
/*  version 1.0                                               */
/*  Eric Majzoub                                              */
/*  Washington University                                     */
/*  27 April 2000                                             */
/*                                                            */
/*  Modified on 10 October 2001:                              */
/*  Added capability for doing tetragonal lattices            */
/*  version 1.1                                               */
/*                                                            */
/*  Modified on 11 January 2002:                              */
/*  version 1.2                                               */
/*  Added Debye-Waller factor for reduced intensity at large  */
/*  q-values, instead of the 1/q^2 I was using before.        */
/*  Cleaned up the code a bit and made some options available */
/*  on the command line instead of compile-time options.      */
/*                                                            */
/*  Modified on 29 April 2002:                                */
/*  version 1.3                                               */
/*  Added b-axis option for full othorhombic lattices.        */
/*                                                            */
/**************************************************************/
/*                                                            */
/*  Modified on 12 May 2003:                                  */
/*  version 1.5                                               */
/*  New option to print out hkl for each peak.                */
/*                                                            */
/*  This works by taking the first hkl for the new 'bin'.     */
/*  The code counts downward from NMAX so the first 'bin' may */
/*  have a large number for the hkl...                        */
/*                                                            */
/**************************************************************/
/*                                                            */
/*  Modified on 14 May 2003:                                  */
/*  version 2.0                                               */
/*  New option for diffraction conditions/geometry            */
/*                                                            */
/*  -Bragg-Brentano geometry factor                           */
/*  -Lorentz polarization factor                              */
/*                                                            */
/*                                                            */
/**************************************************************/
/*                                                            */
/*  Modified on 27 May 2003:                                  */
/*  version 3.0                                               */
/*                                                            */
/*  Added full triclinic hkl generation.                      */
/*                                                            */
/**************************************************************/
/*                                                            */
/*  Version 4.0                                               */
/*  04 June 2003                                              */
/*  Modified on 04 June 2003:                                 */
/*  add switch for primitive lattices, will read in           */
/*  primitive lattice vectors from file                       */
/*                                                            */
/*  - a check on the correctness of the output... used the    */
/*    primitive cell of diamond.  Calculation gives           */
/*    identical output to conventional cell with scaling      */
/*    factor to correct for differences in the number of      */
/*    atoms in the different unit cells.                      */
/*                                                            */
/*                                                            */
/**************************************************************/
/*                                                            */
/*  Version 4.1                                               */
/*  10 Jan 2005                                               */
/*  Eric Majzoub                                              */
/*  Sandia National Laboratories                              */
/*                                                            */
/*  - fixed printing of nan values when tth exceeds 180 deg   */
/*  - removed printing of hkl when prim lattice is used       */
/*    Note: I could probably fix this so it prints the        */
/*          correct values.  As it stands the hkl values for  */
/*          a prim fcc lattice give hkl's with mixed indices  */
/*          and the like.  This is not technically wrong      */
/*          since the lattice vectors are not conventional    */
/*          cubic, but it may confuse me if I forget about    */
/*          this.  I don't want to quote wrong values, so for */
/*          now, I will just refuse to print them with a      */
/*          primitive lattice.                                */
/*                                                            */
/**************************************************************/
/*                                                            */
/*                                                            */
/* Version 4.2                                                */
/* 13 September 2005                                          */
/* Eric Majzoub                                               */
/*                                                            */
/* -changed printing of output to be %.6e, so that very small */
/*  peaks will be seen against the background of larger ones. */
/*  This is only important if somebody sets the -p flag to 0. */
/*                                                            */
/**************************************************************/

/*
 *  gcc ~/src/powder.c ~/src/ericlib.c -lm -Wall -O3 -o ~/bin/powder
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define VERSION "5.2.2"
#define VERSION_DATE "Fri 22 Jan 2016"
#define NMAX 60    /* max h,k,l index */
#define MAXTYP 100 /* max number of types of atoms */
#define PI 3.14159265358979
#define Q_VALS_MAX 5000 /* allow only this many peaks */
#define TOLERANCE 0.005
#define EIGHTPISQUARED 78.9568352042
#define R2D 57.2957795131

/* Global Variables */
int XRAY=0,NEUT=0;
int debug=0;
int zocc_count=0;
int zocc_Z[MAXTYP]={0};

double zocc[MAXTYP]={0.0};
double zocc_pct[MAXTYP]={1.0};

extern char *optarg;
extern int optind,opterr,optopt;


FILE *input=NULL;

/*******************/
/* Data structures */
/*******************/
struct peak{
  double q;
  double intensity;
  int h,k,l;
  int multiplicity;
};

struct bisoval{
  int Z;
  double Biso;
};

/* Function prototypes */
int main(int argc, char *argv[]);
int getopt(int argc, char * const argv[], const char *optstring);
int equal(int,int,int);
int tol(double,double,double);
int cmp(double*,double*);

double q(double,double,double,double,double,double,int,int,int);
double xfact(double,double);
double nfact(double);
double StrFactor(double *,int,int,int,int,double,int);


/**********************************/
/*                                */
/*          MAIN                  */
/*                                */
/**********************************/
int main(int argc, char *argv[])
{
  char *latfile="data.dat";

  char occsep[] = ":";
  char *occtoken;
  char dot[] = ".";
  char *dottoken;
  char tmpocc[100],oswitch[100];

  int i,j,h,k,l;
  int file_info=0;
  int dot_count=0;
  int prim_lat=0;
  int quiet=0;
  int mult;
  int optnum;
  int errval,numatoms;
  int q_n,peak_flag;
  int TWOTHETA=0;
  int LIST_HKL=0;
  int D_COR=0;
  int have_remaining_bins = -1;
  int bins_used=0;
  int num_data=0;
  int ndat_atoms_only=0,ndat_lat_vecs_only=0,ndat_lat_parms_only=0;
  int ndat_biso_only=0,ndat_lat_parms_biso=0,ndat_lat_vecs_biso=0;
  int ntypat=0,match_flag=0;
  extern int XRAY,NEUT,debug;

  double QMAX=8,QMIN=0.1,lat_a=1,lat_c=-1,lat_b=-1;
  double alpha=90, beta=90, gamma=90;
  double qval=0;
  double intensity;
  double *atom;
  double biso=0;
  double tmpz,tZ,killatom=0,kill_all_except=0;
  double tth,ttr,lambda=-1;
  double ptol=0.005;
  double qtol=TOLERANCE;
  double max_intensity=0;
  double output_scale=1;
  double lat_scale=1.0;
  double z_shift=0.0;
  double v[9]={0},R[3]={0},dot01=0,dot02=0,dot12=0,a01=0,a02=0,a12=0;

  struct peak pk[Q_VALS_MAX];
  struct bisoval *bisolist;

  if ( argc == 1 ) {
      printf("\nTry 'powder -h' for help\n\n");
      exit(0);
  }

  while ( (optnum = getopt(argc, argv,
			   "a:A:b:Bc:Cd:Ef:FhK:lLNo:Pp:Q:q:s:S:t:W:Xz:"
			   )) != -1 ){
    switch(optnum) {
    case 'h':
      printf("##################################################\n"
	     "\tpowder, version %s\n"
	     "\tEric Majzoub, %s\n"
	     "\tSandia National Laboratories\n"
             "##################################################\n\n",
	     VERSION,VERSION_DATE);
      printf("Source file: " __FILE__ "\n");
      printf("Compile date: " __DATE__ "\n");
      printf("Compile time: " __TIME__ "\n");
      printf("\nCommand line options (* indicates required):\n"
	     "  -h --- help\n"
	     "\n"
	     "           **********************\n"
	     "                Input options\n"
	     "           **********************\n"
	     "  -f -*- input filename\n"
	     "  -F --- list file formats\n"
	     "\n"
	     "  -L --- read a,b,c and alpha,beta,gamma from file\n"
	     "         (default off) usually input from command line\n"
	     "  -P --- read prim lat input (default conventional cell)\n"
	     "  -B --- read Biso values for each atom (default 1 ang^2)\n"
	     "\n"
	     "  -Q --- qmax (default 8)\n"
	     "  -q --- qmin (default %f)\n"
	     "  -W --- q bin size tolerance (default 0.005)\n"
	     "\n"
	     "  -X --- simulate XRAY powder diffraction\n"
	     "  -N --- simulate NEUTRON powder diffraction\n"
	     "  -C --- use diffraction corrections (caution!) (default off)\n"
	     "         Bragg-Brentano and Lorentz Polarization correction\n"
	     "\n"
	     "  -a -*- a lattice constant [angstroms] (default 1)\n"
	     "  -b -*- b lattice constant (default a)\n"
	     "  -c -*- c lattice constant (default a)\n"
	     "\n"
	     "  -s --- scale all lattice vectors by factor (default 1.0)\n"
	     "  -K -*- make atom Z a non-scattering atom (default 0)\n"
	     "  -A -*- make all atoms EXCEPT Z non-scattering (default 0)\n"
	     "  -o -*- modify occupany of atom type Z (default off)\n"
	     "         use: -o \"Z1.nnn:Z2.nnn\"\n"
	     "         example: -o \"12.099:3.400\" means Z=12 has 99 per cent occupancy,\n"
	     "         and Z=3 has 400 per cent occupancy.\n"
	     "\n"
	     "           **********************\n"
	     "               Output options\n"
	     "           **********************\n"
	     "  -E --- quiet output (no header)\n"
	     "  -t -*- output in two-theta (must enter a lambda (Cu=1.5405))\n"
	     "  -z -*- z-shift (position correction in [mm]) (default 0.0)\n"
	     "         (only works in two theta output mode)\n"
	     "  -p -*- peak intensity tolerance (default 0.005)\n"
	     "         as percentage of maximum peak intensity\n"
	     "  -S -*- output peak intensity scale factor (default 1)\n"
	     "  -l --- generate (hkl, multiplicity) list with peaks\n"
	     "  -d -*- debug level (integer) (default 0..4=everything)\n"
	     "\n\n",QMIN);
      exit(0);
      
    case ':':
      printf("option needs a value\n");
      exit(0);

    case 'E':
      quiet = 1;
      break;
      
    case 'f':
      latfile  = optarg;
      break;

    case 'F':
      printf("\n"
	     "  There are 3 file formats, and an option to add Biso\n"
	     "  values for each element (*) to the end of the file.\n"
	     "  Notes:\n"
	     "  a) all file formats must have atoms in\n"
	     "     reduced (fractional) coordinates.\n"
	     "  b) if alpha, beta, or gamma are not 90 degrees\n"
	     "     they must be input from the file.\n"
	     "  c) default Biso value is 1 for all atoms if they are\n"
	     "     not given in the input file.\n"
	     "\n\n"
	     "1) lattice parameters given on command line:\n"
	     "   (no special command line options)\n\n"
	     "   number of atoms\n"
	     "   x1 y1 z1 Z1\n"
	     "   x2 y2 z2 Z2\n"
	     "   ...\n"
	     "   xn yn zn Zn\n\n"
	     "2) lattice parameters given in file:\n"
	     "   (use -L switch)\n\n"
	     "   number of atoms\n"
	     "   lat_a lat_b lat_c [angstroms]\n"
	     "   alpha beta gamma [degrees]\n"
	     "   x1 y1 z1 Z1\n"
	     "   x2 y2 z2 Z2\n"
	     "   ...\n"
	     "   xn yn zn Zn\n\n"
	     "3) primitive lattice vectors given in file:\n"
	     "   (use -P switch)\n\n"
	     "   number of atoms\n"
	     "   a0 a1 a2 [angstroms]\n"
	     "   b0 b1 b2 [angstroms]\n"
	     "   c0 c1 c2 [angstroms]\n"
	     "   x1 y1 z1 Z1\n"
	     "   x2 y2 z2 Z2\n"
	     "   ...\n"
	     "   xn yn zn Zn\n\n"
	     "*) adding Biso values for each element:\n"
	     "   (add biso values AFTER atoms and add -B switch)\n\n"
	     "   ...\n"
	     "   xn yn zn Zn\n"
	     "   Z1 B1\n"
	     "   Z2 B2\n"
	     "   ...\n"
	     "   Zm Bm\n\n"
	     );
      exit(0);
      break;

    case 'd':
      debug  = atoi(optarg);
      break;

    case 'o':

      //printf("optarg= %s\n",optarg);
      //printf("strlen= %d\n",strlen(optarg));
      /* have to copy this b/c strtok changes the string */
      strcpy(tmpocc,optarg);
      strcpy(oswitch,optarg);
      //printf("tmpocc= %s\n",tmpocc);

      occtoken = strtok( optarg , occsep );
      while ( occtoken != NULL ){
	//printf("token= %s\n",occtoken); fflush(stdout);
	zocc_count++;
	tmpz = atof(occtoken);
	zocc_Z[zocc_count] = (int)floor( tmpz );
	zocc_pct[zocc_count] = ( tmpz - zocc_Z[zocc_count] ) * 10;
	occtoken = strtok( NULL , occsep );
      }

      /* do a check to see if the input is OK */
      /* should be Z.nnn:Z.nnn:Z.nnn */
      /* counts . and : and checks for consistency */
      //printf("tmpocc= %s\n",tmpocc);
      dottoken = strtok( tmpocc , dot );
      while ( dottoken != NULL ){
	//printf("dottoken= %s\n",dottoken); fflush(stdout);
	dot_count++;
	dottoken = strtok( NULL , dot );
      }
      
      if ( dot_count != zocc_count+1 ){
	printf("** Error: Check the command line -o switch settings. Found:\n");
	printf("** dot_count= %d  zocc_count=%d\n", dot_count,zocc_count);
	printf("** Exiting\n"); exit(0);
      }
      // printf("zocc-zocc_Z= %f\n",zocc-zocc_Z);
      // printf("zocc_Z= %d, zocc_pct= %f\n",zocc_Z,zocc_pct);
      break;

    case 'L':
      file_info  = 1;
      break;

    case 'P':
      prim_lat  = 1;
      break;

    case 'C':
      D_COR  = 1;
      break;

    case 'K':
      killatom  = atof(optarg);
      break;

    case 'p':
      ptol  = atof(optarg);
      break;

    case 'W':
      qtol  = atof(optarg);
      break;

    case 'S':
      output_scale  = atof(optarg);
      break;

    case 's':
      lat_scale  = atof(optarg);
      break;

    case 'A':
      kill_all_except  = atof(optarg);
      break;

    case 't':
      TWOTHETA = 1;
      lambda  = atof(optarg);
      break;

    case 'z':
      z_shift  = atof(optarg);
      break;

    case 'X':
      XRAY  = 1;
      break;

    case 'N':
      NEUT  = 1;
      break;

    case 'l':
      LIST_HKL  = 1;
      break;

    case 'B':
      biso  = 1;
      break;

    case 'Q':
      QMAX = atof(optarg);	
      break;
      
    case 'q':
      QMIN = atof(optarg);	
      break;
      
    case 'c':
      lat_c = atof(optarg);	
      break;

    case 'b':
      lat_b = atof(optarg);	
      break;
      
    case 'a':
      lat_a  = atof(optarg);
      break;
      
    case '?':
      printf("unknown option: %c (or option requires a value!)\n",optopt);
      exit(0);
    }
  }


  /* if not tetragonal, assign c length as a */
  if ( lat_b < 0 ) lat_b = lat_a;
  if ( lat_c < 0 ) lat_c = lat_a;

  /* Error checking on input options */
  if ( !XRAY && !NEUT ) {
    printf("Must specify XRAY or NEUTRON diffraction!\n"
	   "  -X --- simulate xray powder diffraction\n"
	   "  -N --- simulate neutron powder diffraction\n");
    exit(0);
  }
  if ( prim_lat == 1 && LIST_HKL == 1 ) {
    printf("\n");
    printf("The H K L listing with primitive lattice vectors is incorrect.\n");
    printf("They won't follow the 'conventional' rules of 'no mixed', etc.\n");
    printf("You must use the conventional cell to get the H K L values.\n");
    printf("Either remove '-l' option, or switch to conventional lattice.\n");
    printf("\n");
    exit(0);
  }

  /* initialize q-vals as -1 */
  for ( q_n=0; q_n<Q_VALS_MAX; q_n++ ) {
    pk[q_n].q = -1;
    pk[q_n].intensity=0;
  }

  /* Open input file and get unit cell information */
  input = fopen(latfile,"r");
  if ( input == NULL )
    {
      printf("Error opening input file\n"); exit(0);
    }

  /* get number of atoms */
  errval = fscanf(input," %d ",&numatoms);
  if ( errval != 1 ) { printf("Bad latfile type.\n"); exit(0); }
  
  /* allocate space for data points */
  /* field : value */
  /* 0 : x */
  /* 1 : y */
  /* 2 : z */
  /* 3 : Z (atomic number as a float)*/
  /* 4 : Biso */
  atom = (double *)malloc( numatoms * 5 * sizeof(double) );
  if ( !atom ) {
      printf("%s\n","Not enough memory for atoms.\n");
      exit(0);
    }
  bisolist = (struct bisoval *)malloc( numatoms * sizeof(struct bisoval) );
  if ( !bisolist ) {
      printf("%s\n","Not enough memory for bisolist.\n");
      exit(0);
    }

  /* initialize bisolist values */
  for(i=0;i<numatoms;i++){
    bisolist[i].Z = -1;
    bisolist[i].Biso = -1;
  }
  
  /* if file_info flag is set, read in alpha, beta, gamma from file */
  if ( file_info == 1 ){
    errval = fscanf(input," %lf %lf %lf ",&lat_a,&lat_b,&lat_c);
    if ( errval != 3 ) {printf("Err reading a, b, c.\n"); exit(0);}
    
    errval = fscanf(input," %lf %lf %lf ",&alpha,&beta,&gamma);
    if ( errval != 3 ) {printf("Err reading alpha, beta, gamma.\n"); exit(0);}
  }
  /* if prim_lat flag is set, read in prim lat vecs from file */
  if ( prim_lat == 1 ){
    for (i=0; i<9; i++) errval = fscanf(input," %lf ",&v[i]);
    if ( errval != 1 ) {printf("Err reading prim lat vecs.\n"); exit(0);}
    
    /* calculate alpha beta and gamma from the prim vecs */
    for (i=0; i<3; i++){
      R[i] = sqrt( v[i*3+0]*v[i*3+0] + v[i*3+1]*v[i*3+1] + v[i*3+2]*v[i*3+2] );
      dot01 += v[0*3+i]*v[1*3+i];
      dot02 += v[0*3+i]*v[2*3+i];
      dot12 += v[1*3+i]*v[2*3+i];
    }
    a01 = atan2(sqrt(R[0]*R[0]*R[1]*R[1]-dot01*dot01),dot01);
    a02 = atan2(sqrt(R[0]*R[0]*R[2]*R[2]-dot02*dot02),dot02);
    a12 = atan2(sqrt(R[1]*R[1]*R[2]*R[2]-dot12*dot12),dot12);
    
    alpha = R2D*a12;
    beta  = R2D*a02;
    gamma = R2D*a01;
    
    lat_a = R[0];
    lat_b = R[1];
    lat_c = R[2];
  }

  /* read in the basis atoms */
  for (i=0; i<numatoms; i++)
    {
      errval = fscanf(input," %lf %lf %lf %lf ",
		      (atom+5*i+0),
		      (atom+5*i+1),
		      (atom+5*i+2),
		      &tZ);
      /* calculate the number of types of atoms */
      for(j=0,match_flag=0;j<numatoms;j++){
	if ( (int)tZ == bisolist[j].Z ) match_flag=1;
      }
      if ( match_flag==0 ) bisolist[ntypat++].Z=tZ;

      /* this code changes the atom scattering factors to zero, etc. */
      if ( tZ == killatom && killatom ) { *(atom+5*i+3)=0; }
      else if ( tZ != kill_all_except && kill_all_except ) { *(atom+5*i+3)=0; }
      else { *(atom+5*i+3) = tZ; }
      if (errval != 4)
	{
	  printf("Bad file type at line %d.\n",i);
	  printf("Check atom count (among other things).\n");
	  free(atom);
	  exit(0);
	}
    }

  /* read in the biso values */
  if ( biso==1 ){
    for(j=0;j<ntypat;j++){
      errval = fscanf(input," %d %lf ",
		      &bisolist[j].Z,
		      &bisolist[j].Biso);
      if (errval != 2 ){
	printf("Error reading Biso values at j=%d\n",j);
	free(atom);
	free(bisolist);
	exit(0);
      }
      if ( debug>0 ){
	printf("# Z = %d   Biso = %f\n",bisolist[j].Z,bisolist[j].Biso);
      }
    }
  }

  for(i=0;i<numatoms;i++){
    if ( (int)*(atom+5*i+3) != 0 ){
      for(j=0;j<ntypat;j++){
	if ( (int)*(atom+5*i+3) == bisolist[j].Z )
	  *(atom+5*i+4)=bisolist[j].Biso;
      }
    }
  }

  /* scale the lattice constants */
  lat_a *= lat_scale;
  lat_b *= lat_scale;
  lat_c *= lat_scale;


  if ( !quiet ) {
    printf("##################################################\n"
	   "#\tpowder, version %s\n"
	   "#\tEric Majzoub\n"
	   "#\tSandia National Laboratories\n"
	   "##################################################\n",VERSION);
    printf("#  Source file: " __FILE__ "\n");
    printf("#  Compile date: " __DATE__ "\n");
    printf("#  Compile time: " __TIME__ "\n");
    printf("#  commandline = ");
    for(i=0; i<argc; i++) printf("%s ",argv[i]); printf("\n");
    printf("#  NOTE: full -o switch arg is: %s\n",oswitch);
    printf("#  input file = %s\n",latfile);
    printf("#  ntypat = %d\n",ntypat);
    printf("#  qmin = %.2f , qmax = %.2f [inv ang]\n",QMIN,QMAX);
    printf("#  Q bin size tolerance is %8.4f\n",qtol);
    if ( file_info == 0  && prim_lat == 0)
      printf("#  Reading lattice parameters from command line.\n");
    if ( file_info == 1 )
      printf("#  Reading lattice parameters from input file.\n");
    if ( prim_lat == 1 )
      printf("#  Primitive lattice: reading vectors from file.\n");
    if ( biso == 1 )
      printf("#  Reading Biso for atoms from file.\n");
    if ( biso == 0 )
      printf("#  Using Biso = 1 for all atoms.\n");
    if ( prim_lat == 0 )
      printf("#  lattice  constants = %12.6f %12.6f %12.6f [angstroms]\n",
	     lat_a, lat_b, lat_c);
    if ( prim_lat == 1 )
      printf("#  prim. vec. lengths = %12.6f %12.6f %12.6f [angstroms]\n",
	     lat_a, lat_b, lat_c);
    printf("#  alpha, beta, gamma = %12.6f %12.6f %12.6f [degrees]\n",
	   alpha, beta, gamma);
    printf("#  lat_scale = %10.6f\n",lat_scale);
    if ( zocc_count > 0 ) {
      printf("#  %8s%12s  (zocc_count= %d)\n","Z","occup",zocc_count);
      for ( i=1; i < zocc_count+1 ; i++ ){
	printf("#  %8d%12.3f\n",zocc_Z[i],zocc_pct[i]);
      }

    }
    if (XRAY) printf("#  Simulating XRAY diffraction\n");
    if (NEUT) printf("#  Simulating NEUTRON diffraction\n");
    printf("#  Peak tolerance cutoff is %5.4f of maximum\n",ptol);
    if (D_COR==0) printf("#  RAW structure factors (no corrections)\n");
    if (D_COR==1) printf("#  Using B-B and LP corrections\n");
    if (killatom) printf("#  Atom Z=%.0f is non-scattering\n",killatom);
    if (TWOTHETA) printf("#  lambda = %5.4f [angstroms]\n",lambda);
    if (debug) printf("#  Debug level is %d\n",debug);
    if (z_shift != 0 && TWOTHETA) printf("#  z-shift = %4.3f [mm]\n",z_shift);
    if (TWOTHETA) printf("#  Output is (tth,intensity)\n");
    if (!TWOTHETA) printf("#  Output is ( q , d , intensity )\n");
    if (LIST_HKL) printf("%1s%19s%7s%5s%5s%18s\n",
			 "#","Output (integers)","H","K","L","Multiplicity");
    printf("#############################################################\n");
  }
  
  /* Error checking on input file ... */
  /* check to see if a b c, alpha, beta, gamma are set in file, etc... */
  rewind(input);
  for ( i=0; fscanf(input, " %lf ", &tZ)>0 ; i++);
  num_data=i;
  ndat_atoms_only=1+4*numatoms;
  ndat_lat_vecs_only=1+4*numatoms+9;
  ndat_lat_parms_only=1+4*numatoms+6;
  ndat_biso_only=1+4*numatoms+2*ntypat;
  ndat_lat_parms_biso=1+4*numatoms+6+2*ntypat;
  ndat_lat_vecs_biso=1+4*numatoms+9+2*ntypat;

  if ( debug > 1 ) {
    printf("#\n");
    printf("# num_data                = %5d\n", num_data);
    printf("# 1+4*numatoms            = %5d (atoms only)\n",
	   1+4*numatoms);
    printf("# 1+4*numatoms+9          = %5d (lat vecs only)\n",
	   1+4*numatoms+9);
    printf("# 1+4*numatoms+6          = %5d (lat parms only)\n",
	   1+4*numatoms+6);
    printf("# 1+4*numatoms+2*ntypat   = %5d (Biso only)\n",
	   1+4*numatoms+2*ntypat);
    printf("# 1+4*numatoms+6+2*ntypat = %5d (lat parms, and Biso)\n",
	   1+4*numatoms+6+2*ntypat);
    printf("# 1+4*numatoms+9+2*ntypat = %5d (lat vecs, and Biso)\n",
	   1+4*numatoms+9+2*ntypat);
    printf("#\n");
  }

  if ( debug>0 && num_data == ndat_atoms_only ){
    printf("# Looks like file has atoms only.\n"
	   "# You cannot then use the L,P, or B switches.\n");
  }
  if ( debug>0 && num_data == ndat_lat_vecs_only ){
    printf("# Looks like file has atoms and lattice vectors.\n"
	   "# Are you using the -P switch?\n");
  }
  if ( debug>0 && num_data == ndat_lat_parms_only ){
    printf("# Looks like file has atoms and lattice parameters.\n"
	   "# Are you using the -L switch?\n");
  }
  if ( debug>0 && num_data == ndat_biso_only ){
    printf("# Looks like file has atoms and Biso values only.\n"
	   "# Are you using the -B switch?\n");
  }
  if ( debug>0 && num_data == ndat_lat_parms_biso ){
    printf("# Looks like file has atoms, lat parms, and Biso values.\n"
	   "# Are you using the -L and -B switches?\n");
  }
  if ( debug>0 && num_data == ndat_lat_vecs_biso ){
    printf("# Looks like file has atoms, lat vecs, and Biso values.\n"
	   "# Are you using the -P and -B switches?\n");
  }

  /* most basic error checking */
  if ( num_data > numatoms*4+1 && file_info==0 && prim_lat==0 && biso==0 ) {
    printf("Are lattice parameters or primitive vectors given in the file?\n"
	   "(see -L and -P flags on input)\n");
    printf("# num_data = %d\n", num_data);
    printf("# 4*numatoms+1 = %d\n", 4*numatoms+1);
    printf("# difference = %d\n",num_data - (4*numatoms+1));
    exit(0);
  }
  /* error checks for -L flag */
  if ( num_data != (numatoms*4+1+6) && file_info==1 && biso==0 ){
    printf("#\n");
    printf("# WARNING: Check atom count! (numatoms=%d in input file)\n"
	   "# File looks like it has %d atoms.\n"
	   "# Only using first %d atoms.\n",numatoms,(num_data-7)/4,numatoms);
  }
  if ( num_data != (numatoms*4+1+9) && prim_lat==1 && biso==0 ){
    printf("#\n");
    printf("# WARNING: Check atom count! (numatoms=%d in input file)\n"
	   "# File looks like it has %d atoms.\n"
	   "# Only using first %d atoms.\n",numatoms,(num_data-10)/4,numatoms);
  }
  if ( num_data != (numatoms*4+1) && file_info==0 && prim_lat==0 && biso==0 ){
    printf("#\n");
    printf("# WARNING: Check atom count! (numatoms=%d in input file)\n"
	   "# File looks like it has %d atoms.\n"
	   "# Only using first %d atoms.\n",numatoms,(num_data-1)/4,numatoms);
  }


  if ( debug > 1 ) {
    printf("# Angles:\n");
    printf("# radians    01 02 12     %.6f  %.6f  %.6f\n",a01,a02,a12);
    printf("# alpha beta gamma (rads) %.6f  %.6f  %.6f\n",a12,a02,a01);
    printf("# alpha beta gamma (degs) %.3f  %.3f  %.3f\n",
	   R2D*a12,R2D*a02,R2D*a01);
  }
  if ( debug > 0 ) {
    printf("#\n");
    printf("# Number of atoms in input file = %d\n",numatoms);
    printf("# Atom input list :\n");
    printf("#%8s%12s%12s%12s%8s%15s\n","atom","x","y","z","Z","Biso");
    for (i=0; i<numatoms; i++)
      {
	printf("#%8d%12.6f%12.6f%12.6f%8.0f%15.5f\n",
	       i,
	       *(atom+5*i+0),
	       *(atom+5*i+1),
	       *(atom+5*i+2),
	       *(atom+5*i+3),
	       *(atom+5*i+4));
      }
    printf("#\n");
  }
  fclose(input);




  /*
   *
   *  The algorithm:  it runs over h,k,l and if the qvalue is
   *  within the range allowed, calculates the intensity and
   *  puts it in the correct qvalue bin.
   *
   *
   */

  if ( debug > 1 ) printf("# Beginning hkl loops... \n");
  for (h=NMAX; h>-(NMAX+1); h--) {
    for (k=NMAX; k>-(NMAX+1); k--) {
      for (l=NMAX; l>-(NMAX+1); l--) {

	if ( h==0 && k==0 && l==0 ) continue;
	
	qval = q(lat_a,lat_b,lat_c,alpha,beta,gamma,h,k,l);
	
	if ( qval > QMAX || qval < QMIN ) continue;

	intensity = StrFactor(atom,numatoms,h,k,l,qval,biso);

	/* intensity *= exp( -Biso*qval*qval/EIGHTPISQUARED ); */

	if ( debug > 2 && intensity > ptol ) {
	  printf("#  h  k  l    ( intensity > TOLERANCE )\n");
	  printf("# %2d %2d %2d : Q = %12.6f\t\tI = %12.6f\n\n",
		 h,k,l,qval,intensity);
	}
	if ( debug > 2 && intensity <= ptol ) {
	  printf("#  h  k  l   ( intensity <= TOLERANCE )\n");
	  printf("# %2d %2d %2d : Q = %12.6f\t\tI=%12.6f\n\n",
		 h,k,l,qval,intensity);
	}

	/***********************************************/
	/*  this code bins the intensities   */
	for ( q_n=0,peak_flag=0; q_n<Q_VALS_MAX; q_n++ ) {
	  /* check if this qval will fit in an existing bin */
	  /* peak intensity will go into first existing bin!! */
	  if ( tol( qval , pk[q_n].q , qtol ) ) {
	    /* if so, add the intensity and set the peak flag */
	    peak_flag=1;
	    pk[q_n].intensity += intensity;
	    pk[q_n].multiplicity++;
	    q_n=Q_VALS_MAX; /* break the loop */
	  }
	}
	/* if the qval has no bin, we must find a new one */
	/* these may overlap with previous bins */
	if ( peak_flag == 0 ) {
	  for ( q_n=0; q_n<Q_VALS_MAX; q_n++ ) {
	    /* continue through q_n until an empty bin is found */
	    if ( pk[q_n].q != -1 ) continue;
	    if ( pk[q_n].q == -1 ) {
	      /* these are initialization values for a new peak */
	      pk[q_n].q = qval;
	      pk[q_n].intensity = intensity;
	      pk[q_n].multiplicity = 1;
	      if ( LIST_HKL ){
		pk[q_n].h = h;
		pk[q_n].k = k;
		pk[q_n].l = l;
	      }
	    }
	    break;
	  }
	}
	/* end code that bins intensities */
	/*************************************************/

      }
    }
  }


  /* check to see if we ran out of peak bins */
  for ( q_n=0; q_n<Q_VALS_MAX; q_n++ ) {
    if ( pk[q_n].q == -1 ) have_remaining_bins = 1;
    if ( pk[q_n].q != -1 ) bins_used++;
  }
  if ( have_remaining_bins != 1 ) {
    printf("Ran out of peak bins!!\n");
    exit(0);
  }

  /* sort the peaks from low to high */
  qsort(pk,Q_VALS_MAX,sizeof(struct peak),
	(int (*)(const void *, const void *)) cmp);

  /* Output scaling and finding maximum intensity to cut off peaks */
  /* below (intensity/max_intensity) value   */
  for ( q_n=0; q_n<Q_VALS_MAX; q_n++ ) {

    pk[q_n].intensity *= output_scale;

    if ( pk[q_n].intensity > max_intensity ) 
      max_intensity=pk[q_n].intensity;
  }
  if ( debug > 2 ) printf("max_intensity = %f\n",max_intensity);  
  if ( debug > 2 ) printf("bins_used = %d\n",bins_used);  


  /* final output */
  for ( q_n=0; q_n<Q_VALS_MAX; q_n++ ) {

    if ( pk[q_n].q != -1 && pk[q_n].intensity/max_intensity > ptol ) {

      ttr = 2*asin(lambda*pk[q_n].q/(4*PI));
      tth = (180/PI)*ttr;
      /* the z-shift is derived in Fundamentals book 1, p 72 */
      /* the arm length on the Scintag is roughly 300 mm */
      tth += (180/PI)*(2*z_shift/300)*cos(tth*PI/360);

      /* if the tth value passes 180 degrees the calculation of tth
	 will return nan */
      if ( !finite(tth) ) continue;

      intensity = pk[q_n].intensity;
      h = pk[q_n].h;
      k = pk[q_n].k;
      l = pk[q_n].l;
      mult = pk[q_n].multiplicity;
      qval = pk[q_n].q;

      if ( D_COR == 1 ) intensity *= 
			(1+cos(ttr)*cos(ttr))/(2*sin(ttr/2)*2*PI*sin(ttr));

      if ( TWOTHETA && LIST_HKL ) {
	printf("%8.4f%20.6e%5d%5d%5d%12d\n",
	       tth,intensity,h,k,l,mult);
      } else if ( TWOTHETA && !LIST_HKL ) { 
	if ( debug > 2 ) printf("tth=%f Q=%f\n",tth,qval);
	printf("%8.4f%20.6e\n",tth,intensity);
      } else if ( LIST_HKL && !TWOTHETA ) {
	printf("%8.4f%20.6e%5d%5d%5d%12d\n",
	       qval,intensity,h,k,l,mult);
      } else {
	printf("%15.10f%15.10f%20.6e\n",qval,(2*PI/qval),intensity);
      }
    }
  }  
  

  return(0); 
}


/*******************************/
/*     Functions               */
/*******************************/
double StrFactor(
		 double *atom,
		 int numatoms,
		 int h,int k,int l,
		 double qval,
		 int biso)
{
  int i,osum;
  extern int XRAY,NEUT,debug;

  double qdotr;
  double I=0;
  double realsum=0,imagsum=0;
  double fact=0;
  double X;
  double x,y,z,Zat,Biso;
  double b_fact,zocc_fact;
  
  for ( i=0; i<numatoms; i++ ) {
    zocc_fact=1; /* unless set differently below */
    x    = *(atom+5*i+0);
    y    = *(atom+5*i+1);
    z    = *(atom+5*i+2);
    Zat  = *(atom+5*i+3);
    Biso = *(atom+5*i+4);

    /* see Warren page 37 for the biso factor */
    if ( biso == 0 ) b_fact = exp( - ( qval * qval / ( 16 * PI * PI ) ) );
    else b_fact = exp( - (Biso * qval * qval / ( 16 * PI * PI ) ) );

    qdotr = h * x + k * y + l * z;
    X = 2*PI*qdotr;

    /* if Zat was set to zero, it's scattering factor is set to zero */
    if ( Zat==0 ) { fact=0; }
    else {
      if ( XRAY ) fact = xfact(Zat,qval);
      if ( NEUT ) fact = nfact(Zat);
      if ( fact == -1 ) {
	printf("Can't find neutron scattering factor for Z = %6.3f\n",Zat);
	exit(0);
      }
    }
    for ( osum=1; osum < zocc_count+1 ; osum++ ){
      if ( zocc_Z[osum] == (int)Zat ) {
	//printf("osum=%d  Zat=%d  zocc_Z=%d  zocc_pct=%f\n",osum, (int)Zat, zocc_Z[osum],zocc_pct[osum]);
	zocc_fact = zocc_pct[osum];
      }
    }

    realsum += cos( X ) * fact * b_fact * zocc_fact;
    imagsum += sin( X ) * fact * b_fact * zocc_fact;

    if ( debug > 3 )
      printf("atom [%3d]: real = %11.6f\t   "
	     "imag = %11.6f\tq.r = %11.6f\n",i,cos(X)*fact,sin(X)*fact,X);
  }
  if ( debug > 1 ) printf("Totals:     REAL = %12.6f\t   IMAG = %12.6f\n",
			  realsum,imagsum);
  I = realsum*realsum + imagsum*imagsum;
  return (I);
}

double xfact(double z, double g)
{
  double a0 = 0.529, r, ans;
  
  r = a0 / pow( z, 1.0/3.0 );
  ans = z / ( 1.0 + g*g*r*r );
  
  return ans;
}

double nfact(double z)
{
  if ( z == -1  ) return (6.671);     /*  Deuterium   */
  if ( z == 1   ) return (-3.7406);     /*  H   */
  if ( z == 2   ) return (3.26);
  if ( z == 3   ) return (-1.90);       /*  Li  */
  if ( z == 4   ) return (7.79);
  if ( z == 5   ) return (6.65);
  if ( z == 6   ) return (6.65);
  if ( z == 7   ) return (9.36);        /*  N   */
  if ( z == 8   ) return (5.803);
  if ( z == 9   ) return (5.654);
  if ( z == 10  ) return (4.566);
  if ( z == 11  ) return (3.63);        /*  Na  */  
  if ( z == 12  ) return (5.375);
  if ( z == 13  ) return (3.45);        /*  Al  */
  if ( z == 14  ) return (4.15);
  if ( z == 15  ) return (5.13);
  if ( z == 16  ) return (2.847);
  if ( z == 17  ) return (9.577);
  if ( z == 18  ) return (1.909);
  if ( z == 19  ) return (3.67);
  if ( z == 20  ) return (4.70);
  if ( z == 21  ) return (12.29);
  if ( z == 22  ) return (-3.438);      /*  Ti  */
  if ( z == 23  ) return (-0.3824);
  if ( z == 24  ) return (3.636);
  if ( z == 25  ) return (-3.73);
  if ( z == 26  ) return (9.45);
  if ( z == 27  ) return (2.49);
  if ( z == 28  ) return (10.3);
  if ( z == 29  ) return (7.718);
  if ( z == 30  ) return (5.680);
  if ( z == 31  ) return (7.288);
  if ( z == 32  ) return (8.185);
  if ( z == 33  ) return (6.58);
  if ( z == 34  ) return (7.97);
  if ( z == 35  ) return (6.795);
  if ( z == 36  ) return (7.81);
  if ( z == 37  ) return (7.09);
  if ( z == 38  ) return (7.02);
  if ( z == 39  ) return (7.75);
  if ( z == 40  ) return (7.16);
  if ( z == 41  ) return (7.054);
  if ( z == 42  ) return (6.715);
  if ( z == 43  ) return (6.8);
  if ( z == 44  ) return (7.03);
  if ( z == 45  ) return (5.88);
  if ( z == 46  ) return (5.91);
  if ( z == 47  ) return (5.922);
  if ( z == 48  ) return (4.87);
  if ( z == 49  ) return (4.065);
  if ( z == 50  ) return (6.225);
  if ( z == 51  ) return (5.57);
  if ( z == 52  ) return (5.80);
  if ( z == 53  ) return (5.28);
  if ( z == 54  ) return (4.92);
  if ( z == 55  ) return (5.42);
  if ( z == 56  ) return (5.07);
  if ( z == 57  ) return (8.24);
  if ( z == 58  ) return (4.84);
  if ( z == 59  ) return (4.58);
  if ( z == 60  ) return (7.69);
  if ( z == 61  ) return (12.6);
  if ( z == 62  ) return (0.80);
  if ( z == 63  ) return (7.22);
  if ( z == 64  ) return (6.5);
  if ( z == 65  ) return (7.38);
  if ( z == 66  ) return (16.9);
  if ( z == 67  ) return (8.01);
  if ( z == 68  ) return (7.79);
  if ( z == 69  ) return (7.07);
  if ( z == 70  ) return (12.43);
  if ( z == 71  ) return (7.21);
  if ( z == 72  ) return (7.7);
  if ( z == 73  ) return (6.91);
  if ( z == 74  ) return (4.86);
  if ( z == 75  ) return (9.2);
  if ( z == 76  ) return (10.7);
  if ( z == 77  ) return (10.6);
  if ( z == 78  ) return (9.6);
  if ( z == 79  ) return (7.63);
  if ( z == 80  ) return (12.692);
  if ( z == 81  ) return (8.776);
  if ( z == 82  ) return (9.405);
  if ( z == 83  ) return (8.532);
  if ( z == 88  ) return (10.0);
  if ( z == 90  ) return (10.31);
  if ( z == 91  ) return (9.1);
  if ( z == 92  ) return (8.417);
  if ( z == 93  ) return (10.55);
  if ( z == 95  ) return (8.3);
  return (-1);
}


/* this is for a full triclinic lattice */
double q(
	 double a,
	 double b,
	 double c,
	 double alpha,
	 double beta,
	 double gamma,
	 int h,
	 int k,
	 int l)
{
  static int been_called = 0;

  static double arads,brads,grads;
  static double sa,ca,sb,cb,sg,cg,denom;

  double term1,term2,term3,term4,term5,term6;
  double sqdinv,q;

  
  if ( been_called == 0 )
    {
      arads = PI * alpha / 180;
      brads = PI * beta / 180;
      grads = PI * gamma / 180;
      sa = sin(arads);
      ca = cos(arads);
      sb = sin(brads);
      cb = cos(brads);
      sg = sin(grads);
      cg = cos(grads);
      denom = 1/(1+2*ca*cb*cg-ca*ca-cb*cb-cg*cg);
      
      been_called = 1;
    }
  
  term1 = h*h*sa*sa/(a*a);
  term2 = k*k*sb*sb/(b*b);
  term3 = l*l*sg*sg/(c*c);
  term4 = 2*h*k*(ca*cb-cg)/(a*b);
  term5 = 2*k*l*(cb*cg-ca)/(b*c);
  term6 = 2*l*h*(cg*ca-cb)/(a*c);
  
  sqdinv = ( term1 + term2 + term3 + term4 + term5 + term6 ) * denom;
  q = 2*PI*sqrt(sqdinv);
  return(q);
}

int cmp (double *q1, double *q2)
{
  if (*q1 < *q2) return -1;
  else if (*q1 == *q2) return 0;
  return 1;
}
