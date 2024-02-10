/*
 * Common C library
 * Eric Majzoub
 * Sandia National Laboratories
 * Livermore, CA
 * 10 October 2001
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include "ericlib.h"

#define PI 3.14159265358979
#define RESCALE_TOL 1e-10

/*
 * Tol: tolerance function.
 * Returns: 1 -- if inside tolerance
 *          0 -- outside tolerance
 */
int tol(double a,double b,double t)
{
  if ( fabs(a-b) < t ) return (1);
  else return(0);
}

/*
 * Random sign function.
 * Returns +1 or -1 randomly
 *
 */
int rsign(void)
{
  if ( drand48() < 0.5 ) { return(1); }
  else return (-1);
}

/* dmod: double modulus function */
double dmod( double a, double b)
{

  if ( a<0.0 ) a = -a;
  if ( b<0.0 ) b = -b;

  while( (a-b)>b ){
    a -= b;
  }

  if ( (a-b) < 0.0 ) return(a);
  else if ( (a-b) > 0.0 ) return( (a-b) );
  else return(-1);

}

/************************************/
/*                                  */
/*        Complex Number            */
/*         Functions                */
/*                                  */
/************************************/

/*
 * cmag: magnitude
 *
 */
double cmag(struct cmplx a)
{
  return ( a.real*a.real + a.imag*a.imag );
}

/*
 * cmult: complex multiply
 *
 */
struct cmplx cmult(struct cmplx a,struct cmplx b)
{
  struct cmplx c;
  c.real = a.real*b.real - a.imag*b.imag;
  c.imag = a.real*b.imag + a.imag*b.real;
  return(c);
}

/*
 * cadd: complex add
 *
 */
struct cmplx cadd(struct cmplx a,struct cmplx b)
{
  struct cmplx c;
  c.real = a.real + b.real;
  c.imag = a.imag + b.imag;
  return(c);
}

/*
 * cadd: complex sub
 *
 */
struct cmplx csub(struct cmplx a,struct cmplx b)
{
  struct cmplx c;
  c.real = a.real - b.real;
  c.imag = a.imag - b.imag;
  return(c);
}

/*****************************/
/*                           */
/*     Vector                */
/*     Operations            */
/*                           */
/*****************************/


/*
 * makevec: make a vector out of doubles
 */
struct vector makevec(double a, double b, double c)
{
  struct vector v;

  v.x = a;
  v.y = b;
  v.z = c;

  return(v);
}

/*
 * vmag: vector magnitude
 */
double vmag(struct vector a)
{
  return ( sqrt(a.x*a.x + a.y*a.y + a.z*a.z) );
}
/*
 * vmagsq: square of the vector magnitude
 */
double vmagsq(struct vector a)
{
  return ( a.x*a.x + a.y*a.y + a.z*a.z );
}


/**************************************/
/*                                    */
/*            vsub                    */
/*                                    */
/**************************************/
struct vector vsub(struct vector a, struct vector b)
  {
  struct vector c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  return(c);
  }

/**************************************/
/*                                    */
/*        vadd                        */
/*                                    */
/**************************************/
struct vector vadd(struct vector a, struct vector b)
  {
  struct vector c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  return(c);
  }

/*
 * vector scalar multiply
 * vsmult
 */
struct vector vsmult(double a, struct vector b)
{
  struct vector c;
  
  c.x = a * b.x;
  c.y = a * b.y;
  c.z = a * b.z;
  return(c);
}

/*
 * vdotprod
 */
double vdotprod(struct vector a, struct vector b)
{
  return( a.x * b.x + a.y * b.y + a.z * b.z );
}

/*
 * vcross
 *
 */
struct vector crossprod(struct vector a, struct vector b)
  {
  struct vector c;

  c.x = a.y*b.z - a.z*b.y;
  c.y = a.z*b.x - a.x*b.z;
  c.z = a.x*b.y - a.y*b.x;
  return(c);
  }

/*
 * vnorm
 *
 */
struct vector normvec(struct vector a)
  {
  double d;
  struct vector c;

  d = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
  c.x = a.x/d;
  c.y = a.y/d;
  c.z = a.z/d;

  return(c);
  }

/*
 * vtriple
 */
double vtriple(struct vector a, struct vector b, struct vector c)
{
  double v;
  v = fabs( vdotprod( c, crossprod( a, b ) ) );
  return(v);
}

/*
 * Recip: return reciprocal space vectors
 */
struct basis recip(struct basis bas)
{
  struct basis rec;
  double vol,tpi;

  /*  tpi = 8.0*atan(1.0); */
  tpi = 1.0;
  vol = vtriple( bas.A, bas.B, bas.C );

  /*
  printf("vol=%f\n",vol);
  printf("%f %f %f\n",bas.A.x,bas.A.y,bas.A.z);
  printf("%f %f %f\n",bas.B.x,bas.B.y,bas.B.z);
  printf("%f %f %f\n",bas.C.x,bas.C.y,bas.C.z);
  */

  rec.A = vsmult( tpi/vol, crossprod( bas.B, bas.C ) );
  rec.B = vsmult( tpi/vol, crossprod( bas.C, bas.A ) );
  rec.C = vsmult( tpi/vol, crossprod( bas.A, bas.B ) );

  /*
  printf("%f %f %f\n",rec.A.x,rec.A.y,rec.A.z);
  printf("%f %f %f\n",rec.B.x,rec.B.y,rec.B.z);
  printf("%f %f %f\n",rec.C.x,rec.C.y,rec.C.z);
  */

  return(rec);
}

/**************************************/
/*                                    */
/* distance between two points        */
/*                                    */
/**************************************/
double dist(struct vector *a,struct vector *b)
   {
   double c;
   double ax,ay,az,bx,by,bz;

   ax = (*a).x; ay = (*a).y; az = (*a).z;
   bx = (*b).x; by = (*b).y; bz = (*b).z;

   c = sqrt((ax-bx)*(ax-bx)+
	    (ay-by)*(ay-by)+
	    (az-bz)*(az-bz));
   return(c);
   }


/**************************************/
/*                                    */
/*       facegen                      */
/*                                    */
/**************************************/
/*
 * Note: the order of the vectors is important to the direction of the
 *       face vector.  You want it to point outward from the center
 *       of the polygon.  So choose a,b,c such that the right hand
 *       rule has it pointing outward.
 */
struct face facegen(struct vector *a,struct vector *b,struct vector *c)
   {
   double ax,ay,az,bx,by,bz,cx,cy,cz;
   struct vector cent,mp;
   struct vector v1,v2,v3,normal;
   struct face fce;

/*   printf("facegen\n"); */

   ax = (*a).x; ay = (*a).y; az = (*a).z;
   bx = (*b).x; by = (*b).y; bz = (*b).z;
   cx = (*c).x; cy = (*c).y; cz = (*c).z;


   mp.x = ax + (bx-ax)/2;  /* find center of face */
   mp.y = ay + (by-ay)/2;
   mp.z = az + (bz-az)/2;
   cent.x = mp.x + (cx-mp.x)/3;
   cent.y = mp.y + (cy-mp.y)/3;
   cent.z = mp.z + (cz-mp.z)/3;

   /* make vecs from vert pts for normal */
   v1 = vsub(makevec(bx,by,bz),makevec(ax,ay,az));
   v2 = vsub(makevec(cx,cy,cz),makevec(ax,ay,az));

   v3 = crossprod(v1,v2);

   normal = normvec(v3);

   fce.vert1.x = ax;
   fce.vert1.y = ay;
   fce.vert1.z = az;

   fce.vert2.x = bx;
   fce.vert2.y = by;
   fce.vert2.z = bz;

   fce.vert3.x = cx;
   fce.vert3.y = cy;
   fce.vert3.z = cz;
   
   fce.center.x = cent.x;
   fce.center.y = cent.y;
   fce.center.z = cent.z;

   fce.normal.x = normal.x;
   fce.normal.y = normal.y;
   fce.normal.z = normal.z;

   return(fce);
   
   }





/**************************************/
/*                                    */
/*            printface               */
/*                                    */
/**************************************/
void printface(struct face f)
   {

   printf("FACE:\n");
     printf("(%2.2f,%2.2f,%2.2f)"
	    " (%2.2f,%2.2f,%2.2f)"
	    " (%2.2f,%2.2f,%2.2f)"
	    " (%2.2f,%2.2f,%2.2f)\n",
	    f.vert1.x,
	    f.vert1.y,
	    f.vert1.z,

	    f.vert2.x,
	    f.vert2.y,
	    f.vert2.z,

	    f.vert3.x,
	    f.vert3.y,
	    f.vert3.z,
	    
	    f.normal.x,
	    f.normal.y,
	    f.normal.z);

   }


/*********************/
/* MATRIX OPERATIONS */
/*********************/
double det3x3(struct matrix R)
{
  return( R.r11*(R.r22*R.r33-R.r23*R.r32)-
	  R.r12*(R.r21*R.r33-R.r23*R.r31)+
	  R.r13*(R.r21*R.r32-R.r31*R.r22)  );
}

/* create a rotation matrix for rotaions away from the standard */
/* cartesian coordinate system.  These follow Byron and Fuller's convention */
/* ranges: theta:0..Pi   phi:0..2Pi  psi:0..2Pi */
struct matrix R_ptp(double phi,double the,double psi)
{
  double sphi,cphi,sthe,cthe,spsi,cpsi;
  struct matrix R;

  sphi=sin(phi); cphi=cos(phi);
  sthe=sin(the); cthe=cos(the);
  spsi=sin(psi); cpsi=cos(psi);

  R.r11 =  cphi * cthe * cpsi - sphi * spsi;
  R.r12 =  sphi * cthe * cpsi + cphi * spsi;
  R.r13 = -sthe * cpsi;

  R.r21 = -cphi * cthe * spsi - sphi * cpsi;
  R.r22 = -sphi * cthe * spsi + cphi * cpsi;
  R.r23 =  sthe * spsi;

  R.r31 = cphi * sthe;
  R.r32 = sphi * sthe;
  R.r33 = cthe;

  return R;
}
struct matrix R_ptp_inv(struct matrix Rf)
{
  double r[9];
  struct matrix R;

  r[0] = Rf.r11;
  r[1] = Rf.r12;
  r[2] = Rf.r13;

  r[3] = Rf.r21;
  r[4] = Rf.r22;
  r[5] = Rf.r23;

  r[6] = Rf.r31;
  r[7] = Rf.r32;
  r[8] = Rf.r33;

  invt_matrx(3,r);

  R.r11 = r[0];
  R.r12 = r[1];
  R.r13 = r[2];

  R.r21 = r[3];
  R.r22 = r[4];
  R.r23 = r[5];

  R.r31 = r[6];
  R.r32 = r[7];
  R.r33 = r[8];

  return R;
}


/**************************************/
/*                                    */
/*    Matrix on Vector  multiply      */
/*                                    */
/**************************************/
struct vector multiply(struct matrix *R,struct vector *v)
{
  double vx,vy,vz;
  double R11,R12,R13,R21,R22,R23,R31,R32,R33;
  struct vector temp;
  
  vx = (*v).x; vy = (*v).y; vz = (*v).z;
  
  R11 = (*R).r11; R12 = (*R).r12; R13 = (*R).r13;
  R21 = (*R).r21; R22 = (*R).r22; R23 = (*R).r23;
  R31 = (*R).r31; R32 = (*R).r32; R33 = (*R).r33;
  
  temp.x = R11*vx + R12*vy + R13*vz;
  temp.y = R21*vx + R22*vy + R23*vz;
  temp.z = R31*vx + R32*vy + R33*vz;
  
  return(temp);
}

/**************************************/
/*                                    */
/*  NxN Matrix on Matrix  multiply    */
/*                                    */
/**************************************/
/*
 * This function takes:
 *
 * int N, size of NxN matrices
 * double *mat1, pointer to mat1
 * double *mat2, pointer to mat2
 * double *mat3, pointer to mat3
 */
void NxNmult(int N,double *mat1,double *mat2,double *mat3)
{
  int k,l,r;
  double sum;

  printf("Test this function before use:\n");
  for(k=0;k<N;k++){      /* row of mat1 */
    for(l=0;l<N;l++){    /* col of mat2 */
      for(sum=0,r=0;r<N;r++){
	sum += mat1[N*k+r] * mat2[N*r+l];
      }
      mat3[N*k+l]=sum;
    }
  }

  return;
}
/**************************************/
/*                                    */
/*  NxM Matrix on N Vector multiply   */
/*                                    */
/**************************************/
void NxMmat_Nvec_mult(int N,int M,double *mat,double *v)
{
  int i,j;
  double sum;
  double *tmp=NULL;

  tmp = (double *)malloc(N*sizeof(double));
  if ( !tmp ){
    printf("Error mallocing for tmp.\n"); exit(0);
  }

  for(i=0;i<N;i++){
    for(sum=0,j=0;j<M;j++){
      sum += mat[M*i+j] * v[j];
    }
    tmp[i]=sum;
  }

  for(i=0;i<N;i++) v[i]=tmp[i];

  free(tmp);
  return;
}

/**************************************/
/*                                    */
/*  Create Rotation Matrix Forward    */
/*                                    */
/**************************************/
struct matrix matrix_for(struct vector n1,struct vector n2)
{

  /* Note: This requires only two orthogonal vectors, it makes the third
     from the first two by the cross product. */
  
  struct vector e1hat,e2hat,e3hat;
  struct vector n1hat,n2hat,n3hat;
  struct matrix R;

  e1hat.x = 1;   e2hat.x = 0;   e3hat.x = 0; 
  e1hat.y = 0;   e2hat.y = 1;   e3hat.y = 0;
  e1hat.z = 0;   e2hat.z = 0;   e3hat.z = 1;

  /* orientation decided by the input vectors definition of a new frame. */
  /* use grahm-schmidt */
  n1hat = normvec(n1);
  n2hat = normvec( vsub( n2, vsmult( vdotprod(n1,n2), n1 ) ) );
  n3hat = normvec( crossprod(n1hat,n2hat) );  
  
  /* define rotation matrix elements */
  R.r11 = vdotprod(e1hat,n1hat);
  R.r12 = vdotprod(e1hat,n2hat);
  R.r13 = vdotprod(e1hat,n3hat);
  
  R.r21 = vdotprod(e2hat,n1hat);
  R.r22 = vdotprod(e2hat,n2hat);
  R.r23 = vdotprod(e2hat,n3hat);
  
  R.r31 = vdotprod(e3hat,n1hat);
  R.r32 = vdotprod(e3hat,n2hat);
  R.r33 = vdotprod(e3hat,n3hat);
  
  return(R);
}


/*
          This function takes an int, and the pointer to the matrix to
          be inverted.  It returns void.

          This program reduces the given matrix to identity, while at
          the same time performing the same operation to an identity.
          When the given reaches identity, the identity reaches givens
          inverse.
   
          This program only uses partial pivoting (of the rows).
          Full pivoting may be added later along with other
          modifications.
          Remember to free id after you use it in whatever program
          calls this function.
 */
void invt_matrx(int n, double *matrx)
{
    int i, j, k, l, m, p, r, s;
    double *tempm=NULL, *tempi=NULL,*id=NULL;

    tempm = (double *)malloc( n * sizeof( double ) );
    if ( !tempm ) { printf("error: invt_matrx tempm malloc\n"); exit(0); }

    tempi = (double *)malloc( n * sizeof( double ) );
    if ( !tempi ) { printf("error: invt_matrx tempi malloc\n"); exit(0); }

    id = (double *)malloc( n * n * sizeof( double ) );
    if ( !id ) { printf("error: invt_matrx id malloc\n"); exit(0); }

    /* create nxn identity */
    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    if (i == j) {
		id[n * i + j] = 1;
	    } else {
		id[n * i + j] = 0;
	    }
	}
    }

    for (l = 0; l < n; l++) {
	k = l;

	/* in case diagonal element is 0 */
	if (matrx[n * k + l] == 0) {
	    for (r = 0; r < n; r++) {
		if (matrx[n * r + l] != 0) {
		    break;
		}
	    }
	    for (s = 0; s < n; s++) {
		matrx[n * k + s] += matrx[n * r + s];
		id[n * k + s] += id[n * r + s];
	    }
	}
	/* get a 1 in row l, column l */
	for (m = 0; m < n; m++) {
	    tempm[m] = matrx[n * k + m] / matrx[n * k + l];
	    tempi[m] = id[n * k + m] / matrx[n * k + l];
	}
	for (m = 0; m < n; m++) {
	    matrx[n * k + m] = tempm[m];
	    id[n * k + m] = tempi[m];
	}

	/* eliminate other elements in lth column */
	for (i = 0; i < n; i++) {
	    if (i == k && k == n - 1) {
		break;
	    }
	    if (i == k && k != n - 1) {
		i++;
	    }
	    /* mult kth row by ith column elem */
	    for (p = 0; p < n; p++) {
		tempm[p] = matrx[n * k + p] * matrx[n * i + l];
		tempi[p] = id[n * k + p] * matrx[n * i + l];
	    }

	    /* perform the row op */
	    for (p = 0; p < n; p++) {
		matrx[n * i + p] -= tempm[p];
		id[n * i + p] -= tempi[p];
	    }
	}
    }

    /* copy the new matrix back to the old positions and free memory */
    for (i=0; i<n*n; i++) matrx[i] = id[i];

    free(tempi);
    free(tempm);
    free(id);
    return;
}


/*
 * make_tetr
 * 
 * Given the center vector, two vectors defining another frame in E^3, and the
 * distance to the verices from the center, will return a 'tetrahedron'
 * rotated with the bottom face of the 'standard' tetrahedron pointing
 * toward z-hat direction of the new frame.
 *
 * Some of this code is taken from my awk program tetrahedron.awk
 *
 * 1. take a tetrahedron (standard defined in the code)
 * 2. translate the center to the origin
 * 3. scale the center to vertex distances (given d)
 * 3.5 scale the distances of the vertexes to give fractional coords
 * 4. rotate it approprately (given n1,n2)
 * 5. translate it out to the position (given c)
 * 6. return the tetrahedron
 *
 */
struct tetrahedron make_tetr(struct vector c,
			     double phi,
			     double the,
			     double psi,
			     double d,
			     struct cellprm cell)
{
  int i;
  int num_vert=5;
  int num_faces=4;
  int debug=0;    /* only for very low level debugging */

  struct face fce[4];
  struct vector w[5];
  struct tetrahedron tetr;
  
  struct vector centroid;
  
  struct matrix R;

  /* input debug */
  /* printf(" c = %lf %lf %lf\n",c.x,c.y,c.z); */
  
  /* initialize the default tetrahedron, from diffraction book 1, 1997 */
  w[1].x = 0;
  w[1].y = 1/SRT3;
  w[1].z = 0;
  
  w[2].x = -0.5;
  w[2].y = -1/(2*SRT3);
  w[2].z = 0;
  
  w[3].x = 0.5;
  w[3].y = -1/(2*SRT3);
  w[3].z = 0;
  
  w[4].x = 0;
  w[4].y = 0;
  w[4].z = SRT2/SRT3;
  
  /* the centroid */
  centroid.x = w[0].x = 0;
  centroid.y = w[0].y = 0;
  centroid.z = w[0].z = 1/(2*SRT6);
  
  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_tetr: pre tran atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* translate the centroid (and whole tetr) to the origin */
  for (i=0; i<num_vert; i++) {
    w[i].x -= centroid.x;
    w[i].y -= centroid.y;
    w[i].z -= centroid.z;
  }
  /* scale the tetrahedral vertex lengths (from centroid) */
  for (i=0; i<num_vert; i++) {
    w[i].x *= (d * 2*SRT2/SRT3);
    w[i].y *= (d * 2*SRT2/SRT3);
    w[i].z *= (d * 2*SRT2/SRT3);
  }

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_tetr: pre rot atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }
  
  /*
   * rotate the tetrahedron so the 'bottom' face normal is pointing along
   * the direction of z-hat in the new frame.  This will be the default
   * since the face normal for the 'bottom' face points along -zhat.
  */


  /* orientation decided by random 3 numbers for phi,theta,psi */
  R = R_ptp(phi,the,psi);
  
  /* Rotate the vertex points on the centered polygon into the nhat frame */
  for (i=0; i<num_vert; i++) w[i] = multiply(&R,&w[i]);

  /* scale the tetrahedral vertex lengths (to frac coords) */
  for (i=1; i<num_vert; i++) {
    w[i].x /= cell.a;
    w[i].y /= cell.b;
    w[i].z /= cell.c;
  }

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_tetr: post rot atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* Translate the polygon */
  for (i=0; i<num_vert; i++) w[i] = rezone( vadd(w[i],c) );

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_tetr: post tran atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* create the face vectors */
  fce[0] = facegen(&w[1],&w[2],&w[3]);
  fce[1] = facegen(&w[1],&w[3],&w[4]);
  fce[2] = facegen(&w[1],&w[4],&w[2]);
  fce[3] = facegen(&w[2],&w[3],&w[4]);

  if (debug) for (i=0; i<num_faces; i++){
    printf("f:make_tetr: face [%d] %15.9lf%15.9lf%15.9lf\n",i,fce[i].normal.x,
	   fce[i].normal.y,fce[i].normal.z);
  }

  /* make the tetrahedron structure */
  for (i=0; i<num_vert;  i++) tetr.v[i] = w[i];
  for (i=0; i<num_faces; i++) tetr.f[i] = fce[i];

  return(tetr);

}



/*
 * make_octa
 * 
 * Given the center vector, two vectors defining another frame in E^3, and the
 * distance to the verices from the center, will return an 'octahedron'.
 *
 *
 * 1. take an octahedron (standard defined in the code)
 * 2. translate the center to the origin
 * 3. scale the center to vertex distances (given d)
 * 3.5 scale the distances of the vertexes to give fractional coords
 * 4. rotate it approprately (given n1,n2)
 * 5. translate it out to the position (given c)
 * 6. return the octaahedron
 *
 */
struct octahedron make_octa(struct vector c,
			    double phi,
			    double the,
			    double psi,
			    double d,
			    struct cellprm cell)
{
  int i;
  int num_vert=7;
  int num_faces=8;
  int debug=0;    /* only for very low level debugging */

  struct face fce[8];
  struct vector w[7];
  struct octahedron octa;
  
  struct vector centroid;

  struct matrix R;

  /* input debug */
  /* printf(" c = %lf %lf %lf\n",c.x,c.y,c.z); */
  
  /* initialize the default octahedron */
  /* the centroid */
  centroid.x = w[0].x = 0;
  centroid.y = w[0].y = 0;
  centroid.z = w[0].z = 0;

  w[1].x =  1;
  w[1].y =  0;
  w[1].z =  0;
  
  w[2].x = -1;
  w[2].y =  0;
  w[2].z =  0;
  
  w[3].x =  0;
  w[3].y =  1;
  w[3].z =  0;
  
  w[4].x =  0;
  w[4].y = -1;
  w[4].z =  0;

  w[5].x =  0;
  w[5].y =  0;
  w[5].z =  1;
  
  w[6].x =  0;
  w[6].y =  0;
  w[6].z = -1;
  
  
  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_octa: pre tran atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* translate the centroid (and whole polygon) to the origin */
  for (i=0; i<num_vert; i++) {
    w[i].x -= centroid.x;
    w[i].y -= centroid.y;
    w[i].z -= centroid.z;
  }
  /* scale the polygon vertex lengths (from centroid) */
  for (i=0; i<num_vert; i++) {
    w[i].x *= d;
    w[i].y *= d;
    w[i].z *= d;
  }

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_octa: pre rot atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }
  
  /*
   * rotate the octahedron 
  */

  /* orientation decided by random 3 numbers for phi,theta,psi */
  R = R_ptp(phi,the,psi);
  
  /* Rotate the vertex points on the centered polygon into the nhat frame */
  for (i=0; i<num_vert; i++) w[i] = multiply(&R,&w[i]);

  /* scale the tetrahedral vertex lengths (to frac coords) */
  for (i=1; i<num_vert; i++) {
    w[i].x /= cell.a;
    w[i].y /= cell.b;
    w[i].z /= cell.c;
  }

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_octa: post rot atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* Translate the polygon  */
  for (i=0; i<num_vert; i++) w[i] = rezone( vadd(w[i],c) );

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_octa: post tran atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* create the face vectors */
  fce[0] = facegen(&w[1],&w[3],&w[2]);
  fce[1] = facegen(&w[1],&w[4],&w[3]);
  fce[2] = facegen(&w[1],&w[5],&w[4]);
  fce[3] = facegen(&w[1],&w[2],&w[5]);
  fce[4] = facegen(&w[6],&w[2],&w[3]);
  fce[5] = facegen(&w[6],&w[3],&w[4]);
  fce[6] = facegen(&w[6],&w[4],&w[5]);
  fce[7] = facegen(&w[6],&w[5],&w[2]);

  if (debug) for (i=0; i<num_faces; i++){
    printf("f:make_octa: face [%d] %15.9lf%15.9lf%15.9lf\n",i,fce[i].normal.x,
	   fce[i].normal.y,fce[i].normal.z);
  }

  /* make the octahedron structure */
  for (i=0; i<num_vert;  i++) octa.v[i] = w[i];
  for (i=0; i<num_faces; i++) octa.f[i] = fce[i];

  return(octa);

}





/**************************************************************/
/*                                                            */
/* Functions requiring the use of                             */
/* Periodic Boundary Conditions                               */
/*                                                            */
/**************************************************************/

/**************************************/
/*                                    */
/* Check that atoms are within bounds */
/*                                    */
/**************************************/
void check_cats(struct cation *t,int n)
{
  int i;
  double x,y,z;

  for (i=0; i<n; i++){
    x = (*(t+i)).center.x;
    y = (*(t+i)).center.y;
    z = (*(t+i)).center.z;
    if ( x > 1.0 || x < 0.0 ) printf("cat  [%d]: out of bounds X %15.9f\n",i,x);
    if ( y > 1.0 || y < 0.0 ) printf("cat  [%d]: out of bounds Y %15.9f\n",i,y);
    if ( z > 1.0 || z < 0.0 ) printf("cat  [%d]: out of bounds Z %15.9f\n",i,z);
  }
}

void check_tetr(struct tetrahedron *t,int n)
{
  int i,j;
  int num_vert=5;
  double x,y,z;

  for (i=0; i<n; i++){
    for (j=0; j<num_vert; j++) {
      x = (*(t+i)).v[j].x;
      y = (*(t+i)).v[j].y;
      z = (*(t+i)).v[j].z;
      if ( x > 1.0 || x < 0.0 ) printf("tetr [%d] v[%d]: out of bounds X %15.9f\n",i,j,x);
      if ( y > 1.0 || y < 0.0 ) printf("tetr [%d] v[%d]: out of bounds Y %15.9f\n",i,j,y);
      if ( z > 1.0 || z < 0.0 ) printf("tetr [%d] v[%d]: out of bounds Z %15.9f\n",i,j,z);
    }
  }
}

void check_octa(struct octahedron *t,int n)
{
  int i,j;
  int num_vert=7;
  double x,y,z;

  for (i=0; i<n; i++){
    for (j=0; j<num_vert; j++) {
      x = (*(t+i)).v[j].x;
      y = (*(t+i)).v[j].y;
      z = (*(t+i)).v[j].z;
      if ( x > 1.0 || x < 0.0 ) printf("tetr [%d] v[%d]: out of bounds X %15.9f\n",i,j,x);
      if ( y > 1.0 || y < 0.0 ) printf("tetr [%d] v[%d]: out of bounds Y %15.9f\n",i,j,y);
      if ( z > 1.0 || z < 0.0 ) printf("tetr [%d] v[%d]: out of bounds Z %15.9f\n",i,j,z);
    }
  }
}


/**************************************/
/*                                    */
/* Re-zone an atom                    */
/*                                    */
/**************************************/
struct vector rezone(struct vector v)
{
  if ( v.x > 1.0 ) v.x -= 1.0;
  if ( v.x < 0.0 ) v.x += 1.0;
  if ( v.x < -1.0 || v.x > 2.0 ) {
    printf("rezone out of range!\n");
    exit(0);
  }

  if ( v.y > 1.0 ) v.y -= 1.0;
  if ( v.y < 0.0 ) v.y += 1.0;
  if ( v.y < -1.0 || v.y > 2.0 ) {
    printf("rezone out of range!\n");
    exit(0);
  }

  if ( v.z > 1.0 ) v.z -= 1.0;
  if ( v.z < 0.0 ) v.z += 1.0;
  if ( v.z < -1.0 || v.z > 2.0 ) {
    printf("rezone out of range!\n");
    exit(0);
  }

  return(v);
}

/**************************************/
/*                                    */
/* Re-scale after a lat prm change    */
/*                                    */
/**************************************/
void rescale_tetr(struct tetrahedron *t, int n,struct cellprm cell, struct vector cell_chng)
{
  int i;
  int j;
  int num_vert=5;
  double centx,centy,centz;
  double vertx,verty,vertz;
  double xdiff,ydiff,zdiff;
  double newx=0,newy=0,newz=0;
  static double tolerance=RESCALE_TOL;

  for (i=0; i<n; i++) {
    centx = (*(t+i)).v[0].x;
    centy = (*(t+i)).v[0].y;
    centz = (*(t+i)).v[0].z;

    /* now go over j=1 to all verts and rescale */
    for (j=1; j<num_vert; j++) {
      vertx = (*(t+i)).v[j].x;
      verty = (*(t+i)).v[j].y;
      vertz = (*(t+i)).v[j].z;
      xdiff = (vertx - centx);
      ydiff = (verty - centy);
      zdiff = (vertz - centz);
      newx  = 0;
      newy  = 0;
      newz  = 0;
      
      /* There are three choices for the difference between the vertex
	 position and center position.  If the difference is negative
	 (1) or positive (2), then there are two cases for each.  If
	 there is no difference in position (3) , then scaling does
	 nothing, and you just return the vertex value.
      */

      /* x */
      if ( !tol( cell_chng.x, 0.0, tolerance ) ) {

	if ( xdiff > 0.0 ) {
	  if ( fabs(xdiff) > 0.5 ) { newx = centx - fabs( centx - vertx + 1.0 ) * cell.a / ( cell.a + cell_chng.x ); }
	  if ( fabs(xdiff) < 0.5 ) { newx = centx + fabs( centx - vertx ) * cell.a / ( cell.a + cell_chng.x ); }
	}
	if ( xdiff < 0.0 ) {
	  if ( fabs(xdiff) > 0.5 ) { newx = centx + fabs( vertx - centx + 1.0 ) * cell.a / ( cell.a + cell_chng.x ); }
	  if ( fabs(xdiff) < 0.5 ) { newx = centx - fabs( vertx - centx ) * cell.a / ( cell.a + cell_chng.x ); }
	}
	if ( tol( xdiff, 0.0, tolerance ) ) newx = vertx;
	if ( newx > 1.0 ) newx -= 1.0;
	if ( newx < 0.0 ) newx += 1.0;
	(*(t+i)).v[j].x = newx;

      }
      
      /* y */
      if ( !tol( cell_chng.y, 0.0, tolerance ) ) {

	if ( ydiff > 0.0 ) {
	  if ( fabs(ydiff) > 0.5 ) { newy = centy - fabs( centy - verty + 1.0 ) * cell.b / ( cell.b + cell_chng.y ); }
	  if ( fabs(ydiff) < 0.5 ) { newy = centy + fabs( centy - verty ) * cell.b / ( cell.b + cell_chng.y ); }
	}
	if ( ydiff < 0.0 ) {
	  if ( fabs(ydiff) > 0.5 ) { newy = centy + fabs( verty - centy + 1.0 ) * cell.b / ( cell.b + cell_chng.y ); }
	  if ( fabs(ydiff) < 0.5 ) { newy = centy - fabs( verty - centy ) * cell.b / ( cell.b + cell_chng.y ); }
	}
	if ( tol( ydiff, 0.0, tolerance ) ) newy = verty;
	if ( newy > 1.0 ) newy -= 1.0;
	if ( newy < 0.0 ) newy += 1.0;
	(*(t+i)).v[j].y = newy;

      }

      /* z */
      if ( !tol( cell_chng.z, 0.0, tolerance ) ) {

	if ( zdiff > 0.0 ) {
	  if ( fabs(zdiff) > 0.5 ) { newz = centz - fabs( centz - vertz + 1.0 ) * cell.c / ( cell.c + cell_chng.z ); }
	  if ( fabs(zdiff) < 0.5 ) { newz = centz + fabs( centz - vertz ) * cell.c / ( cell.c + cell_chng.z ); }
	}
	if ( zdiff < 0.0 ) {
	  if ( fabs(zdiff) > 0.5 ) { newz = centz + fabs( vertz - centz + 1.0 ) * cell.c / ( cell.c + cell_chng.z ); }
	  if ( fabs(zdiff) < 0.5 ) { newz = centz - fabs( vertz - centz ) * cell.c / ( cell.c + cell_chng.z ); }
	}
	if ( tol( zdiff, 0.0, tolerance ) ) newz = vertz;
	if ( newz > 1.0 ) newz -= 1.0;
	if ( newz < 0.0 ) newz += 1.0;
	(*(t+i)).v[j].z = newz;
	
      }   
    }
  }
  

  return;
}

void rescale_octa(struct octahedron *t, int n,struct cellprm cell, struct vector cell_chng)
{
  int i;
  int j;
  int num_vert=7;
  double centx,centy,centz;
  double vertx,verty,vertz;
  double xdiff,ydiff,zdiff;
  double newx=0,newy=0,newz=0;
  static double tolerance=RESCALE_TOL;
  
  for (i=0; i<n; i++) {
    centx = (*(t+i)).v[0].x;
    centy = (*(t+i)).v[0].y;
    centz = (*(t+i)).v[0].z;

    /* now go over j=1 to all verts and rescale */
    for (j=1; j<num_vert; j++) {
      vertx = (*(t+i)).v[j].x;
      verty = (*(t+i)).v[j].y;
      vertz = (*(t+i)).v[j].z;
      xdiff = (vertx - centx);
      ydiff = (verty - centy);
      zdiff = (vertz - centz);
      newx  = 0;
      newy  = 0;
      newz  = 0;
      
      /* There are three choices for the difference between the vertex
	 position and center position.  If the difference is negative
	 (1) or positive (2), then there are two cases for each.  If
	 there is no difference in position (3) , then scaling does
	 nothing, and you just return the vertex value.
      */

      /* x */
      if ( !tol( cell_chng.x, 0.0, tolerance ) ) {

	if ( xdiff > 0.0 ) {
	  if ( fabs(xdiff) > 0.5 ) { newx = centx - fabs( centx - vertx + 1.0 ) * cell.a / ( cell.a + cell_chng.x ); }
	  if ( fabs(xdiff) < 0.5 ) { newx = centx + fabs( centx - vertx ) * cell.a / ( cell.a + cell_chng.x ); }
	}
	if ( xdiff < 0.0 ) {
	  if ( fabs(xdiff) > 0.5 ) { newx = centx + fabs( vertx - centx + 1.0 ) * cell.a / ( cell.a + cell_chng.x ); }
	  if ( fabs(xdiff) < 0.5 ) { newx = centx - fabs( vertx - centx ) * cell.a / ( cell.a + cell_chng.x ); }
	}
	if ( tol( xdiff, 0.0, tolerance ) ) newx = vertx;
	if ( newx > 1.0 ) newx -= 1.0;
	if ( newx < 0.0 ) newx += 1.0;
	(*(t+i)).v[j].x = newx;

      }
      
      /* y */
      if ( !tol( cell_chng.y, 0.0, tolerance ) ) {

	if ( ydiff > 0.0 ) {
	  if ( fabs(ydiff) > 0.5 ) { newy = centy - fabs( centy - verty + 1.0 ) * cell.b / ( cell.b + cell_chng.y ); }
	  if ( fabs(ydiff) < 0.5 ) { newy = centy + fabs( centy - verty ) * cell.b / ( cell.b + cell_chng.y ); }
	}
	if ( ydiff < 0.0 ) {
	  if ( fabs(ydiff) > 0.5 ) { newy = centy + fabs( verty - centy + 1.0 ) * cell.b / ( cell.b + cell_chng.y ); }
	  if ( fabs(ydiff) < 0.5 ) { newy = centy - fabs( verty - centy ) * cell.b / ( cell.b + cell_chng.y ); }
	}
	if ( tol( ydiff, 0.0, tolerance ) ) newy = verty;
	if ( newy > 1.0 ) newy -= 1.0;
	if ( newy < 0.0 ) newy += 1.0;
	(*(t+i)).v[j].y = newy;

      }

      /* z */
      if ( !tol( cell_chng.z, 0.0, tolerance ) ) {

	if ( zdiff > 0.0 ) {
	  if ( fabs(zdiff) > 0.5 ) { newz = centz - fabs( centz - vertz + 1.0 ) * cell.c / ( cell.c + cell_chng.z ); }
	  if ( fabs(zdiff) < 0.5 ) { newz = centz + fabs( centz - vertz ) * cell.c / ( cell.c + cell_chng.z ); }
	}
	if ( zdiff < 0.0 ) {
	  if ( fabs(zdiff) > 0.5 ) { newz = centz + fabs( vertz - centz + 1.0 ) * cell.c / ( cell.c + cell_chng.z ); }
	  if ( fabs(zdiff) < 0.5 ) { newz = centz - fabs( vertz - centz ) * cell.c / ( cell.c + cell_chng.z ); }
	}
	if ( tol( zdiff, 0.0, tolerance ) ) newz = vertz;
	if ( newz > 1.0 ) newz -= 1.0;
	if ( newz < 0.0 ) newz += 1.0;
	(*(t+i)).v[j].z = newz;
	
      }   

    }
  }

  return;
}


/**************************************/
/*                                    */
/* translate a cation                 */
/* Periodic Boundary Conditions       */
/*                                    */
/**************************************/
void trans_cation(struct cation *t,struct vector a)
{
  (*t).center.x += a.x;
  if ( (*t).center.x > 1.0 ) (*t).center.x -= 1.0;
  if ( (*t).center.x < 0.0 ) (*t).center.x += 1.0;
  (*t).center.y += a.y;
  if ( (*t).center.y > 1.0 ) (*t).center.y -= 1.0;
  if ( (*t).center.y < 0.0 ) (*t).center.y += 1.0;
  (*t).center.z += a.z;
  if ( (*t).center.z > 1.0 ) (*t).center.z -= 1.0;
  if ( (*t).center.z < 0.0 ) (*t).center.z += 1.0;
  return;
}

/**************************************/
/*                                    */
/* translate a tetrahedron            */
/* Periodic Boundary Conditions       */
/*                                    */
/**************************************/
void trans_tetr(struct tetrahedron *t,struct vector a)
{
  int i,num_vert=5;

  for (i=0; i<num_vert; i++) {
    (*t).v[i].x += a.x;
    if ( (*t).v[i].x > 1.0 ) (*t).v[i].x -= 1.0;
    if ( (*t).v[i].x < 0.0 ) (*t).v[i].x += 1.0;
    (*t).v[i].y += a.y;
    if ( (*t).v[i].y > 1.0 ) (*t).v[i].y -= 1.0;
    if ( (*t).v[i].y < 0.0 ) (*t).v[i].y += 1.0;
    (*t).v[i].z += a.z;
    if ( (*t).v[i].z > 1.0 ) (*t).v[i].z -= 1.0;
    if ( (*t).v[i].z < 0.0 ) (*t).v[i].z += 1.0;
  }

  return;
}

/**************************************/
/*                                    */
/* translate an octahedron            */
/* Periodic Boundary Conditions       */
/*                                    */
/**************************************/
void trans_octa(struct octahedron *t,struct vector a)
{
  int i,num_vert=7;
  
  for (i=0; i<num_vert; i++) {
  (*t).v[i].x += a.x;
  if ( (*t).v[i].x > 1.0 ) (*t).v[i].x -= 1.0;
  if ( (*t).v[i].x < 0.0 ) (*t).v[i].x += 1.0;
  (*t).v[i].y += a.y;
  if ( (*t).v[i].y > 1.0 ) (*t).v[i].y -= 1.0;
  if ( (*t).v[i].y < 0.0 ) (*t).v[i].y += 1.0;
  (*t).v[i].z += a.z;
  if ( (*t).v[i].z > 1.0 ) (*t).v[i].z -= 1.0;
  if ( (*t).v[i].z < 0.0 ) (*t).v[i].z += 1.0;
  }
  return;
}


/**************************************/
/*                                    */
/* distance between two points        */
/* Periodic Boundary Conditions       */
/*                                    */
/**************************************/
double dist_pbc(struct vector *v,struct vector *w,struct cellprm *cell)
{

  double d;
  double xdiff,ydiff,zdiff;
  
  xdiff =  fabs( (*v).x - (*w).x );  
  if ( xdiff > 0.5 ) xdiff = 1.0 - xdiff;
  xdiff *= (*cell).a;

  ydiff = fabs( (*v).y - (*w).y );
  if ( ydiff > 0.5 ) ydiff = 1.0 - ydiff;
  ydiff *= (*cell).b;

  zdiff = fabs( (*v).z - (*w).z );
  if ( zdiff > 0.5 ) zdiff = 1.0 - zdiff;
  zdiff *= (*cell).c;
  
  d = sqrt( xdiff*xdiff + ydiff*ydiff + zdiff*zdiff );
  return(d);

}

/* the squared distance (faster) */
/* note: this is not the same function at all, it works with
   input in cart coords, NOT fractional!!!!
   This function was specificially written for Energy_rep() */
double distsq_pbc(struct vector *v,struct vector *w,struct cellprm *cell)
{

  double d;
  double xdiff,ydiff,zdiff;
  
  xdiff =  fabs( (*v).x - (*w).x )/(*cell).a;  
  if ( xdiff > 0.5 ) xdiff = 1.0 - xdiff;
  xdiff *= (*cell).a;

  ydiff = fabs( (*v).y - (*w).y )/(*cell).b;
  if ( ydiff > 0.5 ) ydiff = 1.0 - ydiff;
  ydiff *= (*cell).b;

  zdiff = fabs( (*v).z - (*w).z )/(*cell).c;
  if ( zdiff > 0.5 ) zdiff = 1.0 - zdiff;
  zdiff *= (*cell).c;
  
  d = ( xdiff*xdiff + ydiff*ydiff + zdiff*zdiff );
  return(d);

}


/***************************************/
/*                                     */
/*   Rotate a tetrahedron              */
/*                                     */
/***************************************/
void rot_tetr(struct tetrahedron *t,int n,double d,struct cellprm cell,struct matrix R)
{
  
  /* plan:
     1. translate the object to the origin
     2. scale to actual dimensions (makes verts equidistant from cent)
     3. perform rotation
     4. rescale to cell dimensions
     5. move object back to location it came from
     6. rezone the atoms if they went out of bounds
  */

  int i;
  int debug=0;
  int num_vert=5;
  struct vector c;

  c.x = t[n].v[0].x;
  c.y = t[n].v[0].y;
  c.z = t[n].v[0].z;

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: pre trans: atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* translate the object to the origin */
  for (i=0; i<num_vert; i++)  t[n].v[i] = rezone( vsub(t[n].v[i],c) );

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: post trans rezone: atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* scale the lengths to actual values and 'un-rezone' the positions of the atoms */
  for (i=1; i<num_vert; i++) {
    t[n].v[i].x *= cell.a;
    t[n].v[i].y *= cell.b;
    t[n].v[i].z *= cell.c;

    if ( t[n].v[i].x > d ) t[n].v[i].x -= cell.a;
    if ( t[n].v[i].y > d ) t[n].v[i].y -= cell.b;
    if ( t[n].v[i].z > d ) t[n].v[i].z -= cell.c;

    if ( !tol( vmag(t[n].v[i]), d, 1e-8 ) ) {
      printf("# rotation error: rezoning problem vertex %d\n",i);
      exit(0);
    }
  }

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: post scale: atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* Rotate the vertex points on the centered polygon into the nhat frame */
  for (i=1; i<num_vert; i++) t[n].v[i] = multiply(&R,&t[n].v[i]);

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: post rot atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* scale the tetrahedral vertex lengths (to frac coords) */
  for (i=1; i<num_vert; i++) {
    t[n].v[i].x /= cell.a;
    t[n].v[i].y /= cell.b;
    t[n].v[i].z /= cell.c;
  }

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: post scale atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* Translate the polygon */
  for (i=0; i<num_vert; i++) t[n].v[i] = rezone( vadd(t[n].v[i],c) );     

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: post trans rezone atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* create new face vectors */
  t[n].f[0] = facegen(&t[n].v[1],&t[n].v[2],&t[n].v[3]);
  t[n].f[1] = facegen(&t[n].v[1],&t[n].v[3],&t[n].v[4]);
  t[n].f[2] = facegen(&t[n].v[1],&t[n].v[4],&t[n].v[2]);
  t[n].f[3] = facegen(&t[n].v[2],&t[n].v[3],&t[n].v[4]);

  return;
}


/***************************************/
/*                                     */
/*   Rotate an octahedron              */
/*                                     */
/***************************************/
void rot_octa(struct octahedron *t,int n,double d,struct cellprm cell, struct matrix R)
{
  
  /* plan:
     1. translate the object to the origin
     2. scale to actual dimensions (makes verts equidistant from cent)
     3. perform rotation
     4. rescale to cell dimensions
     5. move object back to location it came from
     6. rezone the atoms if they went out of bounds
  */

  int i;
  int debug=0;
  int num_vert=7;
  struct vector c;

  c.x = t[n].v[0].x;
  c.y = t[n].v[0].y;
  c.z = t[n].v[0].z;

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: pre trans: atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* translate the object to the origin */
  for (i=0; i<num_vert; i++)  t[n].v[i] = rezone( vsub(t[n].v[i],c) );

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: post trans rezone: atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* scale the lengths to actual values and 'un-rezone' the positions of the atoms */
  for (i=1; i<num_vert; i++) {
    t[n].v[i].x *= cell.a;
    t[n].v[i].y *= cell.b;
    t[n].v[i].z *= cell.c;

    if ( t[n].v[i].x > d ) t[n].v[i].x -= cell.a;
    if ( t[n].v[i].y > d ) t[n].v[i].y -= cell.b;
    if ( t[n].v[i].z > d ) t[n].v[i].z -= cell.c;

    if ( !tol( vmag(t[n].v[i]), d, 1e-8 ) ) {
      printf("# rotation error: rezoning problem vertex %d\n",i);
      exit(0);
    }
  }

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: post scale: atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* Rotate the vertex points on the centered polygon into the nhat frame */
  for (i=1; i<num_vert; i++) t[n].v[i] = multiply(&R,&t[n].v[i]);

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: post rot atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* scale the tetrahedral vertex lengths (to frac coords) */
  for (i=1; i<num_vert; i++) {
    t[n].v[i].x /= cell.a;
    t[n].v[i].y /= cell.b;
    t[n].v[i].z /= cell.c;
  }

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: post scale atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* Translate the polygon */
  for (i=0; i<num_vert; i++) t[n].v[i] = rezone( vadd(t[n].v[i],c) );     

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:rot_tetr: post trans rezone atom C %15.9lf%15.9lf%15.9lf\n",t[n].v[i].x,t[n].v[i].y,t[n].v[i].z);
  }

  /* create new face vectors */
  t[n].f[0] = facegen(&t[n].v[1],&t[n].v[3],&t[n].v[2]);
  t[n].f[1] = facegen(&t[n].v[1],&t[n].v[4],&t[n].v[3]);
  t[n].f[2] = facegen(&t[n].v[1],&t[n].v[5],&t[n].v[4]);
  t[n].f[3] = facegen(&t[n].v[1],&t[n].v[2],&t[n].v[5]);
  t[n].f[4] = facegen(&t[n].v[6],&t[n].v[2],&t[n].v[3]);
  t[n].f[5] = facegen(&t[n].v[6],&t[n].v[3],&t[n].v[4]);
  t[n].f[6] = facegen(&t[n].v[6],&t[n].v[4],&t[n].v[5]);
  t[n].f[7] = facegen(&t[n].v[6],&t[n].v[5],&t[n].v[2]);

  return;
}


/***************************************/
/*                                     */
/*  PRINTING FUNCTIONS                 */
/*                                     */
/***************************************/

/**************************************/
/*                                    */
/*    print cell frame xbs            */
/*                                    */
/**************************************/
void print_cell_frame_xbs(struct cellprm cell,double FRM_RAD,int mvout)
{

  double HI;
  double LO;
  double CORNER_RAD=0.05;

  if ( mvout ) {
    HI = 100.0;
    LO = 0.00001;
  } else {
    HI = 1.001;
    LO = 0.999;
  }

  printf("* CELL PARAMETERS AND BORDER\n");
  printf("atom 0 %10.5f%10.5f%10.5f\n",0.0,   0.0,   0.0);
  printf("atom 1 %10.5f%10.5f%10.5f\n",cell.a,0.0,   0.0);
  printf("atom 2 %10.5f%10.5f%10.5f\n",cell.a,cell.b,0.0);
  printf("atom 3 %10.5f%10.5f%10.5f\n",0.0,   cell.b,0.0);

  printf("atom 4 %10.5f%10.5f%10.5f\n",0.0,   0.0,   cell.c);
  printf("atom 5 %10.5f%10.5f%10.5f\n",cell.a,0.0,   cell.c);
  printf("atom 6 %10.5f%10.5f%10.5f\n",cell.a,cell.b,cell.c);
  printf("atom 7 %10.5f%10.5f%10.5f\n",0.0,   cell.b,cell.c);

  print_spe_xbs('0',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('1',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('2',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('3',CORNER_RAD,0.0,0.0,0.0);

  print_spe_xbs('4',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('5',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('6',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('7',CORNER_RAD,0.0,0.0,0.0);

  print_bnd_xbs('0','1',cell.a*LO,cell.a*HI,FRM_RAD);
  print_bnd_xbs('3','2',cell.a*LO,cell.a*HI,FRM_RAD);
  print_bnd_xbs('4','5',cell.a*LO,cell.a*HI,FRM_RAD);
  print_bnd_xbs('6','7',cell.a*LO,cell.a*HI,FRM_RAD);

  print_bnd_xbs('1','2',cell.b*LO,cell.b*HI,FRM_RAD);
  print_bnd_xbs('0','3',cell.b*LO,cell.b*HI,FRM_RAD);
  print_bnd_xbs('5','6',cell.b*LO,cell.b*HI,FRM_RAD);
  print_bnd_xbs('4','7',cell.b*LO,cell.b*HI,FRM_RAD);

  print_bnd_xbs('0','4',cell.c*LO,cell.c*HI,FRM_RAD);
  print_bnd_xbs('1','5',cell.c*LO,cell.c*HI,FRM_RAD);
  print_bnd_xbs('2','6',cell.c*LO,cell.c*HI,FRM_RAD);
  print_bnd_xbs('3','7',cell.c*LO,cell.c*HI,FRM_RAD);

  printf("inc %f\n",5.0);
  if ( !mvout ) {
    printf("dup %10.5f%10.5f%10.5f\n",cell.a,0.0,0.0);
    printf("dup %10.5f%10.5f%10.5f\n",0.0,cell.b,0.0);
    printf("dup %10.5f%10.5f%10.5f\n",0.0,0.0,cell.c);
  }

  return;
}

/**************************************/
/*                                    */
/*    print cation                    */
/*                                    */
/**************************************/
void print_cat(struct cation t)
{
  printf("# CATION:\n");
  printf("#     %15.9lf%15.9lf%15.9lf\n",
	 t.center.x,
	 t.center.y,
	 t.center.z);
  return;
}

/**************************************/
/*                                    */
/*    print cation xbs                */
/*                                    */
/**************************************/
void print_cat_xbs(char c,struct cation t, struct cellprm cell)
{
  printf("* CATION:\n");
  printf("atom %c %15.9lf%15.9lf%15.9lf\n",
	 c,
	 cell.a * t.center.x,
	 cell.b * t.center.y,
	 cell.c * t.center.z);
  return;
}

/**************************************/
/*                                    */
/*    print tetrahedron               */
/*                                    */
/**************************************/
void print_tet(struct tetrahedron t)
{
  int i;
  int num_vert=5;
  int num_faces=4;
  
  printf("# TETRAHEDRON:\n");
  printf("# Vertices:\n");
  for (i=0; i<num_vert; i++) {
    printf("# [%d] %15.9lf%15.9lf%15.9lf\n",i,t.v[i].x,t.v[i].y,t.v[i].z);
  }
  printf("# Face normal vectors:\n");
  for (i=0; i<num_faces; i++) {
    printf("# [%d] %15.9lf%15.9lf%15.9lf\n",i,
	   t.f[i].normal.x,
	   t.f[i].normal.y,
	   t.f[i].normal.z);
  }
  return;
}

/**************************************/
/*                                    */
/*    print tetrahedron xbs           */
/*                                    */
/**************************************/
void print_tet_xbs(struct tetrahedron t,struct cellprm cell)
{
  int i;
  int num_vert=5;
  
  printf("* TETRAHEDRON:\n");
  printf("* Vertices:\n");
  printf("atom A %15.9lf%15.9lf%15.9lf\n",
	 cell.a * t.v[0].x,
	 cell.b * t.v[0].y,
	 cell.c * t.v[0].z);
  for (i=1; i<num_vert; i++) {
    printf("atom V %15.9lf%15.9lf%15.9lf\n",
	   cell.a * t.v[i].x,
	   cell.b * t.v[i].y,
	   cell.c * t.v[i].z);
  }
  return;
}


/**************************************/
/*                                    */
/*    print octahedron                */
/*                                    */
/**************************************/
void print_oct(struct octahedron t)
{
  int i;
  int num_vert=7;
  int num_faces=8;
  
  printf("# OCTAHEDRON:\n");
  printf("# Vertices:\n");
  for (i=0; i<num_vert; i++) {
    printf("# [%d] %15.9lf%15.9lf%15.9lf\n",i,t.v[i].x,t.v[i].y,t.v[i].z);
  }
  printf("# Face normal vectors:\n");
  for (i=0; i<num_faces; i++) {
    printf("# [%d] %15.9lf%15.9lf%15.9lf\n",i,
	   t.f[i].normal.x,
	   t.f[i].normal.y,
	   t.f[i].normal.z);
  }
  return;
}


/**************************************/
/*                                    */
/*    print octahedron xbs            */
/*                                    */
/**************************************/
void print_oct_xbs(struct octahedron t,struct cellprm cell)
{
  int i;
  int num_vert=7;
  
  printf("* OCTAHEDRON:\n");
  printf("* Vertices:\n");
  printf("atom A %15.9lf%15.9lf%15.9lf\n",
	 cell.a * t.v[0].x,
	 cell.b * t.v[0].y,
	 cell.c * t.v[0].z);
  for (i=1; i<num_vert; i++) {
    printf("atom V %15.9lf%15.9lf%15.9lf\n",
	   cell.a * t.v[i].x,
	   cell.b * t.v[i].y,
	   cell.c * t.v[i].z);
  }
  return;
}



/**************************************/
/*                                    */
/*  print species specs for xbs       */
/*                                    */
/**************************************/
void print_spe_xbs(char c,double R,double r,double g,double b)
{
  
  printf("spec %c %15.9lf%15.9lf%15.9lf%15.9lf\n",c,R,r,g,b);

  return;
}

/**************************************/
/*                                    */
/*    print bond specs for xbs        */
/*                                    */
/**************************************/
void print_bnd_xbs(char c1,char c2,double m,double M,double R)
{
  /* Note: all bonds will be black here */
  printf("bonds %c %c %12.9lf%15.9lf%15.9lf    %s\n",c1,c2,m,M,R,"Black");

  return;
}

/********************************************************/
/********************************************************/
/*             Ewald summation routines                 */
/********************************************************/
/********************************************************/
/* These routines were adapted from V.Ozolins fortran code */
/* They were not f2c'd. */
double find_Ewald_eta(struct basis bas,double ERRLIM,int nat,struct atom *at)
{
  int i,j;
  int iter;
  int MAXIT=100;
  int debug=0;

  static double eta=0;
  double pi,tpi,tpi2;
  double taumax,rsq;
  double FACT;
  double x,y,rmax,ekmax,etaprev;

  struct vector dtau;
  struct basis rec;
  struct intvector K,R;

  if ( ERRLIM > 1.0 ) {
    printf("ERRLIM in find_Ewald_eta is too large!\n");
    exit(0);
  }

  pi   = 4*atan(1.0);
  tpi  = 2*pi;
  tpi2 = tpi*tpi;

  rec = recip(bas);
  taumax = 0.0;
  for(i=0; i<nat; i++){
    for(j=0; j<=i; j++){
      dtau = makevec( at[i].x - at[j].x,
		      at[i].y - at[j].y,
		      at[i].z - at[j].z );
      rsq = vdotprod( dtau , dtau );
      taumax = fmax( taumax , rsq );
    }
  }

  if ( fabs(eta) < 1e-5 ) eta = 1.0;
  FACT = 2.0;
  etaprev = eta;

  for(iter=0; iter<MAXIT; iter++) {
    x = 2.0 * eta * sqrt( -log( ERRLIM ) );
    ekmax = x*x / 2.0;
    K = gbox( ekmax, 1.0, recip(bas) );

    x = sqrt( -log(ERRLIM) ) / eta;
    y = sqrt(taumax);
    rmax = (x+y)*(x+y);
    ekmax = tpi2 * rmax / 2.0;
    R = gbox( ekmax, 1.0, bas );
    
    if ( K.n1*K.n2*K.n3 > 2*R.n1*R.n2*R.n3 ){
      
      if ( iter == MAXIT ) {
	printf("WARNING: possibly suboptimal value of eta in Ewald sums!\n");
	printf("Real space sums over %d %d %d\n",R.n1,R.n2,R.n3);
	printf("K space sums over %d %d %d\n",K.n1,K.n2,K.n3);
      } else {

	if ( eta <= etaprev ) {
	  x = eta / FACT;
	} else {
	  x = eta - (eta-etaprev) / 3.0;
	}
	etaprev = eta;
	eta = x;
      }
	
    } else if ( 2*K.n1*K.n2*K.n3 < R.n1*R.n2*R.n3 ) {
	
	if ( iter == MAXIT ) {
	  printf("WARNING: possibly suboptimal value of eta in Ewald sums!\n");
	  printf("Real space sums over %d %d %d\n",R.n1,R.n2,R.n3);
	  printf("K space sums over %d %d %d\n",K.n1,K.n2,K.n3);
	} else {
	  
	  if ( eta >= etaprev ) {
	    x = eta * FACT;
	  } else {
	    x = eta + (etaprev-eta) / 3.0;
	  }
	  etaprev = eta;
	  eta = x;
	}
    } else {
      if ( debug ) printf("iter=%d Found an optimal value of eta = %f\n",iter,eta);
    }

  } /* end for loop */
	  
  return(eta);
}

double Eion(struct cellprm cell, struct atom *at,int nat,double ERRLIM,double eta)
{
  int i1,i2,i3;
  int iat1,iat2;
  int n,m;

  double pi,tpi,tpi2;
  double rblen,rmax,ekmax;
  double ecc=0;
  double omega;
  double etasq;
  double logerrlim;
  double x,gsq,rsq;
  double taumax,y;
  double temp;

  struct intvector K,R;
  struct basis qb,rb;
  struct vector v1,v2,v3,vg,vr;
  struct vector dtau;

  pi   = 4*atan(1.0);
  tpi  = 2*pi;
  tpi2 = tpi*tpi;

  /* our cell is orthorhombic, so this is easy */
  rblen = cell.a;
  rblen = fmax(rblen, cell.b);
  rblen = fmax(rblen, cell.c);

  rmax = nat * rblen;
  ekmax = tpi2 * rmax / 2.0;

  rb = cell.bas;
  omega = vtriple( rb.A, rb.B, rb.C );
  qb = recip( rb );
  etasq = eta*eta;
  ecc = 0.0;

  logerrlim = log(ERRLIM);

  x = 2.0 * eta * sqrt( -logerrlim );
  ekmax = x * x / 2.0;

  K = gbox(ekmax, 1.0, qb);

  for(i1=-K.n1; i1<(K.n1+1); i1++){
    for(i2=-K.n2; i2<(K.n2+1); i2++){
      for(i3=-K.n3; i3<(K.n3+1); i3++){

	if ( i1==0 && i2==0 && i3==0 ) continue;
	v1 = vsmult( i1, qb.A );
	v2 = vsmult( i2, qb.B );
	v3 = vsmult( i3, qb.C );
	vg = vsmult( tpi, vadd( v1, vadd( v2, v3 ) ) );
	gsq = vdotprod( vg, vg );

	if ( gsq > -4.0 * etasq * logerrlim ) continue;

	for(iat1=0; iat1<nat; iat1++){
	  for(iat2=0; iat2<=iat1; iat2++){
	    dtau = makevec( at[iat1].x-at[iat2].x,
			    at[iat1].y-at[iat2].y,
			    at[iat1].z-at[iat2].z );
	    x = cos( vdotprod( vg, dtau ) ) * exp( -gsq/(4.0*etasq) ) / gsq ;
	    if ( iat1 != iat2 ) x=2.0*x;
	    ecc += at[iat1].chrg * at[iat2].chrg * x;
	  }
	}

      }
    }
  }

  ecc = 4.0 * pi * ecc / omega;

  taumax = 0.0;
  for(m=0; m<nat; m++){
    for(n=0; n<=m; n++){
      dtau = makevec( at[m].x-at[n].x,
		      at[m].y-at[n].y,
		      at[m].z-at[n].z );
      rsq = vdotprod( dtau, dtau );
      taumax = fmax( taumax, rsq );
    }
  }

  x = sqrt( -logerrlim ) / eta;
  y = sqrt( taumax );
  rmax = (x+y)*(x+y);
  ekmax = tpi2 * rmax / 2.0;


  R = gbox(ekmax, 1.0, rb);

  /* fix from VO 9 sep 2005 */
  taumax = sqrt( taumax );

  for(i1=-R.n1; i1<(R.n1+1); i1++){
    for(i2=-R.n2; i2<(R.n2+1); i2++){
      for(i3=-R.n3; i3<(R.n3+1); i3++){
	if ( i1==0 && i2==0 && i3==0 ) continue;
	v1 = vsmult( i1, rb.A );
	v2 = vsmult( i2, rb.B );
	v3 = vsmult( i3, rb.C );
	vr = vadd( v1, vadd( v2, v3 ) );

	/* fix from VO 9 sep 2005 */
	/* if ( ( eta * vdotprod( vr, vr ) ) > -logerrlim ) continue;  WRONG! */
	x = sqrt( vdotprod( vr, vr ) ) - taumax;
	temp= eta * x;
	if( temp*temp > -logerrlim ) continue;
	
	for(iat1=0; iat1<nat; iat1++){
	  for(iat2=0; iat2<=iat1; iat2++){
	    dtau = vadd ( makevec( at[iat1].x-at[iat2].x,
				   at[iat1].y-at[iat2].y,
				   at[iat1].z-at[iat2].z ),
			  vr);
	    rsq = sqrt( vdotprod( dtau, dtau ) );
	    x = erfc( eta * rsq ) / rsq;
	    if ( iat1 != iat2 ) x=2.0*x;

	    /*	    
	    if ( rsq<0.01 ) {
	      printf("i1 i2 i3: %d %d %d\n",i1,i2,i3);
	      printf("at1 at2: %d %d\n",iat1,iat2);
	      printf("ch1 ch2: %f %f\n",at[iat1].chrg,at[iat2].chrg);
	      printf("at1: %e %e %e\n",at[iat1].x,at[iat1].y,at[iat1].z);
	      printf("at2: %e %e %e\n",at[iat2].x,at[iat2].y,at[iat2].z);
	      printf("a1-a2: %e %e %e\n",at[iat1].x-at[iat2].x,at[iat1].y-at[iat2].y,at[iat1].z-at[iat2].z);
	      printf("vr: %e %e %e\n",vr.x,vr.y,vr.z);
 	      printf("dtau: %e %e %e\n",dtau.x,dtau.y,dtau.z);
	      printf("eta=%20.18e\n",eta);
	      printf("rsq=%20.18e\n",rsq);
	      printf("erfc(rsq)=%e\n",erfc(rsq));
	      printf("x=%20.18e\n",x);
	      printf("\n");
	    }
	    */
	    ecc += at[iat1].chrg * at[iat2].chrg * x;
	  }
	}
	
      }
    }
  }

  /*  take care of the term with l=0 in the real-space summation */
  for(iat1=0; iat1<nat; iat1++){
    for(iat2=0; iat2<=(iat1-1); iat2++){
      dtau = makevec( at[iat1].x-at[iat2].x,
		      at[iat1].y-at[iat2].y,
		      at[iat1].z-at[iat2].z );
      rsq = sqrt( vdotprod( dtau, dtau ) );
      x = erfc( eta * rsq ) / rsq;
      ecc += 2.0 * at[iat1].chrg * at[iat2].chrg * x;
    }
  }

  for(iat1=0; iat1<nat; iat1++){
    for(iat2=0; iat2<nat; iat2++){
      ecc -= at[iat1].chrg * at[iat2].chrg * pi / ( eta * eta * omega );
    }
    ecc -= 2.0 * eta * at[iat1].chrg * at[iat1].chrg / sqrt(pi);
  }

  return(ecc);
}

/* c============================================================ */
/* c     Find the max. dimensions of the reciprocal-space */
/* c     box which contains all plane waves with Ekin < Emax. */
/* c */
/* c  DATE: Wed Oct 16 14:54:34 MDT 1996 */
/* c============================================================ */
struct intvector gbox(double emax, double alat, struct basis bas)
{
  double pi,tpi,tpiba,Gmax;
  double q11,q12,q13,q22,q23,q33,rvol,range1,range2,range3;
  struct intvector N;
  
  pi    = 4*atan(1.0);
  tpi   = 2*pi;
  tpiba = tpi / alat;
  Gmax  = sqrt( 2 * emax ) / tpiba;

  /* can generalize this to non-orthorhombic later */
  q11 = vdotprod(bas.A,bas.A);
  q12 = vdotprod(bas.A,bas.B);
  q13 = vdotprod(bas.A,bas.C);
  q22 = vdotprod(bas.B,bas.B);
  q23 = vdotprod(bas.B,bas.C);
  q33 = vdotprod(bas.C,bas.C);

  rvol = vtriple( bas.A, bas.B, bas.C );
  
  range1 = rvol / sqrt(q22*q33 - q23*q23);
  range2 = rvol / sqrt(q11*q33 - q13*q13);
  range3 = rvol / sqrt(q11*q22 - q12*q12);
    
  N.n1 = (int)( Gmax / range1 + 1.0 );
  N.n2 = (int)( Gmax / range2 + 1.0 );
  N.n3 = (int)( Gmax / range3 + 1.0 );

  /*
  printf("q11=%f\n",q11);
  printf("q12=%f\n",q12);
  printf("q13=%f\n",q13);
  printf("q22=%f\n",q22);
  printf("q23=%f\n",q23);
  printf("q33=%f\n",q33);
  printf("gmax=%f\n",Gmax);
  printf("N: %d %d %d\n",N.n1,N.n2,N.n3);
  */

  return(N);
}
