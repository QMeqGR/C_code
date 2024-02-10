/* constants */

#define D2R 0.0174532925199433
#define SRT2 1.4142135623731
#define SRT3 1.73205080756888
#define SRT6 2.44948974278318

/* Structures */

struct cmplx{
  double real;
  double imag;
};

struct vector{
  double x;
  double y;
  double z;
};

struct matrix{
  double r11,r12,r13;
  double r21,r22,r23;
  double r31,r32,r33;
};

struct intvector{
  int n1,n2,n3;
};

struct Energy{
  double ecc,eng,ecc_eng,pf,vol,tot;
};

struct cation{
  int Z;
  double R;
  double chrg;
  struct vector center;
};

struct atom{
  double x;
  double y;
  double z;
  double w; /* atomic weight */
  double R; /* radius of atom */
  double chrg; /* for Ewald sums */
  int Z;
  int typ;
};

/* geometrical object structures */

struct basis{
  struct vector A,B,C;
};

struct face{
  struct vector vert1;
  struct vector vert2;
  struct vector vert3;
  struct vector center;
  struct vector normal;
};

struct cellprm{
  struct basis bas;
  double a,b,c;
  double alph,beta,gamm;
};


struct tetrahedron{
  int Zc,Zv;
  double Rc;
  double Rv;
  struct vector v[5];         /* NOTE: tetr[0] is the center */
  struct face f[4];
  struct vector or;          /* orientation: phi,theta,psi, byron and fuller */
};

struct octahedron{
  int Zc,Zv;
  double Rc;
  double Rv;
  struct vector v[9];         /* NOTE: octa[0] is the center */
  struct face f[8];
  struct vector or;          /* orientation: phi,theta,psi, byron and fuller */
};

struct icosahedron{
  struct vector ico[13];         /* NOTE: ico[0] is the center */
};


/* Function prototypes */

void print_cell_frame_xbs(struct cellprm,double,int);
void print_cat(struct cation);
void print_cat_xbs(char c,struct cation,struct cellprm);
void print_tet(struct tetrahedron);
void print_tet_xbs(struct tetrahedron,struct cellprm);
void print_oct(struct octahedron);
void print_oct_xbs(struct octahedron,struct cellprm);
void print_spe_xbs(char c,double,double,double,double);
void print_bnd_xbs(char c1,char c2,double,double,double);
void check_cats(struct cation *,int);
void check_tetr(struct tetrahedron *,int);
void check_octa(struct octahedron *,int);
void invt_matrx(int,double *);

void simplex_init(void);
void simplex_trans(double *);
void trans_cation(struct cation *,struct vector);
void trans_tetr(struct tetrahedron *,struct vector);
void trans_octa(struct octahedron *,struct vector);
void rescale_tetr(struct tetrahedron *,int,struct cellprm,struct vector);
void rescale_octa(struct octahedron *,int,struct cellprm,struct vector);
void rot_tetr(struct tetrahedron *,int,double,struct cellprm,struct matrix);
void rot_octa(struct octahedron *,int,double,struct cellprm,struct matrix);
void amoeba_ehm(double **,double y[],int,double,struct Energy (*funk)(void), int *);
void NxNmult(int N,double *,double *,double *);

int rsign(void);

double dmod(double,double);
double det3x3(struct matrix);
double dist(struct vector *,struct vector *);
double dist_pbc(struct vector *,struct vector *,struct cellprm *);
double distsq_pbc(struct vector *,struct vector *,struct cellprm *);
double vdotprod(struct vector, struct vector);
double vtriple(struct vector, struct vector, struct vector);
double vmag(struct vector);
double vmagsq(struct vector);
double Eion(struct cellprm, struct atom *,int,double,double);
double find_Ewald_eta(struct basis,double,int,struct atom *);
double amotry_ehm(double **,double y[],double psum[],int,struct Energy (*funk)(),int,double);

struct vector makevec(double,double,double);
struct vector normvec(struct vector);
struct vector vsub(struct vector, struct vector);
struct vector vadd(struct vector, struct vector);
struct vector crossprod(struct vector, struct vector);
struct vector vsmult(double, struct vector);
struct vector multiply(struct matrix *,struct vector *);
struct vector rezone(struct vector);

struct intvector gbox(double, double, struct basis);

struct matrix R_ptp(double,double,double);
struct matrix R_ptp_inv(struct matrix Rf);
struct matrix matrix_for(struct vector,struct vector);

struct face facegen(struct vector *, struct vector *, struct vector *);
struct tetrahedron make_tetr(struct vector,double,double,double,double d,struct cellprm);
struct octahedron  make_octa(struct vector,double,double,double,double d,struct cellprm);

/* function prototypes that don't seem to be in math.h */
double fmax(double,double);
