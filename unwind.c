
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#define PI (atan(1.0)*4.0)
#define NVAR 4
#define MAXITER 1048576
#define MAXTABLE 8192
#define ASSERT(cond,...)                                                \
  if(!(cond)){fprintf(stderr,"ERROR: ");fprintf(stderr,__VA_ARGS__);exit(2);}


//static int EnableStrings = 1;
static double CompactifiedVolume; // V5, compactified volume to get down to 5d
                                  // theory
static double TenDimensionalPlanckMass; // M10
static double FourDimensionalPlanckMass; // M4
static double StringCoupling;
static double StringMass;
static double LengthOfSmallExtraDimensions;
static double LengthOfLargeExtraDimensions;
static double MuParameter;
static double EpsilonParameter;
static double Lambda0;
static double SizeOfInitialBubble;
static double StartingHubbleConstant;
static double ImpactParameter;
static double KernelSigma; // width of Gaussian collision kernel
static int NumberOfFluxes; // Q0
static int DimensionOfBrane;

static int integer_part(double x);
static double F_function(double chi);
static double potential_function(double z);
static double potential_derivative(double z);
static double kernel_function(double z);
static double S1(double z);
static double S2(double z);
static double numerically_integrate(double *x, double *y, int N);
static void rungekutta4(double *X, double dt);
static void dXdt(double *X, double *dXdt);

static double TimeStep = 1.0e-2;
static double TrajectoryRecord[MAXITER][4]; // store the previous history
static int IterationNumber = 0;
static int DumpCadence = 500;
static int RestartFromDump = 0;
static char DumpFilename[] = "unwind-restart.bin";

static double PolylogTableX[MAXTABLE];
static double PolylogTableY[MAXTABLE];
static int PolylogTableN;
static int LastLookUp;
static int WindowSize;

double numerically_integrate(double *x, double *y, int N)
{
  int n;
  double S = 0.0;
  for (n=1; n<N; ++n) {
    S += (x[n] - x[n-1]) * (y[n] + y[n-1]);
  }
  return 0.5 * S;
}
void rungekutta4(double *X, double dt)
{
  int n;
  double L1[NVAR], X1[NVAR];
  double L2[NVAR], X2[NVAR];
  double L3[NVAR], X3[NVAR];
  double L4[NVAR], X4[NVAR];
  for (n=0; n<NVAR; ++n) X1[n] = X[n] +         0.0 * dt; dXdt(X1, L1);
  for (n=0; n<NVAR; ++n) X2[n] = X[n] + L1[n] * 0.5 * dt; dXdt(X2, L2);
  for (n=0; n<NVAR; ++n) X3[n] = X[n] + L2[n] * 0.5 * dt; dXdt(X3, L3);
  for (n=0; n<NVAR; ++n) X4[n] = X[n] + L3[n] * 1.0 * dt; dXdt(X4, L4);
  for (n=0; n<NVAR; ++n) {
    X[n] += (L1[n] + 2*L2[n] + 2*L3[n] + L4[n]) * dt/6.0;
  }
}

void dXdt(double *X, double *dXdt)
{
  double z = X[1];
  double g = X[2];
  double a = X[3];

  if (g < 1.0) {
    fprintf(stderr, "Stopping the unwind because gamma=%6.5f < 1.0\n", g);
    exit(0);
  }
 
  int p       = DimensionOfBrane;
  double dVdz = potential_derivative(z); 
  double chi  = atanh(sqrt(g*g - 1)/g);
  double gs   = StringCoupling;
  double dp   = pow(LengthOfSmallExtraDimensions, p-3);
  double ms   = StringMass;
  double sig  = pow(ms, p+1) * dp / (gs * 0.5*pow(2*PI, p-1));
  double M4   = FourDimensionalPlanckMass;
  double V    = potential_function(z);

  double S1z = S1(z);
  double S2z = S2(z);

  dXdt[0] = 1.0; // = t'
  dXdt[1] = sqrt(g*g - 1) / g; // = z'
  dXdt[3] = a/(sqrt(3.0)*M4) *
    sqrt(V + sqrt(g*g - 1)/(g*chi) * S1z/pow(a,3) + 2.0*sig*g); // = a'

  double term1 = -3.0 * dXdt[3]/a * (g*g - 1)/g;
  double term2 = -0.5 * dVdz * sqrt(g*g - 1) / (g * sig);
  double term3 = -0.5 * (g*g - 1)/(g*chi) * 1/( sig * pow(a, 3)) * S2z;
  double denom = 1.0 + 0.5/(g*chi) * (1/(g*sqrt(g*g-1)) - 1/chi) *
    S1z/(sig * pow(a, 3));

  dXdt[2] = (term1 + term2 + term3) / denom; // = g'
}

void initialize_constants()
/*
 *  These are initialized when reading into the PolylogTable:
 *
 *  DimensionOfBrane
 *  ImpactParameter // b
 *
 *  Make sure that load_polylog_table is being called before this!
 *
 */
{
  int Q0     = 400;
  int p      = DimensionOfBrane;
  double lm0 = 0.0;
  double ms  = 1.0;
  double L   = 19.7 / ms;
  double d   =  2.0 / ms;
  double dp  = pow(d, p-3); // also called Vwrap in da-code
  double V5  = pow(d, 5);
  double gs  = 1e-2;
  double M10 = ms / pow(gs, 0.25);
  double M4  = sqrt(V5 * L) * pow(M10, 4);
  double mu  = pow(pow(pow(ms, 4) / (gs * 2 * PI), 2) *
                   pow(M10, 2*p-14) * pow(dp, 2) / V5, 0.2);
  double ep  = pow(mu, 5) * (Q0 - 0.5);
  double sig = pow(ms, p+1) * dp / (gs * pow(2*PI, 0.5*(p-1)));

  MuParameter = mu;
  EpsilonParameter = ep;
  CompactifiedVolume = pow(d, 5);
  MuParameter = mu;
  SizeOfInitialBubble =  4 * sig / ep; // r0
  StartingHubbleConstant = sqrt((pow(mu, 5)*Q0*Q0*L/2 + lm0) / (3*M4*M4));
  TenDimensionalPlanckMass = M10;
  FourDimensionalPlanckMass = M4;
  NumberOfFluxes = Q0;
  StringCoupling = gs;
  StringMass = ms;
  Lambda0 = lm0; // background cosmological constant
  CompactifiedVolume = pow(d, 5);
  KernelSigma = 1.0;
  LengthOfSmallExtraDimensions = d;
  LengthOfLargeExtraDimensions = L;
}

void load_polylog_table()
{
  char line[1024];
  int n=0;
  FILE *table = fopen("PolyLog_TABLES/ftempLL-b0over5-p4.dat", "r");

  fgets(line, 1024, table); ImpactParameter = atof(line);
  fgets(line, 1024, table); DimensionOfBrane = atoi(line);

  /* WARNING: this code only works when the ASCII table does *not* end in a new
     line */
  while (!feof(table)) {
    fscanf(table, "%lf %lf\n", &PolylogTableX[n], &PolylogTableY[n]);
    ++n;
    ASSERT(n < MAXTABLE, "polylog table exceeded max table size %d", MAXTABLE);
  }
  PolylogTableN = n;
  fclose(table);
}

int main()
{
  LastLookUp = PolylogTableN+1; // attempting to speed up polylog lookup

  load_polylog_table();
  initialize_constants();

  /* The initial conditions are either created by reading from a file */
  if (RestartFromDump) {
    FILE *restart = fopen(DumpFilename, "r");
    fread(&IterationNumber, sizeof(int), 1, restart);
    fread(&TrajectoryRecord[0][0], sizeof(double), MAXITER*NVAR, restart);
    fclose(restart);
    //    IterationNumber = 38000;
  }
  /* or by populating the first entry of TrajectoryRecord with the initial
     values of X */
  else {
    IterationNumber = 1;
    TrajectoryRecord[0][0] = 0.0; // t
    TrajectoryRecord[0][1] = SizeOfInitialBubble; // z
    TrajectoryRecord[0][2] = 1.0001; // g
    TrajectoryRecord[0][3] = 1.0; // a **** I THINK THIS IS ARBITRARY ******
  }

  int i;
  double dt = TimeStep;
  double X[4];

  for (i=0; i<NVAR; ++i) {
    X[i] = TrajectoryRecord[IterationNumber-1][i];
  }

  while (X[0] < 4000.0) {

    double Q = NumberOfFluxes - 2.0*X[1]/LengthOfLargeExtraDimensions;
    printf("%d %14.12e %14.12e %14.12e %14.12e %14.12e\n", IterationNumber,
	   X[0], X[1], X[2], X[3], Q);

    ASSERT(IterationNumber < MAXITER, "integration exceeded maximum size %d",
           MAXITER);

    rungekutta4(X, dt);
    ++IterationNumber;


    for (i=0; i<NVAR; ++i) {
      TrajectoryRecord[IterationNumber-1][i] = X[i];
    }

    if (IterationNumber % DumpCadence == 0) {
      FILE *restart = fopen(DumpFilename, "w");
      fwrite(&IterationNumber, sizeof(int), 1, restart);
      fwrite(&TrajectoryRecord[0][0], sizeof(double), MAXITER*NVAR, restart);
      fclose(restart);
    }

    double z = X[1];
    double g = X[2];
    double a = X[3];

    int p       = DimensionOfBrane;
    double gs   = StringCoupling;
    double dp   = pow(LengthOfSmallExtraDimensions, p-3);
    double ms   = StringMass;
    double sig  = pow(ms, p+1) * dp / (gs * 0.5*pow(2*PI, p-1));
    double chi  = atanh(sqrt(g*g-1)/g);

    double rhos = sqrt(g*g - 1)/(pow(a,3) * g * chi)* S1(z);
    double kin  = 2.0*g*sig;

    if (rhos > kin) {
      fprintf(stderr, "BRANE REHEATS: energy density in strings "
	      "is greater than kinetic energy\n");
      printf("%14.8e %14.8e %14.8e %14.8e %14.8e %14.8e %14.8e \n",
	     X[0], X[1], X[2], X[3], Q, rhos, kin);
      return 0;
    }
  }

  return 0;
}








/*
 *
 * Physics functions
 *
 */
double Sx(double z, int which)
{
  WindowSize = 4.0*LengthOfLargeExtraDimensions/2.0 /TimeStep;
  
  int i, i0, i1;
  
  i1 = IterationNumber - 1;
  i0 = i1 - WindowSize;
  
  if (i0 < 0) i0 = 0;

  double *integrand = (double*) malloc((i1-i0) * sizeof(double));
  double *times     = (double*) malloc((i1-i0) * sizeof(double));

  for (i=i0; i<i1; ++i) {

    double *Xi = TrajectoryRecord[i];
    double ti = Xi[0];
    double zi = Xi[1];
    double gi = Xi[2];
    double ai = Xi[3];

    double K = 0.0; // sum of kernel over collision points
    double zcol = 0.0; // z-coordinate of collision points

    double L = LengthOfLargeExtraDimensions;
    
    while (zcol < z) {
      zcol += 0.5 * L;
      K += kernel_function(zcol - zi);
    }

    double ms   = StringMass;
    int p       = DimensionOfBrane;
    double dp   = pow(LengthOfSmallExtraDimensions, p-3); // also called Vwrap in da-code
    double chi  = atanh(sqrt(gi*gi-1)/gi);
    
    if (which == 1) {
      integrand[i-i0] = dp * pow(ms,p+2)/pow(2*PI, p) *
	pow(chi,0.5*p)* F_function(chi) * abs(z-zi) * pow(ai,3) * K;
    }
    else if (which == 2) {
      integrand[i-i0] = dp * pow(ms,p+2)/pow(2*PI, p) *
	pow(chi,0.5*p)* F_function(chi) * abs(z-zi)/(z-zi) * pow(ai,3) * K;
    }
    times[i-i0] = ti;
  }

  double definite_integral =
    numerically_integrate(times, integrand, i1-i0);
  free(integrand);
  free(times);

  return definite_integral;  
}

double S1(double z)
{
  return Sx(z, 1);
}
double S2(double z)
{
  return Sx(z, 2);
}

double potential_function(double z)
{
  int Q0     = NumberOfFluxes;
  double mu  = MuParameter;
  double L   = LengthOfLargeExtraDimensions;
  double V0  = pow(mu, 5) * Q0*Q0 * L * 0.5;
  return V0 - 2 * pow(mu, 5) *
    ((Q0 - 0.5)*z - integer_part(2*z/L) * (z - L*(integer_part(2*z/L) + 1)/4.0));
 
}
double potential_derivative(double z)
{
  int Q0     = NumberOfFluxes;
  double mu  = MuParameter;
  double L   = LengthOfLargeExtraDimensions;
  return -2 * pow(mu, 5) * (Q0 - 0.5 - integer_part(2*z/L));
}
double kernel_function(double z)
{
  double s2 = KernelSigma * KernelSigma;
  return pow(2*PI*s2, -0.5) * exp(-0.5 * z*z / s2);
}





/*
 *
 * Helper functions
 *
 */

double F_function(double chi)
{
  double x = log10(chi);
  int i;
  if (LastLookUp == PolylogTableN+1){
    for (i=0; i<PolylogTableN-1; ++i) {
      if (PolylogTableX[i] < x && x <= PolylogTableX[i+1]) {
	double dy = PolylogTableY[i+1] - PolylogTableY[i];
	double dx = PolylogTableX[i+1] - PolylogTableX[i];
	double y  = PolylogTableY[i] + (x - PolylogTableX[i]) * dy/dx;
	LastLookUp = i;
	return pow(10.0, y);
      }
    }
  }
  else {
    if(PolylogTableX[LastLookUp] > x){
      for (i=LastLookUp; i>0; --i){
	if (PolylogTableX[i] > x && x >= PolylogTableX[i-1]) {
	  double dy = PolylogTableY[i] - PolylogTableY[i-1];
	  double dx = PolylogTableX[i] - PolylogTableX[i-1];
	  double y  = PolylogTableY[i-1] + (x - PolylogTableX[i-1]) * dy/dx;
	  LastLookUp = i;
	  return pow(10.0, y);
	}
      }
    }
    else{
      for (i=LastLookUp; i<PolylogTableN-1; ++i) {
	if (PolylogTableX[i] < x && x <= PolylogTableX[i+1]) {
	  double dy = PolylogTableY[i+1] - PolylogTableY[i];
	  double dx = PolylogTableX[i+1] - PolylogTableX[i];
	  double y  = PolylogTableY[i] + (x - PolylogTableX[i]) * dy/dx;
	  LastLookUp = i;
	  return pow(10.0, y);
	}
      }
    }
  }


  ASSERT(0, "polylog value lookup failed on the value chi=%f\n", chi);
}
int integer_part(double x)
{
  if (x > 0) {
    return floor(x);
  } else {
    return ceil(x);
  }
}
