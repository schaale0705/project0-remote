#define VERSION_RT_AKIMA "0108150000"
#define verbosity 1

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <malloc.h>

int main(int, char **);

int    lobatto_coeff(int nx, double *x, double *c);
double firstguess(int n, double term);
double newton(int n, double x);
int    legendre(int nn, double x, double *pl, double *pm, double *pn);
double besselzero(int s);

int    phase_truncation(double *x, double *f, int n, 
			double critratio, double thetatrunc);
int    akima_interpolation(double *x, double *f, int n, 
			   double *x2, double *f2, int n2, int phase_flag);
int    akima_interpolation0(int n, double *x, double *f, double *a, 
			    double *b, double *c, double *d, int phase_flag);

/*-------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
    FILE   *ff;
    int    val=0, i, n, n2, extrapol_type=0;
    double *x, *f, *x2, *f2, dx, pi;

    pi = acos(-1.0);

    n = 2;
    x  = (double *) calloc(n ,sizeof(double));
    f  = (double *) calloc(n ,sizeof(double));

    lobatto_coeff(n, x, f);
    for(i=0;i<n;i++) fprintf(stdout,"-- Lobatto %3d -- : %+e %+e %+e\n",i,x[i],acos(x[i])*180.0/pi,f[i]);
    
    free((void *) x);
    free((void *) f);
    exit(0);

    // Try to open a phase function file
    if((ff=fopen("rt_phase.dat","r")))
	{
	// The original function f(x) @ x
	n=171;
	x  = (double *) calloc(n ,sizeof(double));
	f  = (double *) calloc(n ,sizeof(double));

	// The interpolated values f2(x2) @ x2
	n2=180*2*100+1;
	x2 = (double *) calloc(n2,sizeof(double));
	f2 = (double *) calloc(n2,sizeof(double));

	// Read a phase function file
	ff=fopen("rt_phase.dat","r");
        for(i=0;i<n;i++) fscanf(ff,"%*g %lg %lg %*g %*g %*g\n",&x[i],&f[i]);
	fclose(ff);

	// Get an array with x values where f[] shall be interpolated
	dx = 180.0 / (double) (n2-1);
	for(i=0;i<n2;i++) x2[i]=dx*(double) i;

	//    phase_truncation(x, f, n, 0.2, 10.0);

	// print original values
	ff = fopen("rt_akima_ori","w");
	fprintf(ff,"# ex_norm\n");
	for(i=0;i<n;i++) fprintf(ff,"%3d %+e %+e\n",i,x[i],f[i]);
	fclose(ff);
	system("phase_plot rt_akima_ori rt_akima_ori_xy");

	if(argc > 1) extrapol_type = atoi(argv[1]);
	// Interpolate now
	//-------------------------------------------|
	//        = 0 : Akima standard extrapolation |
	//        = 1 : Use phase function symmetry  |
	val=akima_interpolation(x, f, n, x2, f2, n2, extrapol_type);
	fprintf(stdout,"--A--> return .............. = %d\n",val);

	// print interpolated values
	ff = fopen("rt_akima_int","w");
	fprintf(ff,"# ex_norm\n");
	for(i=0;i<n2;i++) fprintf(ff,"%3d %+e %+e\n",i,x2[i],f2[i]);
	fclose(ff);
	system("phase_plot rt_akima_int rt_akima_int_xy");
	}
    else
	{
	// The original function f(x) @ x
	n=100;
	x  = (double *) calloc(n ,sizeof(double));
	f  = (double *) calloc(n ,sizeof(double));

	// The interpolated values f2(x2) @ x2
	n2=300;
	x2 = (double *) calloc(n2,sizeof(double));
	f2 = (double *) calloc(n2,sizeof(double));

	// Fill in some function values

        // Parabola
	for(i=0;i<n;i++)
	    {
	    x[i] = (double) i;
	    f[i] = 2.0 - 5.0*x[i] + 0.1*x[i]*x[i];
	    }

	// Dirac needles
	for(i=0;i<n;i++)
	    {
	    x[i] = (double) i;
	    f[i] = 0.0;
	    }

	f[n/2]   = 1.0;
	f[n/2+2] = 1.0;

	// Steps
	for(i=0;i<n;i++)
	    {
	    x[i] = (double) i;
	    f[i] = (double) (i/10);
	    }

	// Get an array with x values where f[] shall be interpolated
	for(i=0;i<n2;i++) x2[i]=10.0+0.2*(double) i;

        // print original values
	ff = fopen("rt_akima_ori","w");
	for(i=0;i<n;i++) fprintf(ff,"%3d %+e %+e\n",i,x[i],f[i]);
	fclose(ff);

	// Interpolate now
	//-------------------------------------------|
	//        = 0 : Akima standard extrapolation |
	//        = 1 : Use phase function symmetry  |
	val=akima_interpolation(x, f, n, x2, f2, n2, 0);
	fprintf(stdout,"--A--> return .............. = %d\n",val);

	// print interpolated values
	ff = fopen("rt_akima_int","w");
	for(i=0;i<n2;i++) fprintf(ff,"%3d %+e %+e\n",i,x2[i],f2[i]);
	fclose(ff);
	}

    // Free allocated memory
    free((void *) x);
    free((void *) f);
    free((void *) x2);
    free((void *) f2);

    return(val);
}

/*-------------------------------------------------------------------------*/

int lobatto_coeff(int nx, double *x, double *c)
{
//  Find abcissa x[] and corresponding weights for a
//  Gauss-Lobatto quadrature scheme
    int    n, term, zero, zeros;
    // double pi, lx, xl;
    double lx, xl;
    double pl, pm, pn;
    
    //pi=acos(-1.0);
    n = 2*nx;
    
    x[0] = 1.0;
    c[0] = 2.0/((double) (n*(n-1)));
    zeros = nx-1;

    // Find the roots (zeros) by Newton's iteration scheme
    for(term=1,zero=0;zero<zeros;zero++,term++)
	{
	lx=19262.0;
	xl=firstguess(n,zero);
	while(fabs(xl-lx) > 1.0e-15)
	    {
	    lx = xl;
	    xl -= newton(n, xl);
	    }
	// OK, got it ...
	x[term] = xl;
	// ... now get the corresponding weight
	legendre(n,xl,&pl,&pm,&pn);
	c[term] = 2.0 / ( (double) (n*(n-1))*pm*pm );
	}

    return(0);
}

/*-------------------------------------------------------------------------*/

double firstguess(int n, double term)
{
    double val, pi;
    pi=acos(-1.0);
    val = cos( besselzero(term) /
	  sqrt(  ((double) n - 0.5 ) * ((double) n - 0.5 ) 
		+ (0.25 - 1.0/pi/pi) ));
    return(val);
}

/*-------------------------------------------------------------------------*/

double newton(int n, double x)
{
    double val, pl, pm, pn;

    // Get the Legendre values
    legendre(n-1, x, &pl, &pm, &pn);

    // The MOMO version of Pn''(x) / Pn'(x) ...
    val = (1.0-x*x)*(pm-x*pn) / 
	  ( n*pl + 2.0*(1.0-n)*x*pm + (x*x*(n-1.0) -1.0)*pn );

    // ... and my version ....
    val = (1.0-x*x)*(pm-x*pn) / 
	  ( (n-1.0)*pl + 2.0*(1.5-n)*x*pm + (x*x*(n-1.0) -1.0)*pn );

    return(val);
}

/*-------------------------------------------------------------------------*/

int legendre(int nn, double x, double *pl, double *pm, double *pn)
{
    // nn : order of Legendre polynomial
    // x  : argument
    // pl, pm, pn : Legendre poly values for n-2, n-1 and n

    int n;
    double in;
    
    *pl=1.0;
    *pm=1.0;
    *pn = x;
    
    // To derive this simply use the Legendre recursion
    // formula and substitute 'k' by 'n-1'
    for(n=2;n<=nn;n++)
	{
	in = 1.0/(double) n;
	*pl = *pm;
	*pm = *pn;
	*pn = (2.0-in) *x * *pm - (1.0-in) * *pl;
	}
    return(0);
}

/*-------------------------------------------------------------------------*/

double besselzero(int s)
{
    double pi, val, b;
    double a1,a2,a3,a4,a5;
    
    a1 =    -0.607927101854027;
    a2 =     6.159589352810601e-02;
    a3 =    -0.981080301612648;
    a4 =    11.75078175473267;
    a5 = -2731.249782035937;
    
    pi=acos(-1.0);


    b = 4*(double)(s+1)+1;    

    val = 0.25 * pi * b 
	+ 0.25 * pi /b * (a1 + (a2 + (a3 + (a4 + a5/b/b)/b/b)/b/b)/b/b);

    return(val);
}

/*-------------------------------------------------------------------------*/

int phase_truncation(double *x, double *f, int n, double critratio, double thetatrunc)
{
    int i, itrunc;
    double m, a, c;

    // Find the cut-off angle
    i=0;
    while(x[i] < thetatrunc && i < n && f[0]/f[i] < 1.0/critratio) i++;
    itrunc=i;

    fprintf(stdout,"-----> itrunc = %d, ratio = %e\n",itrunc,f[0]/f[itrunc]);

    if(f[0]/f[itrunc] >= critratio && i < n)
	{
	// Truncate it !!
	// Use a parabola approach
	m = (f[itrunc+1] - f[itrunc])/(x[itrunc+1] - x[itrunc]);
	a = m / 2.0 / x[itrunc];
	c = f[itrunc] - m*x[itrunc]/2.0;

	for(i=0;i<itrunc;i++) f[i] = a * x[i] * x[i] + c;
	return(0);
	}
    
    return(1);
}

/*-------------------------------------------------------------------------*/

int akima_interpolation(double *x,  double *f,  int n,   // Original
			double *x2, double *f2, int n2,	 // Interpolated
			int phase_flag)
{

//    Routines from formulas as given in 
//    c't magazin fuer computertechnik, 6(1989)206-214

    int    i, j, ishow=0;
    double *a, *b, *c, *d, xtmp;

    // Allocate arrays for akima coefficients
    a  = (double *) calloc((n+4),sizeof(double));
    b  = (double *) calloc((n+4),sizeof(double));
    c  = (double *) calloc((n+4),sizeof(double));
    d  = (double *) calloc((n+4),sizeof(double));

    // Extrapolate the values left and right of the data space
    // Get the line slopes and calculate the coefficients (a,b,c,d)
    akima_interpolation0(n, x, f, a, b, c, d, phase_flag);

    // Get the interpolated functions values f2 @ x2
    for(i=0;i<n2;i++)
	{
	ishow=1;
	for(j=0;j<n-1;j++)
	    {
	    if((x2[i] >= x[j]) && (x2[i] <= x[j+1]))
		{
		xtmp = x2[i] - x[j];
		f2[i] = a[j] + xtmp * ( b[j] + xtmp * (c[j] + d[j]*xtmp));
		ishow = 0;
		}
	    }
	if(ishow) fprintf(stdout,"ALARM : x2[%d] is out of range ....\n",i);
	}

    // Free allocated memory
    free((void *) a);
    free((void *) b);
    free((void *) c);
    free((void *) d);

    return(ishow);
}

/*-------------------------------------------------------------------------*/

int akima_interpolation0(int n, double *x, double *f, double *a, 
			 double *b, double *c, double *d, int phase_flag)
{
    int    i, val=0;
    double *x1, *f1, *st, *t;

    x1 = (double *) calloc((n+4),sizeof(double));
    f1 = (double *) calloc((n+4),sizeof(double));
    st = (double *) calloc((n+4),sizeof(double));
    t  = (double *) calloc((n+4),sizeof(double));

    // Remap the input data (shift index up twice)
    for(i=0;i<n;i++)
	{
	x1[i+2] = x[i];
	f1[i+2] = f[i];
	}

    if(phase_flag)
	{
	// ----------------------------------------------------
	// Get 'extrapolated' values from known
	// properties of a phase function. No need to guess !
	// ----------------------------------------------------
	
	// Guess x1 and f1(x1) values left ...
	x1[1]   = x1[2] - (x1[3]-x1[2]);
	f1[1]   = f1[3];
	x1[0]   = x1[2] - (x1[4]-x1[2]);
	f1[0]   = f1[4];

	// ... and right of data space
	x1[n+2] = x1[n+1] + (x1[n+1] - x1[n]);
	f1[n+2] = f1[n];
	x1[n+3] = x1[n+1] + (x1[n+1] - x1[n-1]);
	f1[n+3] = f1[n-1];

	// Calculate slopes for known points
	for(i=0;i<n+3;i++) st[i] = (f1[i+1] - f1[i])/(x1[i+1]-x1[i]);

	// for(i=0;i<n+4;i++) fprintf(stdout,"%3d %+e %+e %+e\n",i,x1[i],f1[i],st[i]);
	}
    else
	{
	// Guess x values left ...
	x1[1]   = x1[2]   + x1[3]   - x1[4];
	x1[0]   = x1[1]   + x1[2]   - x1[3];

        // ... and right of data space
	x1[n+2] = x1[n+1] + x1[n]   - x1[n-1];
	x1[n+3] = x1[n+2] + x1[n+1] - x1[n];

	// Calculate slopes for known points
	for(i=2;i<n+1;i++) st[i] = (f1[i+1] - f1[i])/(x1[i+1]-x1[i]);

	// Get extrapolated function values and slopes left ...
	f1[1]   = (x1[2] - x1[1]) * (st[3] - 2.0 * st[2]) + f1[2];
	st[1]   = (f1[2] - f1[1]) / (x1[2] - x1[1]);

	f1[0]   = (x1[1] - x1[0]) * (st[2] - 2.0 * st[1]) + f1[1];
	st[0]   = (f1[1] - f1[0]) / (x1[1] - x1[0]);

	// ... and right of data space
	f1[n+2] = (2.0 * st[n]   - st[n-1]) * (x1[n+2] - x1[n+1]) + f1[n+1];
	st[n+1] = (f1[n+2] - f1[n+1]) / (x1[n+2] - x1[n+1]);

	f1[n+3] = (2.0 * st[n+1] - st[n])   * (x1[n+3] - x1[n+2]) + f1[n+2];
	st[n+2] = (f1[n+3] - f1[n+2]) / (x1[n+3] - x1[n+2]);
	}

    // Get polynom slopes
    for(i=2;i<n+2;i++)
	{
	if(st[i-2] == st[i-1] && st[i] == st[i+1])
	    t[i] = 0.5 * (st[i+1] +st[i]);
	else
	    t[i] = (fabs (st[i+1] - st[i]) * st[i-1] + fabs (st[i-1] - st[i-2]) * st[i]) / 
		   (fabs (st[i+1] - st[i]) + fabs (st[i-1] - st[i-2]));
	}

    // Shift index back
    for(i=0;i<n;i++)
	{
	x1[i] = x1[i+2];
	f1[i] = f1[i+2];
	t[i]  = t[i+2];
	}

    // Calculate polynom coefficients
//    for(i=0;i<n-2;i++)
    for(i=0;i<n;i++)
	{
	a[i] = f1[i];
	b[i] = t[i];
        c[i] = (3.0 * st[i+2] - 2.0 * t[i] - t[i+1]) / (x1[i+1] - x1[i]);
	d[i] = (t[i] + t[i+1] - 2.0 * st[i+2]) / ((x1[i+1] - x1[i]) * (x1[i+1] - x[i]));
//	fprintf(stdout,"%4d %+e %+e %+e %+e\n",i,a[i],b[i],c[i],d[i]);
	}

    // Free buffer space
    free((void *) x1);
    free((void *) f1);
    free((void *) st);
    free((void *) t);

    return(val);
}

/*-------------------------------------------------------------------------*/
