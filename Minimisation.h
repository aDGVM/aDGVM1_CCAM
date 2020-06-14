#pragma once

#include "MyNrutil.h"

#define ITMAX 500
#define TOL 2.0e-8

// Minimization of a function func of n variables. Input consists of an initial starting point
// p[1..n]; an initial matrix xi[1..n][1..n], whose columns contain the initial set of directions
// (usually the n unit vectors); and ftol, the fractional tolerance in the function value
// such that failure to decrease by more than this amount on one iteration signals doneness. On
// output, p is set to the best point found, xi is the then-current direction set, fret is the returned
// function value at p, and iter is the number of iterations taken. The routine linmin is used.

void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret, double (*func)(double[]))
{
	void linmin(double p[], double xi[], int n, double *fret, double (*func)(double[]));
	int i, ibig, j;
	double del, fp, fptt, t, *pt, *ptt, *xit;

	pt = vector(1, n);
	ptt = vector(1, n);
	xit = vector(1, n);
	*fret = (*func)(p);

	for (j = 1; j <= n; j++)
		pt[j] = p[j];
	for (*iter = 1;; ++(*iter))
	{
		fp = (*fret);
		ibig = 0;
		del = 0.0;

		for (i = 1; i <= n; i++)
		{
			for (j = 1; j <= n; j++)
			{
				xit[j] = xi[j][i];

			}

			fptt = (*fret);
			linmin(p, xit, n, fret, func); //minimize along it,
			if (fptt - (*fret) > del)
			{
				del = fptt - (*fret);
				ibig = i;
			}
		}
		if (2.0 * (fp - (*fret)) <= ftol * (fabs(fp) + fabs(*fret)) + TINY)
		{
			free_vector(xit, 1, n);
			free_vector(ptt, 1, n);
			free_vector(pt, 1, n);
			return;
		}
		if (*iter == ITMAX) nrerror("powell exceeding maximum iterations.");
		for (j = 1; j <= n; j++)
		{
			ptt[j] = 2.0 * p[j] - pt[j];
			xit[j] = p[j] - pt[j];
			pt[j] = p[j];
		}
		fptt = (*func)(ptt);
		if (fptt < fp)
		{
			t = 2.0 * (fp - 2.0 * (*fret) + fptt) * SQR(fp - (*fret) - del) - del * SQR(fp - fptt);
			if (t < 0.0)
			{
				linmin(p, xit, n, fret, func);
				for (j = 1; j <= n; j++)
				{
					xi[j][ibig] = xi[j][n];
					xi[j][n] = xit[j];
				}
			}
		}
	}
}

int ncom;

double *pcom, *xicom, (*nrfunc)(double[]);

void linmin(double p[], double xi[], int n, double *fret, double (*func)(double[]))
// Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and
// resets p to where the function func(p) takes on a minimum along the direction xi from p,
// and replaces xi by the actual vector displacement that p was moved. Also returns as fret
// the value of func at the returned location p. This is actually all accomplished by calling the
// routines mnbrak and brent.
{
	double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
	double f1dim(double x);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double));
	int j;
	double xx, xmin, fx, fb, fa, bx, ax;
	ncom = n;
	pcom = vector(1, n);
	xicom = vector(1, n);
	nrfunc = func;
	for (j = 1; j <= n; j++)
	{
		pcom[j] = p[j];
		xicom[j] = xi[j];
	}
	ax = 0.0;
	xx = 1.0;
	mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
	*fret = brent(ax, xx, bx, f1dim, TOL, &xmin);
	for (j = 1; j <= n; j++)
	{
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom, 1, n);
	free_vector(pcom, 1, n);
}

extern int ncom;
extern double *pcom, *xicom, (*nrfunc)(double[]);
double f1dim(double x)
{
	int j;
	double f, *xt;
	xt = vector(1, ncom);
	for (j = 1; j <= ncom; j++)
		xt[j] = pcom[j] + x * xicom[j];
	f = (*nrfunc)(xt);
	free_vector(xt, 1, ncom);
	return f;
}

// #define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin)
// Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
// between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
// the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
// the minimum is returned as xmin, and the minimum function value is returned as brent, the
// returned function value.
{
	int iter;
	double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	double e = 0.0;

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = w = v = bx;
	fw = fv = fx = (*f)(x);
	for (iter = 1; iter <= ITMAX; iter++)
	{
		xm = 0.5 * (a + b);
		tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
		if (fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
		{
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1)
		{
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0) p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) d = CGOLD * (e = (x >= xm ? a - x : b - x));
			else
			{
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2) d = SIGN(tol1, xm - x);
			}
		}
		else
		{
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		}
		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu = (*f)(u);
		if (fu <= fx)
		{
			if (u >= x) a = x;
			else
				b = x;
			SHFT(v, w, x, u)
			SHFT(fv, fw, fx, fu)
		}
		else
		{
			if (u < x) a = u;
			else
				b = u;
			if (fu <= fw || w == x)
			{
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || v == x || v == w)
			{
				v = u;
				fv = fu;
			}
		}
	}
	nrerror("Too many iterations in brent");
	*xmin = x;
	return fx;
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double))
{
	double ulim, u, r, q, fu, dum;
	*fa = (*func)(*ax);
	*fb = (*func)(*bx);
	if (*fb > *fa)
	{
		SHFT(dum, *ax, *bx, dum)
		SHFT(dum, *fb, *fa, dum)
	}
	*cx = (*bx) + GOLD * (*bx - *ax);
	*fc = (*func)(*cx);
	while (*fb > *fc)
	{
		r = (*bx - *ax) * (*fb - *fc);
		q = (*bx - *cx) * (*fb - *fa);
		u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) / (2.0 * SIGN(FMAX(fabs(q-r),TINY), q - r));
		ulim = (*bx) + GLIMIT * (*cx - *bx);
		if ((*bx - u) * (u - *cx) > 0.0)
		{
			fu = (*func)(u);
			if (fu < *fc)
			{
				*ax = (*bx);
				*bx = u;
				*fa = (*fb);
				*fb = fu;
				return;
			}
			else if (fu > *fb)
			{
				*cx = u;
				*fc = fu;
				return;
			}
			u = (*cx) + GOLD * (*cx - *bx);
			fu = (*func)(u);
		}
		else if ((*cx - u) * (u - ulim) > 0.0)
		{
			fu = (*func)(u);
			if (fu < *fc)
			{
				SHFT(*bx, *cx, u, *cx+GOLD*(*cx-*bx))
				SHFT(*fb, *fc, fu, (*func)(u))
			}
		}
		else if ((u - ulim) * (ulim - *cx) >= 0.0)
		{
			u = ulim;
			fu = (*func)(u);
		}
		else
		{
			u = (*cx) + GOLD * (*cx - *bx);
			fu = (*func)(u);
		}
		SHFT(*ax, *bx, *cx, u)
		SHFT(*fa, *fb, *fc, fu)
	}
}

