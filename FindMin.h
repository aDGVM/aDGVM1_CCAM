#pragma once


// Calculate the minimum of a function, needed in Leaf.h

#define GOLD 1.618034 
#define GLIMIT 100.0 
#define TINY 1.0e-20 
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define R_FINDMIN 0.61803399
#define C_FINDMIN (1.0-R_FINDMIN) 
#define SHFT2(a,b,c) (a)=(b);(b)=(c); 
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
// static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))



void mnbrak(double (*func)(double, double *), double *aa, double *ax, double *bx, double *cx, double *fa, double *fb, double *fc );

double golden(double (*f)(double, double *), double *aa, double ax, double bx, double cx, double tol, double *xmin);


// Here GOLD is the default ratio by which successive intervals are magnified; 
// GLIMIT is the maximum magnification allowed for a parabolic-fit step.




// Given a function func, and given distinct initial points ax and bx, this routine 
// searches in the downhill direction (defined by the function as evaluated at the 
// initial points) and returns new points ax, bx, cx that bracket a minimum of the
// function. Also returned are the function values at the three points, fa, fb, and fc.

void mnbrak(double (*func)(double, double *), double *aa, double *ax, double *bx, double *cx, double *fa, double *fb, double *fc ) 
{ 
	double ulim,u,r,q,fu,dum; 
	
	*fa=(*func)(*ax,aa); 
	*fb=(*func)(*bx,aa); 
	
	// Switch roles of a and b so that we can go downhill in the direction from a to b.
	if (*fb > *fa) 
	{  
		SHFT(dum,*ax,*bx,dum) 
		SHFT(dum,*fb,*fa,dum) 
	} 
	
	// First guess for c. 
	*cx=(*bx)+GOLD*(*bx-*ax); 
	*fc=(*func)(*cx,aa); 
	
	// Keep returning here until we bracket.
	while (*fb > *fc) 
	{  
		// Compute u by parabolic extrapolation from a, b, c. 
		// TINY is used to prevent any possible division by zero.
		r=(*bx-*ax)*(*fb-*fc);  
		q=(*bx-*cx)*(*fb-*fa); 
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx); 
		
		// We won t go farther than this. Test various possibilities: 
		if ((*bx-u)*(u-*cx) > 0.0) 
		{ 
			// Parabolic u is between b and c: try it. 
			fu=(*func)(u,aa); 
			if (fu < *fc) 
			{ 
				//Got a minimum between b and c. 
				*ax=(*bx); 
				*bx=u; 
				*fa=(*fb); 
				*fb=fu; 
				return; 
			} 
			else if (fu > *fb) 
			{ 
				//Got a minimum between between a and u. 
				*cx=u; 
				*fc=fu; 
				return; 
			} 
			// Parabolic fit was no use. Use default magnification.
			u=(*cx)+GOLD*(*cx-*bx);  
			fu=(*func)(u,aa); 
		} 
		else if ((*cx-u)*(u-ulim) > 0.0) 
		{ 
			// Parabolic fit is between c and its allowed limit. 
			fu=(*func)(u,aa); 
			if (fu < *fc) 
			{ 
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx)) 
				SHFT(*fb,*fc,fu,(*func)(u,aa)) 
			} 
		} 
		else if ((u-ulim)*(ulim-*cx) >= 0.0) 
		{ 
			// Limit parabolic u to maximum allowed value. 
			u=ulim; 
			fu=(*func)(u,aa); 
		} 
		else 
		{ 
			//Reject parabolic u, use default magnification. 
			u=(*cx)+GOLD*(*cx-*bx); 
			fu=(*func)(u,aa); 
		} 
		// Eliminate oldest point and continue.
		SHFT(*ax,*bx,*cx,u)  SHFT(*fa,*fb,*fc,fu) 
	} 
} 


// Given a function f, and given a bracketing triplet of abscissas ax, bx, cx 
// (such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)), 
// this routine performs a golden section search for the minimum, isolating it to a 
// fractional precision of about tol. The abscissa of the minimum is returned as xmin,
// and the minimum function value is returned as golden, the returned function value.
 
double golden(double (*f)(double, double *), double *aa, double ax, double bx, double cx, double tol, double *xmin)
{ 
	double f1,f2,x0,x1,x2,x3; 
	
	x0=ax; 
	//At any given time we will keep track of four points, x0,x1,x2,x3. 
	x3=cx; 
	
	if (fabs(cx-bx) > fabs(bx-ax)) 
	{ 
		// Make x0 to x1 the smaller segment, 
		x1=bx; 
		x2=bx+C_FINDMIN*(cx-bx); 
		//and fill in the new point to be tried. 
	} 
	else 
	{ 
		x2=bx; 
		x1=bx-C_FINDMIN*(bx-ax); 
	} 
	f1=(*f)(x1,aa); 
	
	// The initial function evaluations. Note that we never need 
	// to evaluate the function at the original endpoints. 
	f2=(*f)(x2,aa);
	while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) 
	{ 
		if (f2 < f1) 
		{ 
			// One possible outcome, its housekeeping, and a new function evaluation.
			SHFT3(x0,x1,x2,R_FINDMIN*x1+C_FINDMIN*x3)  
			SHFT2(f1,f2,(*f)(x2,aa)) 
		} 
		else 
		{ 
			// The other outcome, and its new function evaluation. 
			SHFT3(x3,x2,x1,R_FINDMIN*x2+C_FINDMIN*x0) 
			SHFT2(f2,f1,(*f)(x1,aa)) 
		}
	} 
	// Back to see if we are done. 
	
	if (f1 < f2) 
	{ 
		// We are done. Output the best of the two current values. 
		*xmin=x1; 
		return f1; 
	} 
	else 
	{ 
		*xmin=x2; return f2; 
	}
}




























