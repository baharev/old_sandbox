/*
Program for calculating minimal detectable differences for general ANOVA models
Copyright (C) 2008  Ali Baharev, All rights reserved.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


/*
=========================================================================
This program is documented in the paper of A. Baharev, S. Kemeny,
On the computation of the noncentral F and noncentral beta distribution;
Statistics and Computing, 2008, 18 (3), 333-340.
http://dx.doi.org/10.1007/s11222-008-9061-3

Preprint of this paper is available at http://reliablecomputing.eu

My e-mail address is: (my first name) dot (my last name) at gmail dot com

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
#include "nmath.h"
#include "dpq.h" 
*/

/* This program depends on the following functions of R */
double pf(double x, double n1, double n2, int lower_tail, int log_p);
double qbeta(double alpha, double p, double q, int lower_tail, int log_p);
double dbeta(double x, double a, double b, int give_log);
double pbeta(double x, double pin, double qin, int lower_tail, int log_p);
double dpois(double x, double lambda, int give_log);
double qpois(double p, double lambda, int lower_tail, int log_p);

/*
2-moment central F approximation; Patnaik P. B. 1949.
The non-central chi-square and F distribution and their applications;
Biometrika; 36: 202-232.
*/
static double patnaik2(double x, double nu1, double nu2, double lambda) {

  return pf( x/(1+lambda/nu1), (nu1+lambda)*((nu1+lambda)/(nu1+2*lambda)), nu2, 1, 0);
}

/*
This function gives an initial value of lambda for the Newton iteration.
First, the lambda value is bracketed, then bisection is used to find a 
better approximation. This function uses the 2-moment central F 
approximation of Patnaik.
*/
static double guess(double prob, double y, double nu1, double nu2) {

  double x;
  double lambdal;
  double lambdam;
  double lambdau;
  double fl;
  double fm;
  double fu;
  int itr_cnt;

  /* FIXME: cancellation ? */
  x = nu2*y/(nu1*(1.0-y));
  lambdal = 0.0;
  lambdau = 1.0;

  fl = pf(x, nu1, nu2, 1, 0);

  /* In this case there is no solution
  FIXME: how this error is handled properly in R ? */
  if (fl < prob) {
    fprintf( stderr, "No solution (probably incorrect input)!\n");
    exit(127);
  }

  fu = patnaik2(x, nu1, nu2, lambdau);

  /* Bracketing lambda: lambdal <= lambda <= lambdau */
  for (itr_cnt=1; ((fl-prob)*(fu-prob)>0.0)&&itr_cnt<=17; ++itr_cnt)
  {
    fl = fu;
    lambdal = lambdau;
    lambdau = 2.0*lambdau;
    fu = patnaik2(x,nu1,nu2,lambdau);
  }

  /* FIXME: how this error is handled properly in R ? */
  if (itr_cnt == 18) {
    printf("\nlambdau = %f\n", lambdau);
    printf("\nfl = %f\n", fl);
    printf("\nfu = %f\n", fu);
    fprintf( stderr, "Failed to bracket lambda!\n");
    exit(127);
  }

  /* find a better approximation of lambda by bisection */

  lambdam = (lambdal + lambdau)/2.0;

  for (itr_cnt=1; ((lambdau-lambdal)>1.0e-4*lambdau)&&itr_cnt<=29; ++itr_cnt)
  {

    fm = patnaik2(x, nu1, nu2, lambdam);

    if ((fm-prob)*(fu-prob) < 0.0) {
      fl = fm;
      lambdal = lambdam;
    }
    else {
      fu = fm;
      lambdau = lambdam;
    }

    lambdam = (lambdal + lambdau)/2.0;

  }

  /* FIXME: how this error is handled properly in R ? */
  if (itr_cnt == 30) {
    fprintf( stderr, "Failed to find initial guess!\n");
    exit(127);
  }

  return lambdam;
}

/*
Given prob, x, a and b, this function returns the corresponding 
noncentrality parameter of the noncentral beta distribution.

I.e. the following equation

I_x(a, b, lambda) = prob

is solved for lambda with Newton iteration.

This function works just fine when supplied with meaningful input
data (and from practically meaningful range) but may easily crash
if not. Please be nice.
*/
double ncbeta(double prob, double x, double a, double b) {

  double ql;
  double qu;
  double c;
  double d;
  double p;
  double lambda;
  double lambda_new;
  double k;
  double f;
  double g;
  double mu;
  double eps;
  double eps2;
  int itr_cnt;

  lambda_new = guess(prob, x, 2.0*a, 2.0*b);

  /* FIXME: are these tolerances OK ?  */
  eps  = 1.0e-7;
  eps2 = 1.0e-6;

  itr_cnt = 0;

  do {

    lambda = lambda_new;

    mu = lambda/2.0;

    ql = qpois(eps, mu, 1, 0);

    qu = qpois(eps, mu, 0, 0);

    k = qu;

    c = pbeta(x, a+k, b, 1, 0);

    d = x*(1.0-x)/(a+k-1.0)*dbeta(x, a+k-1, b, 0);

    p = dpois(k, mu, 0);

    f=p*c;

    p = k/mu*p;

    g = p*d;

    for (k = qu-1; k >= ql; --k) {

      c=c+d;

      d=(a+k)/(x*(a+k+b-1))*d;

      f=f+p*c;

      p=k/mu*p;

      g=g+p*d;

    }

    /* Newton step */
    lambda_new = lambda+2.0*(f-prob)/g;

    ++itr_cnt;
  }
  while ((fabs(lambda_new-lambda) > eps2*lambda_new)&&(itr_cnt<=10));

  /* FIXME: how this error is handled properly in R ? */
  if (itr_cnt == 11) {
    fprintf( stderr, "Newton iteration failed in ncbeta()!\n");
    exit(127);
  }

  return lambda_new;

}

int main() {


  double x, a, b, param;
  double nu1[1+ 9];
  double nu2[1+26];
  int i, j;

  printf("\nProgram for calculating minimal detectable differences ");
  printf("for general ANOVA models\n");
  printf("Copyright (C) 2008, Ali Baharev, All rights reserved.\n");
  printf("This program comes with ABSOLUTELY NO WARRANTY.\n");
  printf("This is free software, and you are welcome to redistribute it\n");
  printf("under certain conditions (GNU GPL).\n\n");
  printf("The output of this program is the corrected table of \n");
  printf("Lorenzen and Anderson (1993) Appendix 12, p. 374\n\n");

  for (i=1; i<=6; ++i)
    nu1[i] = i;

  nu1[7] = 10;
  nu1[8] = 20;
  nu1[9] = 50;

  for (i=1; i<=8; ++i)
    nu2[i] = i;

  for (i=9; i<=19; ++i)
    nu2[i] = 2*i-8;

  for (i=20; i<=23; ++i)
    nu2[i] = 20*(i-18);

  nu2[24] = 200;
  nu2[25] = 500;
  nu2[26] = 1000;

  for (j=1; j<=26; ++j) {
    for (i=1; i<=9; ++i) {

      a = nu1[i]/2.0;
      b = nu2[j]/2.0;
      /* Type  I error probability 0.05
         Type II error probability 0.10 */

      x = qbeta(0.95, a, b, 1, 0);
      param = sqrt( ncbeta( 0.10, x, a, b)/nu1[i]);

      printf("%f\t", param);
    }
    printf("\n");
  }


  return 0;

}

