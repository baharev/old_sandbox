/*********************************************************************
*  Rigorous lower and upper bounds in linear programming
*
*  Copyright (C) 2008 Ali Baharev, All rights reserved. 
*  http://reliablecomputing.eu
*
*  This program is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

/* A very rough implementation of the algorithm of Neumaier and Shcherbina, 
*  Math. Programming A 99 (2004), 283-296. 
*
*  WARNING!
*
*  Only tested with Microsoft Visual Studio 2005 on Win32 XP
*  and gcc 4.1 on Ubuntu 7.10.
*  If you have problems with compilation check functions
*  rounddown(), roundup(), roundnear() and basicHexImage().
*  Other functions should work properly.
*
*  In the previous version (last updated on 16 March 2008), there was a serious
*  bug in function validate_obj, the rounding mode round_near() was set BEFORE 
*  calculating the returned value.
*
*  Last updated: 07 Aug 2008
*  
*/

#ifdef _MSC_VER
#define _CRTDBG_MAP_ALLOC
#endif

#include <stdlib.h>
#include <float.h>

#ifdef _MSC_VER
#pragma fenv_access (on)
#endif

#include <stdio.h>
#include <glpk.h>

#ifdef _MSC_VER
#include <crtdbg.h>
#endif

#ifdef __GNUC__
#include <fenv.h>
#pragma STDC FENV_ACCESS ON
#endif

void rounddown() {
	int err;
#ifdef _MSC_VER
	unsigned int control_word;
	err = _controlfp_s(&control_word, _RC_DOWN, _MCW_RC);
#endif
#ifdef __GNUC__
	err = fesetround(FE_DOWNWARD);
#endif
	if ( err ) {
		fprintf(stderr, "Failed to change rounding mode!\n");
		exit(EXIT_FAILURE);
	}
}

void roundup() {    
	int err;
#ifdef _MSC_VER
	unsigned int control_word;
	err = _controlfp_s(&control_word, _RC_UP, _MCW_RC);
#endif
#ifdef __GNUC__
	err = fesetround(FE_UPWARD);
#endif
	if ( err ) {
		fprintf(stderr, "Failed to change rounding mode!\n");
		exit(EXIT_FAILURE);
	}
}

void roundnear() {
    int err;
#ifdef _MSC_VER
	unsigned int control_word;
	err = _controlfp_s(&control_word, _RC_NEAR, _MCW_RC);
#endif
#ifdef __GNUC__
	err = fesetround(FE_TONEAREST);
#endif
	if ( err ) {
		fprintf(stderr, "Failed to change rounding mode!\n");
		exit(EXIT_FAILURE);
	}
}
/* This code building the hexadecimal image of a double data 
*  is copied from the source of FILIB++ and slightly changed,
*  see http://xsc.de    */
typedef union 
{
	double f;

	struct 
	{
		#if defined(__sparc)
			#if defined(BYTE_ORDER)
			#undef BYTE_ORDER
			#endif
			#if defined(LITTLE_ENDIAN)
			#undef LITTLE_ENDIAN
			#endif
			#if defined(BIG_ENDIAN)
			#undef BIG_ENDIAN
			#endif
			#define LITTLE_ENDIAN 1234
			#define BIG_ENDIAN 4321
			#define BYTE_ORDER BIG_ENDIAN
		#endif

		#if !defined(BYTE_ORDER) && defined(__i386__)
		#define BYTE_ORDER LITTLE_ENDIAN
		#endif

		#if !defined(BYTE_ORDER) && defined(WIN32)
		#define BYTE_ORDER LITTLE_ENDIAN
		#endif

		#if !defined(BYTE_ORDER)
		#error "Undefined byte order."
		#endif

		#if (BYTE_ORDER == LITTLE_ENDIAN)
			unsigned int mant1 :32;
			unsigned int mant0 :20;
			unsigned int expo  :11;
			unsigned int sign  : 1;
		#elif (BYTE_ORDER == BIG_ENDIAN)
			unsigned int sign  : 1;
			unsigned int expo  :11;
			unsigned int mant0 :20;
			unsigned int mant1 :32;
		#elif (BYTE_ORDER == PDP_ENDIAN)
			#error "Sorry, byte order for PDP is unimplemented."
		#else
			#error "Sorry, your machines byte order is unknown."
		#endif
	} ieee;
} a_diee;

void basicHexImage(double d) 
{
	a_diee f;

	char sigBuf[2];
	char expBuf[4];
	char mant0Buf[6];
	char mant1Buf[9];
	char delim[2];

	delim[0] = ':';
	delim[1] = '\0';

	f.f = d;
	sigBuf[0] = (f.ieee.sign == 1 ? '1' : '0');
	sigBuf[1] = '\0';

	printf("%s", sigBuf);

	printf("%s", delim);

	sprintf(expBuf,"%03x",f.ieee.expo);

	printf("%s", expBuf);

	printf("%s", delim);

	sprintf(mant0Buf,"%05x",f.ieee.mant0);

	printf("%s", mant0Buf);

	sprintf(mant1Buf,"%08x",f.ieee.mant1);

	printf("%s", mant1Buf);
}

/* A messy method to get the transpose of A
*  in rowwise representation                  */
void build_transpose(glp_prob* lp, int* jt, int* nz_jt, double* at) {

	int n;
	int j;

	n = glp_get_num_cols(lp);

	nz_jt[0] = 0;
	jt[0] = -1;
	at[0] =  0;

	for (j=1; j<=n; ++j)
		nz_jt[j] = nz_jt[j-1] + 
		glp_get_mat_col(lp, j, jt+nz_jt[j-1], at+nz_jt[j-1]);
}

void compute_r(const int n, 
				double* r,
				int* jt, 
				int* nz_jt, 
				double* at, 
				double* ds, 
				double* dr, 
				double* cs)
{

	int i, j, start, end;

	/* r[j] = A_T*d_R+d_S-c_S */
	for (j=1; j<=n; ++j) {
		/* In the preprint, and probably in the final paper, there is a typo,
		   c_S has negative sign, but in the preprint it has positive sign  */ 
		r[j] = ds[j]-cs[j];

		start = nz_jt[j-1]+1;
		end   = nz_jt[j  ];

		for (i=start; i<=end; ++i)
			r[j] += at[i]*dr[jt[i]];
	}
}

double validate_obj(glp_prob* lp) 
{

	/*
	*	n: number of columns
	*	m: number of rows
	*	n_nzeros: number of non-zero elements in the constraint matrix A
	*	jt: indices of non-zero elements in the j'th column
	*	nz_jt: number of non-zero elements in the j'th column
	*	at: values of the non-zero elements in the transpose of A
	*   For further notations, please consult the GLPK reference manual at
	*   checking Karush-Kuhn-Tucker conditions, and the paper of Neumaier and 
	*   Shcherbina.
	*/

	int i, j, m, n, n_nzeros;
	double mu, d_b, r_bs, bs_l, bs_u, z;
	int* jt;
	int* nz_jt;
	double* at;
	double* dr;
	double* ds;
	double* cs;
	double* r_lb;
	double* r_ub;

	n = glp_get_num_cols(lp);
	m = glp_get_num_rows(lp);
	n_nzeros = glp_get_num_nz(lp);

	jt    = (int*) malloc((1+n_nzeros)*sizeof(int));
	nz_jt = (int*) malloc((1+n       )*sizeof(int));
	at = (double*) malloc((1+n_nzeros)*sizeof(double));
	dr = (double*) malloc((1+n_nzeros)*sizeof(double));
	ds = (double*) malloc((1+n_nzeros)*sizeof(double));
	cs = (double*) malloc((1+n_nzeros)*sizeof(double));
	r_lb = (double*) malloc((1+n)*sizeof(double));
	r_ub = (double*) malloc((1+n)*sizeof(double));

	/* Get d_R, d_S, c_S */

	for (i=1; i<=m; ++i)
		dr[i] = glp_get_row_dual(lp, i);

	for (j=1; j<=n; ++j) {
		ds[j] = glp_get_col_dual(lp, j);
		cs[j] = glp_get_obj_coef(lp, j);
	}

	/* A messy method to get the transpose of A */
	build_transpose(lp, jt, nz_jt, at);

	/* From here, the steps at (11) of Neumaier and Shcherbina are performed */
	
	rounddown();

	mu = 0.0;

	for (i=1; i<=m; ++i) {
		if (dr[i] > 0.0)
			mu += dr[i]*glp_get_row_lb(lp, i);
	}

	for (j=1; j<=n; ++j) {
		if (ds[j] > 0.0)
			mu += ds[j]*glp_get_col_lb(lp, j);
	}

	compute_r(n, r_lb, jt, nz_jt, at, ds, dr, cs);

	roundup();
	
	compute_r(n, r_ub, jt, nz_jt, at, ds, dr, cs);

	for (j=1; j<=n; ++j) {
		r_ub[j] = ( (-r_lb[j]) > r_ub[j] ) ? (-r_lb[j]) : r_ub[j] ;
	}

	d_b = 0.0;

	for (i=1; i<=m; ++i) {
		if (dr[i] < 0.0)
			d_b += (-dr[i])*glp_get_row_ub(lp, i);
	}

	for (j=1; j<=n; ++j) {
		if (ds[j] < 0.0)
			d_b += (-ds[j])*glp_get_col_ub(lp, j);
	}

	r_bs = 0.0;

	for (j=1; j<=n; ++j) {
		bs_l = glp_get_col_lb(lp, j);
		bs_u = glp_get_col_ub(lp, j);
		r_bs += r_ub[j]*(((-bs_l)>bs_u)?(-bs_l):bs_u);

	}
	
	/* This is the rigorous lower bound */
	z = -(d_b-mu+r_bs);

	free(jt);
	free(nz_jt);
	free(at);
	free(dr);
	free(ds);
	free(cs);
	free(r_lb);
	free(r_ub);

	roundnear();

	return z;
}

/* Run-time check if rounding works or not */
void test_rounding() {

	double lb, ub;

	printf("\nChecking rounding mode...\n");

	/* The textbook example of 0.1 */
	rounddown();
	lb = ((double) 1)/((double) 10);
	roundup();
	ub = ((double) 1)/((double) 10);
	roundnear();
	if (!(lb<ub)) {
		fprintf(stderr, "Error: rounding mode cannot be set!\n");
		exit(EXIT_FAILURE);
	}

	printf("Textbook example: 0.1 is not exactly representable\n");
	printf("Lower bound on 0.1 in hexadecimal format:\n");
	basicHexImage(lb);
	printf("\nUpper bound on 0.1 in hexadecimal format:\n");
	basicHexImage(ub);
	printf("\nThe bounds differ, rounding seems to be OK!\n\n");

}

void print_objective(glp_prob* lp) {

	double z;

	z = glp_get_obj_val(lp);
	printf("\nThe objective value calculated by glp_simplex in hexadecimal ");
	printf("format is: \n");
	basicHexImage(z);
	printf(" (decimal: %.15f)", z);

	/* The algorithm of Neumaier and Shcherbina is invoked here */
	z = validate_obj(lp);
	
	printf("\nRigorous lower bound on the objective value ");
	printf("in hexadecimal format is: \n");
	basicHexImage(z);
	printf(" (decimal: %.15f)\n\n", z);
}

int main() {

	glp_prob* lp;
	
	/*	We build the LP object through the C API of GLPK 
	*	All rows and columns must be double bounded or fixed!
	*	
	*	n: number of columns
	*	m: number of rows
	*	n_nzeros: number of non-zero elements in the constraint matrix A
	*	ia, ja: indices of non-zero elements (row, column) of A
	*	ar: values of the non-zero elements in A
	*/

	int n, m, n_nzeros;

	int* ia;
	int* ja;
	double* ar;

	printf("\nThis program calculates rigorous lower bound on the ");
	printf("objective value of an\n");
	printf("LP problem if all rows and columns are double-bounded or fixed.\n");
	printf("For details, see\n");
	printf("Neumaier and Shcherbina, ");
	printf("Math. Programming A 99 (2004), 283-296.\n");
	printf("Preprint of this paper is available at\n");
	printf("http://www.mat.univie.ac.at/~neum/papers.html\n\n");
	printf("--------------------------------------------------------------\n");

	test_rounding();

	printf("--------------------------------------------------------------\n");

	n        = 4;
	m        = 3;
	n_nzeros = 8;

	ia    = (int*) malloc((1+n_nzeros)*sizeof(int));
	ja    = (int*) malloc((1+n_nzeros)*sizeof(int));

	ar = (double*) malloc((1+n_nzeros)*sizeof(double));

	lp = glp_create_prob();

	glp_set_obj_dir(lp, GLP_MIN);
	glp_add_cols(lp, n);
	glp_add_rows(lp, m);

	ia[1] = 1;
	ja[1] = 1;
	ar[1] =-2;

	ia[2] = 1;
	ja[2] = 2;
	ar[2] = 2;

	ia[3] = 1;
	ja[3] = 3;
	ar[3] =-2;

	ia[4] = 1;
	ja[4] = 4;
	ar[4] = 2;

	ia[5] = 2;
	ja[5] = 1;
	ar[5] =-3;

	ia[6] = 2;
	ja[6] = 2;
	ar[6] = 4;

	ia[7] = 3;
	ja[7] = 3;
	ar[7] = 2;

	ia[8] = 3;
	ja[8] = 4;
	ar[8] = 4;

	glp_set_col_bnds(lp, 1, GLP_DB, -1, 1);
	glp_set_col_bnds(lp, 2, GLP_DB, -2, 2);
	glp_set_col_bnds(lp, 3, GLP_DB, -1, 3);
	glp_set_col_bnds(lp, 4, GLP_DB, -3, 1);

	glp_set_row_bnds(lp, 1, GLP_DB, -2,  1);
	glp_set_row_bnds(lp, 2, GLP_DB, -3,  2);
	glp_set_row_bnds(lp, 3, GLP_FX, -1, -1);

	glp_load_matrix(lp, n_nzeros, ia, ja, ar);

	free(ia);
	free(ja);
	free(ar);

	glp_set_obj_coef(lp, 1,-1);
	glp_set_obj_coef(lp, 2, 1);
	glp_set_obj_coef(lp, 3, 2);
	glp_set_obj_coef(lp, 4,-1);

	/* lpx_write_cpxlp(lp, "LP_problem.lp"); */

	printf("\nThe following LP problem has been built:\n");
	printf("Minimize:\n");
	printf("  z =   -   x_1 +   x_2 + 2 x_3 -   x_4\n");
	printf("Subject To\n");
	printf("  -2 <= - 2 x_1 + 2 x_2 - 2 x_3 + 2 x_4 <=  1\n");
	printf("  -3 <= - 3 x_1 + 4 x_2                 <=  2\n");
	printf("                          2 x_3 + 4 x_4  = -1\n");
	printf("Bounds\n");
	printf(" -1 <= x_1 <= 1\n");
	printf(" -2 <= x_2 <= 2\n");
	printf(" -1 <= x_3 <= 3\n");
	printf(" -3 <= x_4 <= 1\n\n");

	printf("--------------------------------------------------------------\n");
	printf("\nSolving the LP problem:\n\n");

	glp_simplex(lp, NULL);

	print_objective(lp);

	glp_set_obj_coef(lp, 1, 0.0);

	printf("--------------------------------------------------------------\n");
	printf("\nColumn x_1 is removed from the objective, ");
	printf("re-solving the LP problem:\n\n");

	glp_simplex(lp, NULL);

	print_objective(lp);

	glp_delete_prob(lp);

#ifdef _MSC_VER
	_CrtDumpMemoryLeaks();
#endif

	printf("\nSuccessfully finished!\n");

	return 0;

}
