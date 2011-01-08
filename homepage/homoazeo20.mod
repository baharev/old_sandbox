#################### DESCRIPTION ###############################################
#
# Last updated: 21 Aug 2009
#
# The homoazeo40.mod has three solutions corresponding to the lower 
# concentration steady state (LSS), upper stable steady state or higher steady 
# state (HSS), and unstable steady state (USS). The distillate molar flow rate, 
# denoted by D, can be chosen as bifurcation parameter.
# The homoazeo20.mod has a single solution, and also differs from 
# homoazeo40.mod in the number of theoretical stages and in the 
# feed stage location, besides the bifurcation parameter value.
#
# The unstable steady state is missed by state-of-the-art methods such as 
# the inside-out algorithm or the simultaneous correction procedure, even when 
# supplied with a good initial guess obtained with case study. Bifurcation 
# diagrams were computed with continuation methods earlier, for example:
#
# Arjun Vadapalli, J. D. Seader;
# A generalized framework for computing bifurcation diagrams using process 
# simulation programs;
# Computers and Chemical Engineering 25 (2001) 445–464
#
# A. Kannan, Manish R. Joshi, G. Ramalinga Reddy, and Denish M. Shah;
# Multiple-Steady-States Identification in Homogeneous Azeotropic Distillation 
# Using a Process Simulator;
# Ind. Eng. Chem. Res., 2005, 44 (12), 4386-4399; DOI: 10.1021/ie049443s
#
# Güttinger, Dorn and Morari (see below) used the following software to compute
# the unstable steady state:
#
# Doedel, E. J.; Wang, X. AUTO94: Software for Continuation and 
# Bifurcation Problems in Ordinary Differential Equations; Computer
# Science Department of Concordia University: Montreal, Canada, 1994
#
# Except the number of stages and the feed stage location, the model and its 
# parameters correspond to the "Auto model" discussed in:
#
# Thomas E. Güttinger, Cornelius Dorn, Manfred Morari;
# Experimental Study of Multiple Steady States in Homogeneous Azeotropic
# Distillation;
# Ind. Eng. Chem. Res. 1997, 36, 794-802
#
# Solving this model with interval methods seems to be first discussed in:
#
# A. Baharev, E. Rév;
# A complete nonlinear system solver using affine arithmetic;
# Interval Analysis and Constraint Propagation for Applications (IntCP 2009);
# Workshop held in conjunction with the 15th International Conference on 
# Principles and Practice of Constraint Programming (CP 2009);
# Lisbon, Portugal, September 20th, 2009
#
# Manuscript of this talk is available from:
#
# http://reliablecomputing.eu
#
# The above interval method finds the three solutions in 9 s without any user 
# provided initial guesses.
#
################### PROBLEM STATISTICS #########################################
#
# Presolve eliminates 60 constraints and 60 variables.
# Substitution eliminates 380 variables.
# Adjusted problem:
# 140 variables, all nonlinear
# 140 constraints; 594 nonzeros
#	  120 nonlinear constraints
#	  20 linear constraints
# 0 objectives.
#
################# PARAMETERS FROM SPECIFICATION ################################

# NUMBER OF STAGES
param N := 20;

# FEED STAGE LOCATION
param N_F := 10;

# NUMBER OF COMPONENTS
param C := 3;

# DISTILLATE MOLAR FLOW RATE
param D;

let D := 0.42;

# VAPOR FLOW RATE, FROM SPECIFICATION
param V := 1.38;

# FEED LOCATION AND VALUES ARE GIVEN IN THE DATA SECTION

param F{j in 0..N};
param f{i in 1..C, j in 0..N};

# AUXILIARY PARAMETERS

param L{j in 0..N} = V - D + sum{k in 0..j} F[k];
param B := F[N_F] - D;

################## PARAMETERS ##################################################

param a{1..C};
param b{1..C};
param c{1..C};

param r{1..C, 1..C};
param s{1..C, 1..C};

param P := 100000.0;

# LOWER/UPPER BOUNDS ON THE VARIABLES, INITIAL ESTIMATES (IF NEEDED)
# VALUES ARE GIVEN IN THE DATA SECTION

param x_L{1..C, 1..N}; param x_U{1..C, 1..N}; param x_0{1..C, 1..N};
param K_L{1..C, 1..N}; param K_U{1..C, 1..N}; param K_0{1..C, 1..N};

param T_0{1..N};

############### VARIABLES ######################################################

var x{i in 1..C, j in 1..N} >= x_L[i,j], <= x_U[i,j], := x_0[i,j];

var K{i in 1..C, j in 1..N} >= K_L[i,j], <= K_U[i,j], := K_0[i,j];

var T{j in 1..N} >= 336.3, <= 383.4, := T_0[j];

####### DEFINED VARIABLES (THEY ARE ELIMINATED BY PRESOLVE / SUBTITUTION) ######

var p{i in 1..C, j in 1..N} = exp(a[i]+b[i]/(T[j]+c[i]));

var rcp_T{j in 1..N} = 1.0/T[j];

var Lambda{i1 in 1..C, i2 in 1..C, j in 1..N} = exp(r[i1,i2]+s[i1,i2]*rcp_T[j]);

var sum_xLambda{i in 1..C, j in 1..N} = sum{i1 in 1..C} (x[i1,j]*Lambda[i,i1,j]);

var rcp_sum_xLambda{i in 1..C, j in 1..N} = 1.0/sum_xLambda[i,j];

var gamma{i in 1..C, j in 1..N} =
  exp( -log(sum_xLambda[i,j]) + 1.0 - (sum{i2 in 1..C} (x[i2,j]*Lambda[i2,i,j]*rcp_sum_xLambda[i2,j])) );
  
############## EQUATIONS #######################################################

# AUXILIARY EQUATIONS

E_aux_K{j in 1..N, i in 1..C}:
	K[i,j] - gamma[i,j]*(p[i,j]/P) = 0.0;

# MATERIAL BALANCES

M_tot{i in 1..C}:
	D*(K[i,1]*x[i,1]) + B*x[i,N] - f[i,N_F] = 0.0;

# NOTE THE UNUSUAL FORMULATION
M_eq{j in 1..N-1, i in 1..C}:
	L[j]*x[i,j] + sum{i1 in j+1..N} f[i,i1] - B*x[i,N] - V*(K[i,j+1]*x[i,j+1]) = 0.0;

# SUMMATION EQUATIONS

S_x_eq{j in 1..N}:
	sum{i in 1..C} x[i,j] - 1.0 = 0.0;

################### DATA SECTION ###############################################

data;

let a[1] := 23.4832;
let a[2] := 20.5110;
let a[3] := 20.9064;

let b[1] := -3634.01;
let b[2] := -2664.30;
let b[3] := -3096.52;

let c[1] := -33.768;
let c[2] := -79.483;
let c[3] := -53.668;

let r[1,2] :=  0.7411;
let r[1,3] :=  0.9645;
let r[2,3] := -1.4350;

let r[2,1] := -1.0250;
let r[3,1] := -0.9645;
let r[3,2] :=  2.7470;

let r[1,1] := 0.0;
let r[2,2] := 0.0;
let r[3,3] := 0.0;

let s[1,2] := -477.00;
let s[1,3] := -903.1024;
let s[2,3] :=  768.20;

let s[2,1] :=  72.78;
let s[3,1] := -140.9995;
let s[3,2] := -1419.0;

let s[1,1] := 0.0;
let s[2,2] := 0.0;
let s[3,3] := 0.0;

# LOWER AND UPPER BOUNDS ON THE VARIABLES

for {j in 1..N} {
	let x_L[1,j] := 0.0;
	let x_U[1,j] := 0.9998;
	let x_L[2,j] := 1.0e-4;
	let x_U[2,j] := 0.9999;
	let x_L[3,j] := 1.0e-4;
	let x_U[3,j] := 0.9999;
}

# THIS BOUND SEEMS TO BE REASONABLE FROM ENGINEENERING POINT OF VIEW

let x_L[1,1] := 0.83;

# THESE BOUNDS WERE CALCULATED IN ADVANCE

for {j in 1..N} {
	let K_L[1,j] := 0.9784;
	let K_L[2,j] := 0.2445;
	let K_L[3,j] := 0.2745;
	
	let K_U[1,j] := 40.52;
	let K_U[2,j] := 1.317;
	let K_U[3,j] := 1.975;
}

# COMES FROM x_L[1,1] (K=y/x, y<1) 

let K_U[1,1] := 1.21;

# FEED VALUES, LOCATION

for {i in 1..C, j in 0..N}
	let f[i,j] := 0.0;

let f[1, N_F] := 0.4098370;
let f[2, N_F] := 0.01229769;
let f[3, N_F] := 0.06090665;

for {j in 0..N}
	let F[j] := sum{i in 1..C} f[i,j];

################################################################################

# DUMB INITIAL ESTIMATES (IF NEEDED)

for {i in 1..C, j in 1..N}
	let x_0[i,j] := 0.33;

for {j in 1..N}
	let T_0[j] := 337.0;

for {i in 1..C, j in 1..N}
	let K_0[i,j] := 1.0;

# for {j in N_F+1..N} {
	# let K_0[1,j] := 40.52;
# }

################################################################################

option show_stats 1;
option presolve 10;
option substout 1;
option var_bounds 2;
option nl_comments 0;
#option nl_permute 0;
#option auxfiles cfru;
#write ghomoazeo20;

# option solver ipopt;
solve;

# display x, K, T;
#
# :          x           K        :=
# 1 1    0.885425     0.998633
# 1 2    0.886609     0.998249
# 1 3    0.887662     0.997992
# 1 4    0.888442     0.997942
# 1 5    0.888758     0.998196
# 1 6    0.88836      0.998891
# 1 7    0.886901     1.00022
# 1 8    0.883855     1.00252
# 1 9    0.878327     1.00642
# 1 10   0.868552     1.01332
# 1 11   0.868906     1.01317
# 1 12   0.869448     1.01297
# 1 13   0.87023      1.01271
# 1 14   0.871248     1.01247
# 1 15   0.872282     1.01249
# 1 16   0.872525     1.01344
# 1 17   0.869717     1.01701
# 1 18   0.85776      1.02776
# 1 19   0.816329     1.06461
# 1 20   0.610188     1.35326
# 2 1    0.00184681   0.587416
# 2 2    0.00273454   0.590562
# 2 3    0.00376075   0.593621
# 2 4    0.00494088   0.59632
# 2 5    0.00629679   0.598291
# 2 6    0.0078634    0.599048
# 2 7    0.00970051   0.597945
# 2 8    0.0119149    0.594077
# 2 9    0.0147065    0.586051
# 2 10   0.0184815    0.571425
# 2 11   0.0187737    0.572323
# 2 12   0.01926      0.573736
# 2 13   0.020071     0.575887
# 2 14   0.0214283    0.578985
# 2 15   0.0237153    0.583
# 2 16   0.0276246    0.587066
# 2 17   0.0345271    0.5881
# 2 18   0.0476505    0.577605
# 2 19   0.0774215    0.532747
# 2 20   0.187846     0.3853
# 3 1    0.112729     1.0175
# 3 2    0.110656     1.02415
# 3 3    0.108577     1.03049
# 3 4    0.106617     1.03586
# 3 5    0.104945     1.03938
# 3 6    0.103776     1.03988
# 3 7    0.103399     1.03581
# 3 8    0.104231     1.02502
# 3 9    0.106967     1.00421
# 3 10   0.112967     0.967722
# 3 11   0.11232      0.969563
# 3 12   0.111292     0.972438
# 3 13   0.109699     0.976766
# 3 14   0.107323     0.982859
# 3 15   0.104002     0.990359
# 3 16   0.0998501    0.996767
# 3 17   0.0957556    0.994043
# 3 18   0.0945892    0.961037
# 3 19   0.10625      0.844087
# 3 20   0.201966     0.504426
# ;
#
# T [*] :=
 # 1  336.413
 # 2  336.424
 # 3  336.435
 # 4  336.449
 # 5  336.464
 # 6  336.482
 # 7  336.501
 # 8  336.523
 # 9  336.55
# 10  336.588
# 11  336.591
# 12  336.596
# 13  336.605
# 14  336.62
# 15  336.644
# 16  336.687
# 17  336.759
# 18  336.89
# 19  337.181
# 20  338.479
# ;
#
# END OF FILE
