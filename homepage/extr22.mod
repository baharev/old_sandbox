#################### DESCRIPTION ###############################################
#
# Last updated: 22 Aug 2009
#
# These models have a single solution corresponding to the steady state of an 
# extractive distillation column with 22 and 30 theoretical stages, 
# respectively. The model (MESH equations) and the notion are discusses in 
# greater detail in:
#
# A. Baharev, T. Achterberg, E. Rév;
# Computation of an extractive distillation column with affine arithmetic;
# AIChE Journal, 2009, 55 (7), 1695-1704.
# (preprint available from http://reliablecomputing.eu)
#
# The above cited paper seems to be the first paper on computing distillation 
# columns with interval methods. More recent results are given in:
#
# A. Baharev, E. Rév;
# A complete nonlinear system solver using affine arithmetic;
# Interval Analysis and Constraint Propagation for Applications (IntCP 2009);
# Workshop held in conjunction with the 15th International Conference on 
# Principles and Practice of Constraint Programming (CP 2009);
# Lisbon, Portugal, September 20th, 2009.
# (manuscript available from http://reliablecomputing.eu)
#
#
############## PROBLEM STATISTICS ##############################################
#
# Presolve eliminates 68 constraints and 68 variables.
# Substitution eliminates 463 variables.
# Adjusted problem:
# 387 variables:
#	  242 nonlinear variables
#	  145 linear variables
# 387 constraints; 1326 nonzeros
#	  286 nonlinear constraints
#	  101 linear constraints
# 0 objectives.
#
#
################## PARAMETERS FROM SPECIFICATION ###############################

# NUMBER OF STAGES
param N := 22;

# FEED STAGE LOCATIONS
param N_F1 := 9;
param N_F2 := 16;

# NUMBER OF COMPONENTS
param C := 3;

# REFLUX RATIO
param R := 5.0;

# DISTILLATE MOLAR FLOW RATE
param D := 0.73;

# PURITY RESTRICTION ON THE DISTILLATE COMPOSITION
param x_acetone_min := 0.92;

# LOWER BOUND ON ALL MOL FRACTIONS IN THE ENTIRE COLUMN
param xmin := 0.001;

################## PARAMETERS ##################################################

# FEED LOCATIONS AND VALUES ARE GIVEN IN THE DATA SECTION

param F{j in 0..N};
param f{i in 1..C, j in 0..N};

# FURTHER MODEL PARAMETERS

param lambda{1..C};
param a{1..C};
param b{1..C};
param c{1..C};
param Vm{1..C};
param k{1..C, 1..C};

param RG := 1.98721;
param P := 760.0;

# LOWER/UPPER BOUNDS ON THE VARIABLES, INITIAL ESTIMATES (IF NEEDED)
# VALUES ARE GIVEN IN THE DATA SECTION

param x_L{1..C, 0..N};   param x_U{1..C, 0..N};   param x_0{1..C, 0..N};
param y_L{1..C, 1..N+1}; param y_U{1..C, 1..N+1}; param y_0{1..C, 1..N+1};
param K_L{1..C, 1..N};   param K_U{1..C, 1..N};   param K_0{1..C, 1..N};

param V_L{1..N+1} default (1-0.32)*(R+1)*D;
param V_U{1..N+1} default (1+0.32)*(R+1)*D;
param V_0{1..N+1};

param l_L{1..C, 0..N}; param l_U{1..C, 0..N}; param l_0{1..C, 0..N};

param v_L{1..C, 1..N+1};
param v_U{1..C, j in 1..N+1} default V_U[j];
param v_0{1..C, 1..N+1};

param T_0{1..N};
param Q_L;	param Q_U;	param Q_0;

############### VARIABLES ######################################################

var x{i in 1..C, j in 0..N  } >= x_L[i,j], <= x_U[i,j], := x_0[i,j];

var y{i in 1..C, j in 1..N+1} >= y_L[i,j], <= y_U[i,j], := y_0[i,j];

var K{i in 1..C, j in 1..N  } >= K_L[i,j], <= K_U[i,j], := K_0[i,j];

var V{j in 1..N+1} >= V_L[j], <= V_U[j], := V_0[j];

var v{i in 1..C, j in 1..N+1} >= v_L[i,j], <= v_U[i,j], := v_0[i,j];

var l{i in 1..C, j in 0..N  } >= l_L[i,j], <= l_U[i,j], := l_0[i,j];

var Q >= Q_L, <= Q_U, := Q_0;

var T{j in 1..N} >= 327.0, <= 374.0, := T_0[j];

####### DEFINED VARIABLES (THEY ARE ELIMINATED BY PRESOLVE / SUBTITUTION) ######

var L{j in 0..N} = V[j+1] - D + sum{i1 in 0..j} F[i1];

var H{j in 1..N+1} = sum{i in 1..C} lambda[i]*y[i,j];

var rcp_T{j in 1..N} = 1.0/T[j];

var p{i in 1..C, j in 1..N} = exp(a[i]-b[i]/(T[j]+c[i]));

var Lambda{i1 in 1..C, i2 in 1..C, j in 1..N} =
  Vm[i2]/Vm[i1]*exp((-k[i1,i2]/RG)*rcp_T[j]);

var sum_xLambda{i in 1..C, j in 1..N} = sum{i1 in 1..C} (x[i1,j]*Lambda[i,i1,j]);

var rcp_sum_xLambda{i in 1..C, j in 1..N} = 1.0/sum_xLambda[i,j];

var gamma{i in 1..C, j in 1..N} =
  exp( -log(sum_xLambda[i,j]) + 1.0 - (sum{i2 in 1..C} (x[i2,j]*Lambda[i2,i,j]*rcp_sum_xLambda[i2,j])) );
  
############## EQUATIONS #######################################################

M_aux_y1_x0{i in 1..C}:
	x[i,0] - y[i,1] = 0.0;

M_aux_yNp1_xN{i in 1..C}:
	y[i,N+1] - x[i,N] = 0.0;
	
M_aux_V1:
	V[1] - (R+1.0)*D = 0.0;
	
M_aux_l{i in 1..C, j in 0..N}:
	l[i,j] - L[j]*x[i,j] = 0.0;

M_aux_v{i in 1..C, j in 1..N+1}:
	v[i,j] - V[j]*y[i,j] = 0.0;

M_eq{i in 1..C-1, j in 1..N}:
	l[i,j] + D*y[i,1] - v[i,j+1] - sum{i1 in 1..j} f[i,i1] = 0.0;

E_eq{i in 1..C, j in 1..N}:
	y[i,j] - K[i,j]*x[i,j] = 0.0;
	
S_x_eq{j in 1..N}:
	sum{i in 1..C} x[i,j] - 1.0 = 0.0;

S_y_eq{j in 1..N}:
	sum{i in 1..C} y[i,j] - 1.0 = 0.0;

H_eq{j in 1..N+1}:
	Q - V[j]*H[j] = 0.0;

E_aux_K{i in 1..C, j in 1..N}:
	K[i,j] - gamma[i,j]*(p[i,j]/P) = 0.0;

################### DATA SECTION ###############################################
	
data;

let lambda[1] := 6.960;
let lambda[2] := 8.426;
let lambda[3] := 9.717;

let a[1] := 16.732;
let a[2] := 18.510;
let a[3] := 18.304;

let b[1] := 2975.9;
let b[2] := 3593.4;
let b[3] := 3816.4;

let c[1] := -34.523;
let c[2] := -35.225;
let c[3] := -46.130;

let Vm[1] := 74.050;
let Vm[2] := 40.729;
let Vm[3] := 18.069;

for {i in 1..C, j in 0..N} {

	let x_L[i,j] := xmin;
	let x_U[i,j] := 1.0-(C-1)*xmin;
}

for {i in 1..C, j in 1..N+1} {
	
	let y_L[i,j] := xmin;
	let y_U[i,j] := 1.0-(C-1)*xmin;
}

let y_L[1,1] := x_acetone_min;
let y_U[2,1] := 1.00-y_L[1,1]-y_L[3,1];
let y_U[3,1] := 1.00-y_L[1,1]-y_L[2,1];

let x_L[1,0] := y_L[1,1];
let x_U[2,0] := y_U[2,1];
let x_U[3,0] := y_U[3,1];

for {j in 1..N} {
	let K_L[1,j] := 0.98;
	let K_U[1,j] := 38.97;
	let K_L[2,j] := 0.80;
	let K_U[2,j] := 7.53;
	let K_L[3,j] := 0.26;
	let K_U[3,j] := 1.01;
}

let Q_L := (R+1)*D*min(lambda[1], lambda[2], lambda[3]);
let Q_U := (R+1)*D*max(lambda[1], lambda[2], lambda[3]);

for {i in 1..C, j in 0..N}
	let f[i,j] := 0.0;

let f[3, N_F1] := 2.000;

let f[1, N_F2] := 0.783;
let f[2, N_F2] := 0.217;

for {j in 0..N}
	let F[j] := sum{i in 1..C} f[i,j];

for {i in 1..C, j in 1..N+1}
	let v_L[i,j] := 0.0;

for {i in 1..C, j in 0..N}
	let l_L[i,j] := 0.0;

for {i in 1..C, j in 0..N}
	let l_U[i,j] := sum{i1 in 0..j} f[i,i1] + v_U[i,j+1] - y_L[i,1]*D;

let k[1,2] := -157.981;
let k[1,3] :=  393.27;
let k[2,3] := -52.605;

let k[2,1] := 592.638;
let k[3,1] := 1430.0;
let k[3,2] := 620.63;

let k[1,1] := 0.0;
let k[2,2] := 0.0;
let k[3,3] := 0.0;

# DUMB INITIAL ESTIMATES (IF NEEDED)
	
for {i in 1..C, j in 0..N}
	let x_0[i,j] := (x_L[i,j]+x_U[i,j])/2.0;

for {i in 1..C, j in 1..N+1}
	let y_0[i,j] := (y_L[i,j]+y_U[i,j])/2.0;
	
for {i in 1..C, j in 1..N+1}
	let v_0[i,j] := (v_L[i,j]+v_U[i,j])/2.0;
	
for {i in 1..C, j in 0..N}
	let l_0[i,j] := (l_L[i,j]+l_U[i,j])/2.0;

for {i in 1..C, j in 1..N}
	let K_0[i,j] := (K_L[i,j]+K_U[i,j])/2.0;
	
for {j in 1..N+1}
	let V_0[j] := (V_L[j]+V_U[j])/2.0;
	
for {j in 1..N}
	let T_0[j] := 335.0;

let Q_0 := (Q_L + Q_U)/2.0;

################################################################################
	
option show_stats 1;
option presolve 10;
option substout 1;
option var_bounds 2;
option nl_comments 0;
option nl_permute 0;

solve;

############# SOLUTION #########################################################
#
# display x, y, K, V, v, l, T, Q;
#
# :          x           y           K        :=
# 1 0    0.960874        .           .
# 1 1    0.959312    0.960874     1.00163
# 1 2    0.955497    0.959573     1.00427
# 1 3    0.949355    0.956396     1.00742
# 1 4    0.940351    0.951285     1.01163
# 1 5    0.927161    0.943799     1.01795
# 1 6    0.906743    0.93285      1.02879
# 1 7    0.87122     0.915938     1.05133
# 1 8    0.793406    0.88662      1.11749
# 1 9    0.496192    0.822866     1.65836
# 1 10   0.490846    0.81739      1.66527
# 1 11   0.484328    0.810775     1.67402
# 1 12   0.4763      0.802706     1.68529
# 1 13   0.466426    0.792763     1.69966
# 1 14   0.455152    0.780512     1.71484
# 1 15   0.449474    0.766428     1.70517
# 1 16   0.498428    0.758651     1.52209
# 1 17   0.492984    0.753777     1.52901
# 1 18   0.483645    0.74569      1.54181
# 1 19   0.465843    0.731821     1.57096
# 1 20   0.421191    0.705449     1.67489
# 1 21   0.263107    0.639602     2.43095
# 1 22   0.0359304   0.402027    11.1891
# 1 23       .       0.0359304       .
# 2 0    0.021219        .           .
# 2 1    0.0172358   0.021219     1.2311
# 2 2    0.0145837   0.0179004    1.22743
# 2 3    0.0128556   0.0156926    1.22068
# 2 4    0.0117756   0.0142565    1.21069
# 2 5    0.011166    0.0133623    1.19669
# 2 6    0.0109282   0.0128625    1.177
# 2 7    0.0110392   0.0126762    1.14829
# 2 8    0.0115696   0.0127878    1.10529
# 2 9    0.0120494   0.0132671    1.10106
# 2 10   0.0177664   0.0194769    1.09627
# 2 11   0.024715    0.0269559    1.09067
# 2 12   0.033255    0.0360515    1.08409
# 2 13   0.0438889   0.0472386    1.07632
# 2 14   0.0573555   0.0611806    1.06669
# 2 15   0.0749975   0.0788449    1.0513
# 2 16   0.100742    0.101908     1.01157
# 2 17   0.106487    0.107351     1.00811
# 2 18   0.115962    0.116282     1.00276
# 2 19   0.131648    0.131034     0.99534
# 2 20   0.1568      0.155544     0.991995
# 2 21   0.178607    0.195366     1.09383
# 2 22   0.088771    0.233543     2.63085
# 2 23       .       0.088771        .
# 3 0    0.017907        .           .
# 3 1    0.0234521   0.017907     0.763556
# 3 2    0.0299191   0.0225269    0.752925
# 3 3    0.0377889   0.0279116    0.738618
# 3 4    0.0478734   0.0344586    0.719786
# 3 5    0.061673    0.0428384    0.694605
# 3 6    0.0823285   0.0542874    0.6594
# 3 7    0.117741    0.0713858    0.606295
# 3 8    0.195025    0.100593     0.515794
# 3 9    0.491759    0.163867     0.333226
# 3 10   0.491387    0.163133     0.331985
# 3 11   0.490957    0.162269     0.330515
# 3 12   0.490445    0.161243     0.328769
# 3 13   0.489686    0.159999     0.326738
# 3 14   0.487492    0.158308     0.324739
# 3 15   0.475529    0.154727     0.325379
# 3 16   0.400831    0.139441     0.347881
# 3 17   0.400529    0.138872     0.346722
# 3 18   0.400393    0.138028     0.34473
# 3 19   0.40251     0.137145     0.340726
# 3 20   0.422009    0.139007     0.329392
# 3 21   0.558285    0.165032     0.295605
# 3 22   0.875299    0.36443      0.416349
# 3 23       .       0.875299        .
# ;
#
# V [*] :=
 # 1 4.38       5 4.34469    9 4.14969   13 4.12795   17 4.1115    21 4.00412
 # 2 4.37511    6 4.3259    10 4.14574   14 4.11925   18 4.10561   22 3.71216
 # 3 4.36791    7 4.29764   11 4.14096   15 4.11045   19 4.09515   23 3.24489
 # 4 4.35807    8 4.24985   12 4.13513   16 4.11502   20 4.07294
# ;
#
# :          v           l        :=
# 1 0        .       3.50719
# 1 1    4.20863     3.4968
# 1 2    4.19823     3.47602
# 1 3    4.17745     3.44433
# 1 4    4.14577     3.39908
# 1 5    4.10052     3.33398
# 1 6    4.03542     3.23494
# 1 7    3.93638     3.06657
# 1 8    3.768       2.7132
# 1 9    3.41464     2.68725
# 1 10   3.38868     2.65595
# 1 11   3.35739     2.61785
# 1 12   3.31929     2.57104
# 1 13   3.27248     2.51368
# 1 14   3.21512     2.44892
# 1 15   3.15036     2.42042
# 1 16   3.12186     3.18072
# 1 17   3.09915     3.14307
# 1 18   3.06151     3.07847
# 1 19   2.99691     2.95481
# 1 20   2.87325     2.6426
# 1 21   2.56104     1.57395
# 1 22   1.49239     0.198152
# 1 23   0.11659         .
# 2 0        .       0.0774493
# 2 1    0.0929392   0.0628265
# 2 2    0.0783164   0.0530541
# 2 3    0.068544    0.0466411
# 2 4    0.062131    0.0425651
# 2 5    0.0580549   0.0401519
# 2 6    0.0556418   0.0389878
# 2 7    0.0544776   0.0388564
# 2 8    0.0543463   0.0395645
# 2 9    0.0550543   0.0652562
# 2 10   0.0807461   0.0961335
# 2 11   0.111623    0.133588
# 2 12   0.149078    0.179509
# 2 13   0.194999    0.236528
# 2 14   0.252018    0.308598
# 2 15   0.324088    0.403863
# 2 16   0.419353    0.642885
# 2 17   0.441375    0.67892
# 2 18   0.47741     0.738114
# 2 19   0.536604    0.835032
# 2 20   0.633522    0.983778
# 2 21   0.782268    1.06846
# 2 22   0.866948    0.489562
# 2 23   0.288052        .
# 3 0        .       0.0653605
# 3 1    0.0784326   0.0854854
# 3 2    0.0985575   0.108843
# 3 3    0.121915    0.137101
# 3 4    0.150173    0.173047
# 3 5    0.186119    0.22177
# 3 6    0.234842    0.293719
# 3 7    0.306791    0.414432
# 3 8    0.427504    0.666923
# 3 9    0.679996    2.66324
# 3 10   0.676308    2.65888
# 3 11   0.671949    2.65369
# 3 12   0.66676     2.64739
# 3 13   0.660466    2.63904
# 3 14   0.652109    2.62293
# 3 15   0.635998    2.56073
# 3 16   0.573804    2.5579
# 3 17   0.570973    2.55362
# 3 18   0.566687    2.54856
# 3 19   0.56163     2.55309
# 3 20   0.566165    2.64773
# 3 21   0.660806    3.33975
# 3 22   1.35282     4.82718
# 3 23   2.84025         .
# ;
#
# T [*] :=
 # 1 329.204    5 329.517    9 333.492   13 333.532   17 332.773   21 335.521
 # 2 329.265    6 329.683   10 333.496   14 333.549   18 332.804   22 349.687
 # 3 329.33     7 329.987   11 333.503   15 333.483   19 332.889
 # 4 329.409    8 330.697   12 333.515   16 332.759   20 333.251
# ;
#
# Q = 30.8373
#
# END OF FILE
