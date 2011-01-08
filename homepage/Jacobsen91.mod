#################### DESCRIPTION ###############################################
#
# Last updated: 22 Aug 2009
#
# The Jacobsen91.mod has 5 solutions. They correspond to the steady states of 
# the methanol-propanol column (mass reflux, energy-balances included) discussed
# in:
#
# E. W. Jacobsen, S. Skogestad;
# Multiple Steady States in Ideal Two-Product Distillation;
# AIChE Journal, 1991, 37, 499-511.
#
# One solution is infeasible in practice (negative flow rates). Mass reflux Lw 
# can be chosen as bifurcation parameter. yD and boilup V in the above paper are
# y[1] and V_N, respectively, in the AMPL model. Bifurcation diagrams are given 
# in Figure 8 of that paper.
#
# Solving this problem with the interval methods is discussed in:
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

# The first detailed report on applying interval methods to this problem seems 
# to be:
#
# A. Baharev;
# Application of Interval Methods to Chemical Engineering Problems;
# PhD dissertation; in Hungarian, Budapest University of Technology and 
# Economics, Department of Chemical and Environmental Process Engineering, 2009
#
# High resolution bifurcation diagrams are presented on page 72 of the 
# dissertation. The dissertation is available from:
#
# http://reliablecomputing.eu
#
################## PROBLEM STATISTICS ##########################################
#
# Presolve eliminates 0 constraints and 1 variable.
# Substitution eliminates 24 variables.
# Adjusted problem:
# 19 variables:
#         16 nonlinear variables
#         3 linear variables
# 19 constraints, all nonlinear; 90 nonzeros
# 0 objectives.
#
################# PARAMETERS FROM SPECIFICATION ################################

param Lw  := 96.0;
param V_N := 3.0;

param N   := 8;
param N_F := 5;

param  F{1..N};
param  z{1..N};
param qF{1..N};

param alpha := 3.55;
param M1 := 32.04;
param M2 := 60.10;

param V_L{1..N}; param V_U{1..N};

############### VARIABLES ######################################################

var x{1..N} >= 1.0e-4, <= 1.0;
var V{j in 1..N} >= V_L[j], <= V_U[j];

var D >= 0.0, <= 1.12;
var Q >= 0.0, <= 2.30;
var d >= 0.0, <= 0.52;
var q >= 0.0, <= 0.187;

############# DEFINED VARIABLES (THEY ARE ELIMINATED BY SUBTITUTION) ###########

var y{j in 1..N} = (alpha*x[j])/(1.0+(alpha-1.0)*x[j]);

var HL{j in 0..N-1} =
  if j == 0 then
    0.1667*exp(-1.087*y[1])
  else 
    0.1667*exp(-1.087*x[j]);

var HV{j in 1..N} =
    0.1349*exp(-3.98*x[j])+0.4397*exp(-0.088*x[j]);

############## EQUATIONS #######################################################

CS_d:
  d = D*y[1];

CS_Q:
  Q = V[1]*(HV[1]-HL[0]);

CS_q:
  q = D*HL[0];

M_eq{j in 1..N-1}:
  sum{k in 1..j}F[k]*z[k] + V[j+1]*y[j+1] = d + (sum{k in 1..j}F[k]+V[j+1]-D)*x[j];

M_tot:
  F[N_F]*z[N_F] = d + (F[N_F]-D)*x[N];

H_eq{j in 1..N-1}:
  sum{k in 1..j} qF[k] + V[j+1]*HV[j+1] = Q + q + (sum{k in 1..j}F[k]+V[j+1]-D)*HL[j];

spec_Lw:
  Lw = (V[1]-D)*(M1*y[1]+M2*(1.0-y[1]));

################### DATA SECTION ###############################################

for {j in 1..N} {
  let   F[j] := 0.0;
  let   z[j] := 0.0;
  let  qF[j] := 0.0;
  let V_L[j] := 2.0;
  let V_U[j] := 4.0;
}

let V_L[N] := V_N;
let V_U[N] := V_N;

let  F[N_F] := 1.0;
let  z[N_F] := 0.5;
let qF[N_F] := F[N_F]*0.1667*exp(-1.087*z[N_F]);


################################################################################

data;

######### SOLUTION 1 ###########################################################
# 
# var x :=
# 1  0.995726
# 2  0.986317
# 3  0.958235
# 4  0.880923
# 5  0.708062
# 6  0.545244
# 7  0.315806
# 8  0.132220;
# 
# var V :=
# 1  3.4175
# 2  3.41489
# 3  3.40704
# 4  3.38495
# 5  3.33062
# 6  3.25931
# 7  3.13778
# 8  3.00000;
# 
# var D = 0.424407;
# 
# var Q = 1.193;
# 
# var d = 0.423895;
# 
# var q = 0.0238897;
# 
# 
######### SOLUTION 2 ###########################################################
# 
# var x :=
# 1  0.196794
# 2  0.103881
# 3  0.0783623
# 4  0.0717959
# 5  0.0701344
# 6  0.020268
# 7  0.00565783
# 8  0.00157120;
# 
# var V :=
# 1  3.11562
# 2  3.0452
# 3  3.03525
# 4  3.03396
# 5  3.03374
# 6  2.99562
# 7  2.99835
# 8  3.00000;
# 
# var D = 1.07511;
# 
# var Q = 1.22522;
# 
# var d = 0.500118;
# 
# var q = 0.10809;
# 
# 
######### SOLUTION 3 ###########################################################
# 
# var x :=
# 1  0.997599
# 2  0.99207
# 3  0.97462
# 4  0.922358
# 5  0.787479
# 6  0.710547
# 7  0.537422
# 8  0.292755;
# 
# var V :=
# 1  3.28779
# 2  3.28632
# 3  3.28165
# 4  3.26747
# 5  3.22876
# 6  3.19044
# 7  3.13072
# 8  3.00000;
# 
# var D = 0.293312;
# 
# var Q = 1.14755;
# 
# var d = 0.293114;
# 
# var q = 0.0165009;
# 
# 
######### SOLUTION 4 ###########################################################
# 
# var x :=
# 1  0.934517
# 2  0.81886
# 3  0.602621
# 4  0.35567
# 5  0.192086
# 6  0.074504
# 7  0.0256403
# 8  0.00819998;
# 
# var V :=
# 1  3.45205
# 2  3.41756
# 3  3.34137
# 4  3.21336
# 5  3.09752
# 6  3.01186
# 7  2.9975
# 8  3.00000;
# 
# var D = 0.505736;
# 
# var Q = 1.21114;
# 
# var d = 0.495947;
# 
# var q = 0.0290348;
# 
# 
######### SOLUTION 5 ###########################################################
# 
# var x :=
# 1  0.98717
# 2  0.960178
# 3  0.886708
# 4  0.722446
# 5  0.475227
# 6  0.254145
# 7  0.103815
# 8  0.0355046;
# 
# var V :=
# 1  3.47014
# 2  3.46248
# 3  3.44119
# 4  3.3893
# 5  3.28595
# 6  3.1457
# 7  3.03157
# 8  3.00000;
# 
# var D = 0.483422;
# 
# var Q = 1.21222;
# 
# var d = 0.481659;
# 
# var q = 0.0272839;
# 
# 
################################################################################


option show_stats 1;
option presolve 10;
option substout 1;
option var_bounds 2;
option nl_comments 0;

# option solver ipopt;

solve;

# display x;
# display V;
# display D;
# display Q;
# display d;
# display q;
#
# END OF FILE
