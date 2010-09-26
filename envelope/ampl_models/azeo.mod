
################## PARAMETERS ##################################################

param C := 3;  # NUMBER OF COMPONENTS

# VALUES ARE GIVEN IN THE DATA SECTION

param a{1..C};
param b{1..C};
param c{1..C};
param V{1..C};
param k{1..C,1..C};

param ln_P := log(760.0);


# LOWER/UPPER BOUNDS ON THE VARIABLES

param x_L{1..C}; param x_U{1..C};

# INITIAL ESTIMATES OF THE VARIABLES (IF NEEDED)

param    x_0{1..C};
param ln_K_0{1..C};
param    T_0;

############### VARIABLES ######################################################

var x{i in 1..C} >= x_L[i], <= x_U[i], := x_0[i];

var ln_K{i in 1..C} := ln_K_0[i];

var T >= 330.0, <= 380.0, := T_0;

####### DEFINED VARIABLES (THEY ARE ELIMINATED BY PRESOLVE / SUBTITUTION) ######

var ln_pr{i in 1..C} = a[i]-b[i]/(T+c[i]) - ln_P;

var rT = (1.0/1.9858775)/T;

var Lambda{i in 1..C, j in 1..C} = V[j]/V[i]*exp(-k[i,j]*rT);

var s{i in 1..C} = sum{j in 1..C} (x[j]*Lambda[i,j]);

var t{i in 1..C} = x[i]/s[i];

var u{i in 1..C} = 1.0 - sum{j in 1..C} Lambda[j,i]*t[j];

var ln_gamma{i in 1..C} = -log(s[i]) + u[i];
  
############## EQUATIONS #######################################################

# AUXILIARY EQUATIONS

E_aux_K{i in 1..C}:
  ln_K[i] = ln_gamma[i]+ln_pr[i];

# SUMMATION EQUATIONS

S_x_eq:
  sum{i in 1..C} x[i] - 1.0 = 0.0;

# AZEOTROPE CONDITION

y_equals_x{i in 1..C}:
  x[i]*ln_K[i] = 0.0;

################### DATA SECTION ###############################################

data;

let a[1] := 18.679031;
let a[2] := 16.264448;
let a[3] := 18.584878;

let b[1] := 3667.704902;
let b[2] := 2904.342681;
let b[3] := 3984.922839;

let c[1] := -46.976;
let c[2] := -51.191;
let c[3] := -39.734;

let V[1] := 58.39;
let V[2] := 89.57;
let V[3] := 18.05;

let k[1,1] := 0.0;
let k[1,2] := 694.0825;
let k[1,3] := 393.1971;

let k[2,1] :=-149.7978;
let k[2,2] := 0.0;
let k[2,3] := 6811.3433;

let k[3,1] := 926.2630;
let k[3,2] := 1888.8509;
let k[3,3] := 0.0;

# LOWER AND UPPER BOUNDS ON THE VARIABLES

let x_L[1] := 0.00;
let x_U[1] := 1.00;

let x_L[2] := 0.00;
let x_U[2] := 1.00;

let x_L[3] := 0.00;
let x_U[3] := 1.00;

################################################################################

# INITIAL ESTIMATES (IF NEEDED)

let x_0[1] := 0.333;
let x_0[2] := 0.333;
let x_0[3] := 0.333;

let T_0 := 345.0;

for {i in 1..C}
	let ln_K_0[i] := 0.0;


################################################################################

option show_stats 1;
option presolve 10;
option substout 1;
option var_bounds 2;
#option nl_comments 1;
#option nl_permute 0;
#option auxfiles cf;
#write gazeo;

option solver "/home/ali/ampl/ipopt";
solve;

display x, T;

#print "";
#solexpand {i in 1.._nvars:  _var[i].status == 'sub'} _con[_var[i].defeqn];

#
# END OF FILE
