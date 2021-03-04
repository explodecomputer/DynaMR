# DynaMR

# Starting parameters

Inputs:
Genotype (V)
Random environmental noise (V)
Variances of total levels (C) 
Means of total levels (C) 
Sample size (C) 

Outputs:
Total levels (V)


# Dynamic model

(Standard version)
Inputs:
Total levels (V)
Kinetic parameters (C) 
Time scale (C) 

Outputs:
Levels over time (V)
Steady state values (V)

(Perturbation version)
Inputs:
Total levels (V)
Kinetic parameters (C) 
Time scale (C) 
Modification time (C)
Modification magnitude (C)

Outputs:
Levels over time (V)
Steady state values (V)


# Disease model

Inputs:
Effects on disease (C)
Steady state values (V)
Other genetics (V)
Other noise (V)

Outputs:
Disease state (V)


# MR analysis

Inputs:
Genotypes (V)
Disease state (V)
Protein levels at timepoint t (V)
Outputs:
MR estimates


# True causal effect

Inputs:
Perturbation effect (C)
Disease state (V)

Outputs:
Effects of protein perturbations on disease