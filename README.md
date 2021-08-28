# DynaMR

## Paper outline (note the model numbering here is different to the coding model numbering, see below)

- Model 1: Simple protein phosph/dephsph protein cycle
    - Protein level -> disease
    - MR vs RCT vs Obs
    - Instrument specificity
- Model 2: Protein complex of 3 activated proteins
    - Same questions as model 1
- Model 3: Feedback loop using Baker et al
    - Steiger filtering
    - Temporal effects
    - Instrument specificity
- Model 4: Disease states
    - Non-linearities
    - Genetic interactions
    - Additive genetic variance on disease

---

## Code model numbering (note the model numbering here is different to the paper model numbering, see above)

- Model 1: Simple protein phosph/dephsph protein cycle
    - X -> Y phosphorylated by K and dephosporylated by P
- Model 2: Baker et al model
    - Only the starting conditions for population are simulated (not the parameters)
    - Only dynamical analysis
- Model 3: Baker et al model
    - Both starting conditions and parameters for population are simulated
    - Only dynamical analysis
- Model 4: Baker et al model
    - Both starting conditions and parameters for population are simulated
    - Includes steady-state analysis as in Model 1
- Model 5: Protein complex of 2 activated proteins and complexes, STAT system example
    - The code to run dynamics simulation is /scripts/dynamics_model5_test_run.r

---

## Running

Create a `config.json` file like:

```json
{
    "resultsdir": "/path/to/results/dir"
}
```

Then run `snakemake`

```
module add languages/anaconda3/5.2.0-tflow-1.11
source ~/.bash_profile
snakemake 

```

---

## Starting parameters

Inputs:

- Genotype (V)
- Random environmental noise (V)
- Variances of total levels (C) 
- Means of total levels (C) 
- Sample size (C) 

Outputs:

- Total levels (V)


## Dynamic model

### Standard version

Inputs:

- Total levels (V)
- Kinetic parameters (C) 
- Time scale (C) 

Outputs:

- Levels over time (V)
- Steady state values (V)

### Perturbation version

Inputs:

- Total levels (V)
- Kinetic parameters (C) 
- Time scale (C) 
- Modification time (C)
- Modification magnitude (C)

Outputs:

- Levels over time (V)
- Steady state values (V)

# Disease model

Inputs:

- Effects on disease (C)
- Steady state values (V)
- Other genetics (V)
- Other noise (V)

Outputs:

- Disease state (V)


# MR analysis

Inputs:

- Genotypes (V)
- Disease state (V)
- Protein levels at timepoint t (V)

Outputs:

- MR estimates


# True causal effect

Inputs:

- Perturbation effect (C)
- Disease state (V)

Outputs:

- Effects of protein perturbations on disease

