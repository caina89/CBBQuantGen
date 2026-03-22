## Lambda estimation using REML 
REML uses likelihood function to estimate the Genetic Variance ($\sigma_g^2$) and the Residual Variance ($\sigma_e^2$). Their ratio then defines the $\lambda$ for BLUP/Ridge Regression. 
To demonstrate this we use the p3d or lme4 R package, adapted for genomic data in the R package rrBLUP, as shown in `reml_lambda.R`. 

In this script, the mixed.solve function performs Restricted Maximum Likelihood (REML). It looks at the distribution of the phenotypes $y$ and asks: "Given the relationship structure of these 10 SNPs, what values of $\sigma_g^2$ and $\sigma_e^2$ make this specific pattern of $y$ most likely?" It ignores the fixed effects (the mean) to ensure the variance estimates aren't biased.

### The Relationship to $\lambda$
* If the data is very noisy, REML will estimate a large $\sigma_e^2$.
* A large $\lambda$ automatically "shrinks" your marker effects in the BLUP solution.

### Merits of using REML 
REML is more stable because it uses the mathematical properties of the entire dataset at once to find the "biological" $\lambda$. As long as your assumption of a Normal distribution for SNP effects is mostly true, the REML $\lambda$ will be very close to the "optimal" prediction $\lambda$.

## Elastic net 
In Elastic Net the penalty term is: 

$\text{Penalty} = \lambda\left[\alpha\underbrace{\sum|\beta|}_{L_1\text{(LASSO)}} + (1 - \alpha)\underbrace{\frac{1}{2}\sum\beta^2}{L_2\text{ (Ridge)}}\right]$
* When $\alpha = 1$: The $L_2$ part becomes $0$. You have a Pure LASSO model.
* When $\alpha = 0$: The $L_1$ part becomes $0$. You have a Pure Ridge model.
* When $\alpha = 0.5$: You are giving equal "weight" to both penalties.
  
To compare coefficients ($\beta$) given to SNPs in Ridge Regression, LASSO and Elastic Net, we can use the script shown in `compare_coefficients.R` which uses the glmnet, a package for penalized regression in R.

When you run the script `compare_coefficients.R` you will see that: 

Ridge ($\alpha = 0$)
* Result: $[0.21, 0.19, 0.20, 0.18, \dots]$
* Behavior: Because the SNPs are in LD, Ridge refuses to pick a favorite. It sees 10 SNPs that all look like the causal variant and spreads the credit almost equally among them. This captures the heritability of the region but fails to identify the specific causal SNP.

LASSO ($\alpha = 1$)
* Result: $[1.85, 0, 0, 0, \dots]$
* Behavior: LASSO sees the redundancy and "fires" 9 of the SNPs. It picks the one that fits the data slightly better and gives it a large weight, setting the rest to exactly zero. This handles the kurtosis (sparsity) perfectly but can be unstable if that one SNP has a measurement error.

Elastic Net ($\alpha = 0.5$)
* Result: $[0.85, 0.42, 0.38, 0, 0, \dots]$
* Behavior: This is the compromise. It zeroes out the SNPs that are clearly noise (Selection), but for the SNPs that are highly correlated and close to the signal, it groups them together (Grouping). It keeps a few strong candidates in the model, providing a more stable prediction than LASSO.
