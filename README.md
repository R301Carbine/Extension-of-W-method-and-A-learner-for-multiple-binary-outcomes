# Extension-of-W-method-and-A-learner-for-multiple-binary-outcomes
The file "MainFunctions.R" contains the implementation of the proposed methods: RRR_WB, based on the W-method, and RRR_AB, derived from the A-learner approach. The arguments for these methods are documented within MainFunctions.R, and the outputs of the methods are as follows:
1. obj: A matrix storing the objective function values recorded at each parameter update (rows represent iterations, and columns represent different initializations).
2. V: An estimated $m \times r$ orthogonal matrix representing coefficients from multiple outcomes to latent variables.
3. W: An estimated $p\times r$ matrix representing coefficients from explanatory variables to latent variables.
Heterogeneous treatment effects are represented as $XWV^\prime$, where $r=2$, subgroup structures can be visually interpreted using biplots of $W$ and $V$.

The estimated matrices for explanatory variables and outcomes, referred to as $W$ and $V$, can be visualized using the "Biplot" function.
This function generates a graphical representation that helps interpret the relationships between explanatory variables and multiple outcomes. It is important to ensure that the rows and columns of these matrices are properly labeled before using the function.

In the biplot:

Positive Relationship: If an explanatory variable forms an angle of less than 90 degrees with the vector of a specific outcome, this indicates a positive relationship. In other words, as the value of the explanatory variable increases, the treatment effect on that outcome becomes larger.

Negative Relationship: If an explanatory variable forms an angle greater than 90 degrees with the vector of a specific outcome, this indicates a negative relationship. In this case, as the value of the explanatory variable increases, the treatment effect on that outcome decreases.

This type of visualization is particularly useful for understanding the variability in treatment effects across different subgroups and outcomes. It provides an intuitive way to explore how explanatory variables are associated with specific outcomes.
