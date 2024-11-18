# Extension-of-W-method-and-A-learner-for-multiple-binary-outcomes
The file MainFunctions.R contains the implementation of the proposed methods: RRR_WB, based on the W-method, and RRR_AB, derived from the A-learner approach. The arguments for these methods are documented within MainFunctions.R, and the outputs of the methods are as follows:
1. obj: A matrix storing the objective function values recorded at each parameter update (rows represent iterations, and columns represent different initializations).
2. V: An estimated $m \times r$ orthogonal matrix representing coefficients from multiple outcomes to latent variables.
3. W: An estimated $p\times r$ matrix representing coefficients from explanatory variables to latent variables.
Heterogeneous treatment effects are represented as $XWV^\prime$, where $r=2$, subgroup structures can be visually interpreted using biplots of $W$ and $V$.
