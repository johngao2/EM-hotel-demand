## Current Objective (subject to change):
Implement EM algorithm (van Ryzin, Vulcano 2017) to forecast customer preferences given product assortments at each time step, **and do it quick**. Transaction data is simulated.

## Background:
**Original EM paper:**  
Dempster, A., Laird, N., & Rubin, D. (1977). Maximum Likelihood from Incomplete Data via the EM Algorithm. *Journal of the Royal Statistical Society. Series B (Methodological), 39(1)*, 1-38. Retrieved from http://www.jstor.org/stable/2984875

**Newer, problem-specific EM algorithm (what I'm trying to reproduce and extend):**  
van Ryzin, G., & Vulcano, G. (2017). Technical note: an expectation-maximization method to
estimate a rank-based choice model of demand, *Operations Research 65 (2017), no. 2*, 396–407. Retrieved from https://pubsonline.informs.org/doi/pdf/10.1287/opre.2016.1559

**Helpful tutorials:**  
S. Borman. (2004). The Expectation Maximization Algorithm: a Short Tutorial. Manuscript. Retrieved from http://www.seanborman.com/publications/EM_algorithm.pdf  
J. Bilmes. (1997). A Gentle Tutorial on the EM Algorithm and its Application to Parameter Estimation for
Gaussian Mixture and Hidden Markov Models. *Technical Report ICSI-TR-97-02* University of
Berkeley. Retrieved from http://melodi.ee.washington.edu/people/bilmes/mypapers/em.pdf

## Tools:
IPOPT: Fast open source non-linear programming solver. https://projects.coin-or.org/Ipopt, https://github.com/coin-or/Ipopt

CPPAD: Fast open source package for computing derivatives. Needed since otherwise I'd have to code up gradients and hessians and whatnot D: to feed into IPOPT. https://coin-or.github.io/CppAD/doc/cppad.htm, https://github.com/coin-or/CppAD

