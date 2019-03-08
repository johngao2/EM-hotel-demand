## Project Overview:
~~Step 1: Implement simple multinomial logit go get familiar with likelihood maximization using CppAD and IPOPT packages.~~ **- Done**  

Step 2: Replicate EM algorithm (van Ryzin, Vulcano 2017) on transaction data to forecast customer preferences given product assortments at each time step, **and do it quick**. Use on simulated and real datasets from the original paper, and compare results.  

Step 3: Extend van Ryzin & Vulcano to include seasonality and be fast enough to work on hotel problems where the assortment size is gigantic.

## Background:  
**Newer, problem-specific EM algorithm (what I'm trying to reproduce and extend):**  
van Ryzin, G., & Vulcano, G. (2017). "Technical Note: an Expectation-Maximization Method to
Estimate a Rank-based Choice Model of Demand", *Operations Research 65 (2017), no. 2*, 396–407. Retrieved from https://pubsonline.informs.org/doi/pdf/10.1287/opre.2016.1559

Supplemental material: https://pubsonline.informs.org/doi/suppl/10.1287/opre.2016.1559  
(Test data also comes from here)

**Original EM paper:**  
Dempster, A., Laird, N., & Rubin, D. (1977). "Maximum Likelihood from Incomplete Data via the EM Algorithm". *Journal of the Royal Statistical Society. Series B (Methodological), 39(1)*, 1-38. Retrieved from http://www.jstor.org/stable/2984875

**Helpful tutorials:**  
S. Borman. (2004). "The Expectation Maximization Algorithm: a Short Tutorial". Manuscript. Retrieved from http://www.seanborman.com/publications/EM_algorithm.pdf  

J. Bilmes. (1997). "A Gentle Tutorial on the EM Algorithm and its Application to Parameter Estimation for
Gaussian Mixture and Hidden Markov Models". *Technical Report ICSI-TR-97-02* University of
Berkeley. Retrieved from http://melodi.ee.washington.edu/people/bilmes/mypapers/em.pdf

**Tools:**  
IPOPT: Fast open source non-linear programming solver. https://projects.coin-or.org/Ipopt, https://github.com/coin-or/Ipopt  

CPPAD: Fast open source package for computing derivatives. Needed since otherwise I'd have to code up gradients and hessians and whatnot D: to feed into IPOPT. https://coin-or.github.io/CppAD/doc/cppad.htm, https://github.com/coin-or/CppAD

