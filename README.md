## Current Objective (subject to change):
Run simple logit model from scratch to forecast customer preferences given product assortments at each time step, **and do it quick**. Transaction data is simulated. Hoping to do the same with Multinomial Logit and Expectation-Maximization in the future.

### Tools:
IPOPT: Fast open source non-linear programming solver. https://projects.coin-or.org/Ipopt, https://github.com/coin-or/Ipopt

CPPAD: Fast open source package for computing derivatives. Needed since otherwise I'd have to code up gradients and hessians to feed into IPOPT. https://coin-or.github.io/CppAD/doc/cppad.htm, https://github.com/coin-or/CppAD

