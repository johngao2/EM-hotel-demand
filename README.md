# README

## Objective
Run logit model from scratch to forecast customer preferences given product assortments at each time step, **and do it quick**. Transaction data will be drawn/simulated based on a known distribution with reasonable parameters, or just taken from another paper. Then I'll run the same data through an EM algorithm, and compare.

## Baseline
sklearn's logistic regression will be used a baseline. The main metric is time, but accuracy will also be compared. Accuracy really shouldn't differ though since it's literally the same model.

## Methodology

### Tools:
IPOPT: A fast, open source non-linear programming solver. https://projects.coin-or.org/Ipopt, https://github.com/coin-or/Ipopt

CPPAD: A fast, open source package for computing derivatives. https://coin-or.github.io/CppAD/doc/cppad.htm, https://github.com/coin-or/CppAD

### Data Format
wip

### The Math: Simple Logit
wip

## The Math: EM Algorithm
wip
