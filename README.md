## Project Log:
#### Step 1: Implement simple multinomial logit go get familiar with likelihood maximization using CppAD and IPOPT packages. *(MNL_baseline > simple_mnl.cpp)*  **DONE**  
**Subtasks:**
- Run IPOPT solver on single problem. **DONE**
- Run IPOPT solver repeatedly on changing problem formulations. **DONE** 

#### Step 2: Replicate EM algorithm (van Ryzin, Vulcano 2017) on transaction data to forecast customer preferences given product assortments at each time step, **and do it quick**. Use on simulated and real datasets from the original paper, and compare results. *(censored_EM_replication > hotel5 > EM_hotel_5.cpp)* **DONE**
**Subtasks:** 
- Implement Uncensored EM on van Ryzin's hotel data, compare results. **DONE, BUT WITH ISSUES**
  - Check with closed form maximizer. **DONE; MAXIMIZER ISN'T THE PROBLEM**
  - If bug still isn't solved, run through algo on small by hand and track every number, then compare with algo output. **DONE; RESULTS REPLICATED, BUG SOLVED**
  - Debug vanishing taylor coefficient (scaling probably works) **DONE; NOT REALLY A PROBLEM ANYMORE; HAPPENS WAY AFTER CONVERGENCE**
- Implement censored EM on van Ryzin's simulated data, compare results.
  - Need to first simulate the data myself, following steps outlined in 5.1 (van Ryzin, Vulcano 2017) **DONE**
  - Implement censored algo on my own simulated data **DONE, SEEMS TO RECOVER PARAMS WELL**

#### Step 3: Extend van Ryzin & Vulcano to include seasonality and be fast enough to work on hotel problems where the assortment size is gigantic.
**Subtasks:**
- Create availability table for Cabot data, where each row is either a booking/cancellation, and each col is a possible stay date/room type combination, cell value is number of rooms available
- Create transaction vector, where each row is a subperiod of day (periods/day = max bookings in 1 day):
  - One version will have all bookings at beginning of day, other will have all bookings at end of day
- Create customer types, where each type is a list of possible product choices, including stay dates and room types. Orderings to be decided later

**Sprint 1 goals:**
- Input data:
  - Preprocess original transaction data. **DONE**
  - Build transaction vector: all bookings at beginning of day, each cell is date range and room type. **DONE**
    - simplify by only including room type info in each cell.
  - Build availability table: each row is possible book date, each col is (roomtype, stay date) tuple.
    - simplify by making each col only roomtype.
  - Customer types: Room type dependent only, constains 2, 3, and 4 tuples, cheap to expensive and expensive to cheap, grouped by similarity.
- Model:
  - Add lasso regularization to ipopt code.
  
### Other notes:
- Review R GLM for poission regression, what kind of link functions are reasonable?
  - Think of variables we need to feed into GLM to model seasonality (only need to worry about time related stuff; one seasonality for booking, one for arrival)
- Lambda will depend on when arrival occurs and booking
- Try lasso penalty in ipopt
  - min f(x) + alpha(sum(y)) where -y_i <= x_i <= y_i, use this since solvers HATE abs
- booking:
  - for now, constant lambda during 2 week period with weekly seasonality pattern
  - use multiplicative to preserve signs
- arrival:
  - Start with 2 weekend pair (which can shift) + something with a room types
  - each customer type consists of a pair of weekends (e.g. one weekend vs weekend after)
  - most important thing to consider is type; weekend pair used for tiebreakers

## References:  
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

CPPAD: Fast open source package for computing derivatives. https://coin-or.github.io/CppAD/doc/cppad.htm, https://github.com/coin-or/CppAD

