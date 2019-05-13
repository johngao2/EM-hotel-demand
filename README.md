**Read info.pdf for a description of the project**
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

**Step 3 sprint 1 goals:**
- Data processing:
  - Preprocess original transaction data. **DONE**
  - Build transaction vector: all bookings at beginning of day, each row is a (subperiod of day (periods/day = max bookings in 1 day), arrive date, intraday period) tuple, each cell is a room type. **DONE**
    - make sure to collapse groups **DONE**
  - Build availability table: each row is possible (book date, arrive date) tuple, each col is roomtype. **DONE**
    - don't collapse groups, but add/subtract capacity based on partial fillings **DONE**
    - assume an order placed on a certain day affects the capacity for all of that day's intraday periods **DONE**
    - need to confirm that it matches num of rows with trans vec **DONE**
  - Customer types: Room type dependent only, constains 2, 3, and 4 tuples, cheap to expensive and expensive to cheap,  grouped by similarity. **DONE**
- **Log likelihood at end of sprint: -54573.9**
 
  
**Step 3 sprint 1 goals:**
- Data processing:
  - Group arrivals into weeks
    - Gotta deal with availability downwards bias somehow
  - For avail: assume the first time a room goes negative 1 that the room capacity is increased by 1 (to account for rm 199 and rm 200)
  - Investigate lambda scaling for nonbusy dates, see if it helps
- Model:
  - Add lasso regularization to ipopt code.
  - Code up VR AIC and RMSE metrics for testing
  
## Other notes:
**Sprint 1 Q&A:**
- Group arrive dates into weekly? If so, how to deal with availability, since it'll be biased down?
- Is there a better way to deal with intraday stuff?
  - scale num periods by lambda for sparse dates
- How will sparsity affect performance?
  - want hessians to have few nonzero values
- How will abundant availability affect performance?
  - assume from the first day a value goes negative a room became available
- Evaluation metric? (VR uses AIC and RMSE)
  - this is fine
 
**Sprint 1 Assumptions:**
- No longer looking at arrival and depart dates, but just arrival (can be tightened to weekends only, later)
  - This means both trans and avail datasets will focus on arrival time
  - Of course this still makes the data super sparse, so bad performance is expected in this sprit
- Availability assessed partially, e.g. if 1 person in a 4 person room cancels then add 1/4 to capacity
- Assume all rooms are available at start of year (e.g. nobody booked 2 years in advance)

**misc**
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
- transactions:
  - try putting intraday transactions both at start and end of day

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

