**Read info.pdf for an (incomplete) description of the project**
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

**Step 3 sprint 1 goals:**
- Data processing:
  - Preprocess original transaction data. **DONE**
  - Build transaction vector: all bookings at beginning of day, each row is a (subperiod of day (periods/day = max bookings in 1 day), arrive date, intraday period) tuple, each cell is a room type. **DONE**
    - make sure to collapse groups **DONE**
  - Build availability table: each row is possible (book date, arrive date) tuple, each col is roomtype. **DONE**
    - don't collapse groups, but add/subtract capacity based on partial fillings **DONE**
    - assume an order placed on a certain day affects the capacity for all of that day's intraday periods **DONE**
    - need to confirm that it matches num of rows with trans vec **DONE**
  - Customer types: Room type dependent only, contains 2, 3, and 4 tuples, cheap to expensive and expensive to cheap,  grouped by similarity. **DONE**
- **Sprint 1 AIC: 109163** with stopping criteria: max change in x vector < 1e-3

**Step 3 sprint 2 goals:**
- Data processing:
  - cust types:
    - stay length, day of week arrival are independent, stay length is reduced to 4 values, with 4 being stay for 4 days or more **DONE**
    - customer chooses b/w UNIT, 2 consecutive weeks. Also, they choose between prioritizing week vs prioritizing unit. Price is included in weeks already, dont worry about it. **DONE**
  - arrival:
    - each col is a (unit, week, day of week, stay_len) tuple, rows can be reduced to just look dates **DONE,**
    - when collapsing the multiindexed cols into numbers, that's used as the product number for trans vec **DONE**
  - trans:
    - each product num represents a combination of (unit, week, day of week, stay len tuple) **DONE**
- Model:
  - Add lasso regularization to ipopt code. **DONE**
  - Code up VR AIC and RMSE metrics for testing **DONE AIC**

**Step 3 sprint 3 goals:**
- parse wide csv's (boost tokenizer is probably a good call), https://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c/1595366#1595366 **DONE**
- Add L2 reg **DONE**
- Debug mem leak with valgrind **Done**
- Rerun with independent all factors **4900 types, still too slow, still many ~0 gradient variables**
  - try again without weekday arrivals **still stuck on 2nd deriv checker**
- Edit model such that mu matrix is saved, since building it takes years **DONE**
- Rewrite model such that cust types use a different format than product id's
- use binary search to find lagrange multiplier such that x_i sums to 1 **DONE**
  - do this by subbing prospective values into the quadratic soln' on the board and summing all x until sum of x = 1, with a precision of 1e-8 **DONE**
- move sigma and mu matrix preprocessing to python **tabled**
- **NEW WAY OF OPTIMIZING LL**: **DONE**
  - use binary search to find lambda (the lagrange, not arrival rate) such that x_i sums to 1 **DONE**
    - do this by subbing prospective values into the quadratic soln' on the board and summing all x until sum of x = 1, with a precision of 1e-8 **DONE**
- use python to create sigma matrix (or preference matrix) **CURRENTLY TABLED**
- implement day of week lambdas **DONE WITH STATIC MAP**

**Step 3 sprint 4 goals:**
- implement AIC into C++ code for quick stats **DONE**
- implement day of week lambdas as a formula of coefficients and d, d=day
  - book_arrive_diffs are aggregated on each d using mean
  - times with non-transactions have book_arrive_diffs of max(diff_book_arrive) across the entire dataset
- lambda quality of fit check:
  - take sum of (lambda * x_i for all compatible types) for each t
  - aggregate those sums by week
  - take difference of actual week arrival count and the lambda*xi aggregate (this is our 'residual')
  - graph residual over time, see if patterns emerge
  - repeat this with different aggregations, e.g. by month, day of week, etc.
- test new lambdas with large scale non-independent choice assumptions
- add 1 to available rooms everytime avail goes into negatives (should only affect DD room type, of which 2 were added in 2018)
  - debug weird negative values

## Other notes:
**Backlog:**
- Scale num periods by lambda for sparse dates, might fix sparsity problem if implementable
- Implement RMSE from Van Ryzin paper
- try putting intraday transactions at end of day

**misc notes**
- want hessians to have few nonzero values
- lasso won't work due to log(x) in objective function

**Sprint 1 Assumptions:**
- No longer looking at arrival and depart dates, but just arrival (can be tightened to weekends only, later)
  - This means both trans and avail datasets will focus on arrival time
  - Of course this still makes the data super sparse, so bad performance is expected in this sprit
- Availability assessed partially, e.g. if 1 person in a 4 person room cancels then add 1/4 to capacity
- Assume all rooms are available at start of year (e.g. nobody booked 2 years in advance)

**Sprint 2 Assumptions:**
- Customers decide based on arrival day of week, stay length, unit type, week of arrival
  - all 4 factors independent of each other
  - dow and stay len are independent, unit type and week of arrival are ordered
    - unit type orders are same as sprint 1, week orders are consecutive weeks with sooner preference and later preference types

**Sprint 3 Notes:**
- We have no way of knowing capacity at beginning of 2018

**misc**
- Review R GLM for poission regression, what kind of link functions are reasonable?
  - Think of variables we need to feed into GLM to model seasonality (only need to worry about time related stuff; one seasonality for booking, one for arrival)
- Lambda will depend on when arrival occurs and booking
- booking:
  - for now, constant lambda during 2 week period with weekly seasonality pattern
  - use multiplicative to preserve signs
- arrival:
  - Start with 2 weekend pair (which can shift) + something with a room types
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

