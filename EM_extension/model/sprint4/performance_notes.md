**HOTEL5**: 
    - params: hotel5 dataset, alpha 0.1, stopping criteria 1e-4, 
    - Benchmark:
        - LL: -394.27
    - Closed form:
        - LL: -213.30

**SIM DATA**
    - dataset: 10k times, 15 products, 10 types, sample #1
    - params: alpha 0.1, stopping criteria 1e-4, 
    - True values:
        - True x_vec: [0.148792348, 0.231331468, 0.085121344, 0.004869606, 0.142366357, 0.05919758, 0.070649566, 0.022415964, 0.203436424, 0.031819343]
        - True lambda: 0.8
    - Benchmark:
        - Benchmark LL: -217817.27
        - Benchmark x_vec: [0.13766, 0.17220, 0.08931, 0.09812, 0.07734, 0.09444, 0.07908, 0.05426, 0.18230, 0.01527]
        - Benchmark lambda: 0.73732
        - Benchmark RMSE: 0.04626
    - Closed form Xi's and lambda:
        - LL: -126888.76
        - x_vec: [0.15428, 0.23251, 0.08725, 0.00759, 0.14325, 0.06021, 0.07410, 0.02171, 0.20693, 0.01215]
        - lambda: 0.78511
        - RMSE: 0.007853
    - Closed form Xi, ipopt lambda:
        - LL: -210467.46
        - x_vec: [0.15212, 0.22927, 0.08604, 0.00747, 0.14126, 0.05937, 0.07307, 0.02140, 0.20404,
        0.02594]
        - lambda: 
        - RMSE: 0.008960


**SPRINT 1 DATA**
    - dataset: 3.5m times, 7 products, 8 types
    - params: alpha 0.1, stopping criteria 1e-3, 
    - Benchmark:
        - LL: -54573.92
    - Closed form:
        - LL: -16405.24

**Independent model:**
    - dataset: 24k times, 4900 products, 4900 types
    - params: alpha 0.1, stopping criteria 1e-9, 
    - Benchmark:
        - too slow, unrunnable
    - Closed form:
        - LL: -50089.51
    - Closed form + weekly lambdas:
        - LL -59845.03
        - lambdas: [0.28830, 0.30427, 0.27572, 0.27369, 0.24937, 0.18999, 0.16983]