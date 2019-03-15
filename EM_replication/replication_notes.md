### Log-Likelihood results (buy-up model):

| Hotel | My result | Their Result | Diff (num) |
|-------|-----------|--------------|------------|
| 1     | -2427.35  | -2194        | 233.35     |
| 2     | -440.561  | -318         | 122.561    |
| 3     | -1576.22  | -1376        | 200.22     |
| 4     | -519.77   | -375         | 144.77     |
| 5     | -493.95   | -371         | 122.95     |

Possible troubleshooting tips:
- For taylor poly error: possible some value is too big or too small, look at value of solutions when this crashes
- First try to test with closed form x_i, confirm number
- Also try calculating small handmade dataset by hand and confirming that numbers match when running through the algo
