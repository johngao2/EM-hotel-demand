# declare the variables and their bounds, 
# set notation could be used, but this is straightforward
var x1 >= 1, <= 5; 
var x2 >= 1, <= 5; 
var x3 >= 1, <= 5; 
var x4 >= 1, <= 5;

# specify the objective function
minimize obj:
    x1 * x4 * (x1 + x2 + x3) + x3;
        
# specify the constraints
s.t.
  inequality:
    x1 * x2 * x3 * x4 >= 25;
                
  equality:
    x1^2 + x2^2 + x3^2 +x4^2 = 40;

# specify the starting point            
let x1 := 1;
let x2 := 5;
let x3 := 5;
let x4 := 1;

# solve the problem
solve;

# print the solution
display x1;
display x2;
display x3;
display x4;