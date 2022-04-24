# EquilibriumScaling
To run the program first Gurobi C++ API should be ready 

Scaling.cpp:
 - Select .lp file to feed the program with your LP/IP or MIP and it gives you Row and Column Scaling Factors 
 - Optionally You can reach scaled A, B and C matrices, for that you need to uncomment certain part of the code. For printing A, B, C, printing to console is not recommended since it makes program slowed down, so it is better to write them to the file  
 - Use can use the scaling factors (R and C) or A, B and C for your own needs, for example, they can be feed into the program again and make program solve scaled program
