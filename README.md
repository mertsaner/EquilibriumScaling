## Equilibrium Scaling

### Prerequisite
- To run the program firstly Gurobi C++ API should be ready. To install Gurobi Solver [click here.](https://www.gurobi.com/free-trial/ "Gurobi Solver Website")
- The program require Eigen library to read large .lp files. Hence you need to include Eigen library to the program. For Eigen library download instructions [click here.](https://eigen.tuxfamily.org/ "Eigen library download")

### Usage
- To configure the program according your needs, you need to edit `ScalingFactors.cpp` file.
- In `ScalingFactors.cpp` edit lines 13-16 to feed the program with your LP/IP/MIP.
- Running the file will give you Row and Column Scaling Factors.
 
### Notes 
- Optionally You can reach scaled A, B and C matrices, for that you need to uncomment certain parts of the code. While printing A, B, C matrices, printing to console is not recommended since it makes program slow down. So it is better to write them to the file.  
- You can use the scaling factors (R and C) or A, B and C for your own needs, for example, they can be fed into the program again and make the program solve already scaled program.
