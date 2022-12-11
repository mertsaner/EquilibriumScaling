## Equilibrium Scaling

### Prerequisites
- To run the program firstly Gurobi C++ API should be ready. To install Gurobi Solver [click here.](https://www.gurobi.com/free-trial/ "Gurobi Solver Website")
- The program requires the Eigen library to read large .lp files. Hence you need to include the Eigen library in the program. For Eigen library download instructions [click here.](https://eigen.tuxfamily.org/ "Eigen library download")

### Usage
- To configure the program according to your needs, you need to edit the `ScalingFactors.cpp` file.
- In `ScalingFactors.cpp` edit lines 13-16 to feed the program with your LP/IP/MIP.
- Running the file will give you Row and Column Scaling Factors.
- `Small_Investment_Model_presolved.lp`, `my.lp` and `trial.lp` are three sample .lp files to test the program.
 
### Notes 
- `Scaling.h`is a header file that compares the result of the scaled program with the original program. To do that it uses a .xml file instead of .lp files like ScalingFactors.cpp.
- Optionally You can reach scaled A, B and C matrices, for that you need to uncomment certain parts of the code. While printing A, B, and C matrices, printing to the console is not recommended since it makes the program slow down. So it is better to write them to the file.  
- You can reuse the scaling factors (R and C) or A, B and C depending on your own needs. For example, they can be fed into the program again and make the program solve the already scaled program. You can edit/update `Scaling.h` to compare result of scaled and unscaled versions.
