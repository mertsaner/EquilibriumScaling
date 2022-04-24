#pragma once
#include "External.h"
#include "ClassIncludes.h"
#include </home/mert/snap/eigen-3.4.0/Eigen/Eigen>
#include </home/mert/snap/eigen-3.4.0/Eigen/Sparse>

using namespace std;

using TripletVector = vector<Eigen::Triplet<double>>;
using SparseMatrix  = Eigen::SparseMatrix<double>;
using Vector        = vector<double>;

tuple<SparseMatrix, Vector, Vector> read_Abc2(GRBModel& model) {
    // number of variables, number of constraints, number of nonzeros in A
    size_t n = model.get(GRB_IntAttr_NumVars); //number of variables
    size_t m = model.get(GRB_IntAttr_NumConstrs); //number of constraints
    size_t num_nnz = model.get(GRB_IntAttr_NumNZs); //number of non zeros

    // allocate space
    TripletVector triplet(num_nnz);
    Vector b(m);
    Vector c(n);

    // Read the objective coefficients
    auto obj_expr = model.getObjective().getLinExpr();

    for ( size_t j = 0; j < n; ++j) {
        c[j] = obj_expr.getCoeff(j);
    }

    // Read the coefficient matrix A and RHS vector b
    size_t k = 0;
    for (size_t i = 0; i < m; ++i) {
        auto con = model.getConstr(i);
        auto row = model.getRow(con);

        for (size_t r = 0; r < row.size(); ++r) {
            auto col_idx = row.getVar(r).index();
            auto coeff   = row.getCoeff(r);
            triplet[k++] = Eigen::Triplet<double>(i, col_idx, coeff);
        }
        // read the constraint RHS
        b[i] = con.get(GRB_DoubleAttr_RHS);
    }
    // Construct the sparse matrix.
    SparseMatrix A(m, n);
    A.setFromTriplets(triplet.begin(), triplet.end());
    return {A, b, c};
}

void afterScaling () {
        ptree root;
        string p = (string) getenv("PROJECT_PATH") + "/data/ShortHeat2V1S1.xml";
        initialize_ptree(root,  p);
    // initialize DetailedModelGRB
    DetailedModelGRB d(root);
    GRBModel model = d.get_model(); // original model
    GRBEnv env = model.getEnv();
    int n = model.get(GRB_IntAttr_NumVars); // number of variables (columns)
    int m = model.get(GRB_IntAttr_NumConstrs); // number of constraints (rows)

    //First Apply Row Scaling
    auto[A, b, c] = read_Abc2(model);
    double row_max[m];     //Declaring One dimensional scaling matrices for rows like int arr[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double row_multi[m];
    double col_max[n];     //Declaring One dimensional scaling matrices for columns
    double col_multi[n];
    double max = 0; //Maximum value in the row
    int cnt = 0;    //Counter
    bool exitLoop = false;  //Boolean for zero row control

    //Apply Row Scaling
    for (int i = 0; i < m; i++) {
        row_max[i] = 0;   //Creating [m][1] zero matrices
        row_multi[i] = 0; //Creating [m][1] zeros matrices
        for (int j = 0; j < n; j++) {
            if (A.coeff(i, j) == 0.0) {         //Finding if a row contains at least one nonzero element
                cnt = cnt + 1;
                if (cnt == n) //Since we are only checking rows if cnt==n then there is a row that is full of zeros
                {
                    exitLoop = true;
                    break; //Break the for loop if all rows are zero
                }
            }
            if (!exitLoop) { //Scaling Procedure Started for Matrices whose rows have at least one non-zero element
                if (abs(A.coeff(i, j)) > max) { //Absolute value of the element in the row
                    max = abs(A.coeff(i, j));
                    row_max[i] = max; //Adding the biggest abs value to row_max (1D matrix)
                }
            }
        }

        //if (col_max[j] != 0 && (col_max[j] > 0.000000000000001  || col_max[j] < (-0.000000000000001)))
        if (row_max[i] != 0) {
            row_multi[i] = 1.0 / row_max[i]; //1.0 Written for fractional Division
        }

        if (row_multi[i] == 0)
        {
            row_multi[i] = 1;
        }

        // b[i] = b[i] * row_multi[i];

        // Multiply each element of A with 1/the biggest value in that row (row multi)
        for (int j = 0; j < n; j++) {
            A.coeffRef(i, j) = A.coeffRef(i, j) * row_multi[i];
        }

        cnt = 0; //Row counter reset
        max = 0; //Row max reset
        exitLoop = false; //Boolean for zero row control reset
    }


    vector<double> R;
    for (double ia: row_multi) {
        R.push_back(ia);
    }

    //Apply Column Scaling
    cnt = 0; //Row counter reset
    max = 0; //Row max reset
    exitLoop = false;  //Boolean for zero column/row control reset

    for (int j = 0; j < n; j++) {
        col_max[j] = 0;   //Creating [1][n] zero matrices
        col_multi[j] = 0; //Creating [1][n] zeros matrices
        for (int i = 0; i < m; i++) {
            if (A.coeff(i, j) == 0.0) {         //Finding if a column contains at least one nonzero element
                cnt = cnt + 1;
                if (cnt == m) //Since we are only checking columns if cnt==m then there is a column that is full of zeros
                {
                    // cout<<"Entered 2" << endl;
                    exitLoop = true;
                    break; //Break the for loop if all columns are zero
                }
            }

            if (!exitLoop) { //Scaling Procedure Started for Matrices whose rows have at least one non-zero element
                if (abs(A.coeff(i, j)) > max) { //Absolute value of the element in the row
                    max = abs(A.coeff(i, j));
                    col_max[j] = max; //Adding the biggest abs value to row_max (1D matrix)
                }
            }
        }

        //Calculation of specific row scaling factor
        if (col_max[j] != 0) {
            col_multi[j] = 1.0 / col_max[j]; //1.0 Written for fractional Division
        }

        if (col_multi[j] == 0)
        {
            col_multi[j] = 1;
        }

        //   c[j] = c[j] * col_multi[j];

        // Multiply each element of A with 1/the biggest value in that row (row multi)
        //   for (int i = 0; i < m; i++) {
        //       A.coeffRef(i, j) = A.coeffRef(i, j) * col_multi[j];
        //    }

        cnt = 0; //Row counter reset
        max = 0; //Row max reset
        exitLoop = false; //Boolean for zero column control reset
    }

    // create an empty vector and push all elements
    vector<double> C;
    for (double jb: col_multi) {
        C.push_back(jb);
    }

//End of Scaling
//================================================
    assert(R.size() == model.get(GRB_IntAttr_NumConstrs));
    for (double ra : R){
        assert(ra > 0);
    }

    assert(C.size()==model.get(GRB_IntAttr_NumVars));
    for (double ca : C){
        assert(ca >= 0);
    }

    GRBModel model_scaled(env); // scaled model which is filled in the following

// create scaled variables
    for (int j = 0; j < model.get(GRB_IntAttr_NumVars); j++){
        GRBVar v = model.getVar(j); // original variable
        double UB = v.get(GRB_DoubleAttr_UB); // original upper bound
        double LB = v.get(GRB_DoubleAttr_LB); // original lower bound
        double Obj = v.get(GRB_DoubleAttr_Obj); // original obj coeff
        string v_name = v.get(GRB_StringAttr_VarName); // original variable name
        char v_type = v.get(GRB_CharAttr_VType);

        if (v_type != 'C') C.at(j) = 1; // integer columns are not scaled
        model_scaled.addVar(1/C.at(j)*LB, 1/C.at(j)*UB, Obj*C.at(j), v_type, v_name + "_scaled");
    }
    model_scaled.update();

// create scaled constraints
    for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); i++){ // i: row index
        GRBConstr constr = model.getConstr(i);
        string constr_name = constr.get(GRB_StringAttr_ConstrName);

        GRBLinExpr lin_expr = model.getRow(constr);
        GRBLinExpr lin_expr_scaled = 0; // new linear expression

        for (int k = 0; k < lin_expr.size(); k++){ // k: local index in linear expression
            int j = lin_expr.getVar(k).index(); // column index
            GRBVar v_scaled = model_scaled.getVar(j); // scaled variable

            double coeff = lin_expr.getCoeff(k); // = Coefficient_Matrix(i,j)
            double scaled_coeff = R.at(i) * coeff * C.at(j); // calculate scaled matrix coefficient

            lin_expr_scaled += scaled_coeff * v_scaled; // add term to linear expression
        }
        double rhs = constr.get(GRB_DoubleAttr_RHS); // original right-hand side
        double rhs_scaled = R.at(i) * rhs; // scaled right-hand side

        // model constraints
        char sense = constr.get(GRB_CharAttr_Sense); // get constraint sense
        if (sense == '>'){
            model_scaled.addConstr(lin_expr_scaled >= rhs_scaled, constr_name + "_scaled");
        }
        if (sense == '='){
            model_scaled.addConstr(lin_expr_scaled == rhs_scaled, constr_name + "_scaled");
        }
        if (sense == '<'){
            model_scaled.addConstr(lin_expr_scaled <= rhs_scaled, constr_name + "_scaled");
        }
    }
    model_scaled.update();

// both should produce the same results
    model.optimize();
    model_scaled.optimize();

    cout << model.get(GRB_DoubleAttr_ObjVal) << endl;
    cout << model_scaled.get(GRB_DoubleAttr_ObjVal) << endl;
}
