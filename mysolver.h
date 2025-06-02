// please add your code in this file
struct MySolver {
    // add struct array
};

struct MySolver_complex {
    // add struct array
};

#ifdef DIRECT_SOLVER
//! real system
void preprocess(struct MySolver* solver, const int n, const int* row_ptr, const int* col_idx)
{
    // add the code
}

void direct_solver(struct MySolver* solver, const int n, const double* val)
{
    // add the code
}

void solve(struct MySolver* solver, const int n, const double* x, const double* b)
{
    // add the code
}

//! complex system
void preprocess_complex(struct MySolver_complex* solver, const int n, const int* row_ptr, const int* col_idx)
{
    // add the code
}

void direct_solver_complex(struct MySolver_complex* solver, const int n, const double* val, const double* val_im)
{
    // add the code
}

void solve_complex(struct MySolver_complex* solver, const int n, const double* x, const double* x_im, const double* b, const double* b_im)
{
    // add the code
}
#endif

#ifdef ITERATIVE_SOLVER
//! real system
void analyse(struct MySolver* solver, const int n, const int* row_ptr, const int* col_idx)
{
    // add the code
}

void preprocess(struct MySolver* solver, const int n, const double* val)
{
    // add the code
}

void iterative_solver(struct MySolver* solver, const int n, const double* x, const double* b)
{
    // add the code
}

//! complex system
void analyse_complex(struct MySolver_complex* solver, const int n, const int* row_ptr, const int* col_idx)
{
    // add the code
}

void preprocess_complex(struct MySolver_complex* solver, const int n, const double* val, const double* val_im)
{
    // add the code
}

void iterative_solver_complex(struct MySolver_complex* solver, const int n, const double* x, const double* x_im, const double* b, const double* b_im)
{
    // add the code
}
#endif
