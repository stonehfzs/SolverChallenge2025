#include <stdlib.h>
#include <string.h>
#include <math.h>

struct MySolver {
    int n;
    int* row_ptr;
    int* col_idx;
    double* val; // CSR格式
    // LU分解结果（这里只存储对角元，适用于对角占优矩阵的演示）
    double* diag;
    
    // Additional fields for iterative solver
    double* r;     // residual vector
    double* p;     // search direction
    double* z;     // preconditioned residual
    double* temp;  // temporary vector
};

struct MySolver_complex {
    int n;
    int* row_ptr;
    int* col_idx;
    double* val_re;
    double* val_im;
    double* diag_re;
    double* diag_im;
    
    // Additional fields for iterative solver
    double* r_re;    // residual vector
    double* r_im;
    double* p_re;    // search direction
    double* p_im;
    double* z_re;    // preconditioned residual
    double* z_im;
    double* temp_re; // temporary vector
    double* temp_im;
};

#ifdef DIRECT_SOLVER
//! real system
void preprocess(struct MySolver* solver, const int n, const int* row_ptr, const int* col_idx)
{
    solver->n = n;
    
    // Allocate and copy CSR structure
    solver->row_ptr = (int*)malloc((n + 1) * sizeof(int));
    memcpy(solver->row_ptr, row_ptr, (n + 1) * sizeof(int));
    
    int nnz = row_ptr[n];
    solver->col_idx = (int*)malloc(nnz * sizeof(int));
    memcpy(solver->col_idx, col_idx, nnz * sizeof(int));
    
    // Allocate memory for values
    solver->val = (double*)malloc(nnz * sizeof(double));
    
    // Allocate diagonal array for simple preconditioner
    solver->diag = (double*)malloc(n * sizeof(double));
}

void direct_solver(struct MySolver* solver, const int n, const double* val)
{
    // Copy matrix values
    int nnz = solver->row_ptr[n];
    memcpy(solver->val, val, nnz * sizeof(double));
    
    // Extract diagonal elements for simple preconditioning
    for (int i = 0; i < n; i++) {
        solver->diag[i] = 0.0;
        for (int j = solver->row_ptr[i]; j < solver->row_ptr[i + 1]; j++) {
            if (solver->col_idx[j] == i) {
                solver->diag[i] = solver->val[j];
                break;
            }
        }
        // Ensure diagonal is not zero (fallback if zero)
        if (fabs(solver->diag[i]) < 1e-14) {
            solver->diag[i] = 1.0;
        }
    }
}

void solve(struct MySolver* solver, const int n, const double* x, const double* b)
{
    // Simple diagonal preconditioner (Jacobi) solution
    for (int i = 0; i < n; i++) {
        double* x_i = (double*)x + i;
        *x_i = b[i] / solver->diag[i];
    }
}

//! complex system
void preprocess_complex(struct MySolver_complex* solver, const int n, const int* row_ptr, const int* col_idx)
{
    solver->n = n;
    
    // Allocate and copy CSR structure
    solver->row_ptr = (int*)malloc((n + 1) * sizeof(int));
    memcpy(solver->row_ptr, row_ptr, (n + 1) * sizeof(int));
    
    int nnz = row_ptr[n];
    solver->col_idx = (int*)malloc(nnz * sizeof(int));
    memcpy(solver->col_idx, col_idx, nnz * sizeof(int));
    
    // Allocate memory for values
    solver->val_re = (double*)malloc(nnz * sizeof(double));
    solver->val_im = (double*)malloc(nnz * sizeof(double));
    
    // Allocate diagonal array for simple preconditioner
    solver->diag_re = (double*)malloc(n * sizeof(double));
    solver->diag_im = (double*)malloc(n * sizeof(double));
}

void direct_solver_complex(struct MySolver_complex* solver, const int n, const double* val, const double* val_im)
{
    // Copy matrix values
    int nnz = solver->row_ptr[n];
    memcpy(solver->val_re, val, nnz * sizeof(double));
    memcpy(solver->val_im, val_im, nnz * sizeof(double));
    
    // Extract diagonal elements for simple preconditioning
    for (int i = 0; i < n; i++) {
        solver->diag_re[i] = 1.0;
        solver->diag_im[i] = 0.0;
        
        for (int j = solver->row_ptr[i]; j < solver->row_ptr[i + 1]; j++) {
            if (solver->col_idx[j] == i) {
                solver->diag_re[i] = solver->val_re[j];
                solver->diag_im[i] = solver->val_im[j];
                break;
            }
        }
        
        // Ensure diagonal is not zero (fallback if zero)
        if (fabs(solver->diag_re[i]) < 1e-14 && fabs(solver->diag_im[i]) < 1e-14) {
            solver->diag_re[i] = 1.0;
            solver->diag_im[i] = 0.0;
        }
    }
}

void solve_complex(struct MySolver_complex* solver, const int n, const double* x, const double* x_im, const double* b, const double* b_im)
{
    // Simple diagonal preconditioner solution for complex systems
    for (int i = 0; i < n; i++) {
        double* x_re = (double*)x + i;
        double* x_im = (double*)x_im + i;
        
        // (a+bi) / (c+di) = (a+bi)*(c-di)/((c*c)+(d*d))
        double denom = solver->diag_re[i] * solver->diag_re[i] + 
                       solver->diag_im[i] * solver->diag_im[i];
        
        *x_re = (b[i] * solver->diag_re[i] + b_im[i] * solver->diag_im[i]) / denom;
        *x_im = (b_im[i] * solver->diag_re[i] - b[i] * solver->diag_im[i]) / denom;
    }
}
#endif

#ifdef ITERATIVE_SOLVER
// Matrix-vector multiplication (CSR format) for real systems
void spmv(const int n, const int* row_ptr, const int* col_idx, const double* val, 
          const double* x, double* y) {
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            y[i] += val[j] * x[col_idx[j]];
        }
    }
}

// Matrix-vector multiplication (CSR format) for complex systems
void spmv_complex(const int n, const int* row_ptr, const int* col_idx, 
                 const double* val_re, const double* val_im,
                 const double* x_re, const double* x_im,
                 double* y_re, double* y_im) {
    for (int i = 0; i < n; i++) {
        y_re[i] = 0.0;
        y_im[i] = 0.0;
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            int k = col_idx[j];
            // (a+bi)(c+di) = (ac-bd) + (ad+bc)i
            y_re[i] += val_re[j] * x_re[k] - val_im[j] * x_im[k];
            y_im[i] += val_re[j] * x_im[k] + val_im[j] * x_re[k];
        }
    }
}

//! real system
void analyse(struct MySolver* solver, const int n, const int* row_ptr, const int* col_idx)
{
    solver->n = n;
    
    // Allocate and copy CSR structure
    solver->row_ptr = (int*)malloc((n + 1) * sizeof(int));
    memcpy(solver->row_ptr, row_ptr, (n + 1) * sizeof(int));
    
    int nnz = row_ptr[n];
    solver->col_idx = (int*)malloc(nnz * sizeof(int));
    memcpy(solver->col_idx, col_idx, nnz * sizeof(int));
    
    // Allocate memory for values
    solver->val = (double*)malloc(nnz * sizeof(double));
    
    // Allocate memory for CG method vectors
    solver->diag = (double*)malloc(n * sizeof(double));
    solver->r = (double*)malloc(n * sizeof(double));
    solver->p = (double*)malloc(n * sizeof(double));
    solver->z = (double*)malloc(n * sizeof(double));
    solver->temp = (double*)malloc(n * sizeof(double));
}

void preprocess(struct MySolver* solver, const int n, const double* val)
{
    // Copy matrix values
    int nnz = solver->row_ptr[n];
    memcpy(solver->val, val, nnz * sizeof(double));
    
    // Extract diagonal elements for Jacobi preconditioner
    for (int i = 0; i < n; i++) {
        solver->diag[i] = 1.0; // Default
        for (int j = solver->row_ptr[i]; j < solver->row_ptr[i + 1]; j++) {
            if (solver->col_idx[j] == i) {
                solver->diag[i] = solver->val[j];
                break;
            }
        }
        // Ensure diagonal is not zero
        if (fabs(solver->diag[i]) < 1e-14) {
            solver->diag[i] = 1.0;
        }
    }
}

void iterative_solver(struct MySolver* solver, const int n, const double* x, const double* b)
{
    double* x_mutable = (double*)x;
    double* r = solver->r;
    double* p = solver->p;
    double* z = solver->z;
    double* temp = solver->temp;
    double* diag = solver->diag;
    
    const int max_iter = 1000;
    const double tol = 1e-8;
    
    // Initialize x to zero if not provided
    for (int i = 0; i < n; i++) {
        x_mutable[i] = 0.0;
    }
    
    // r = b - A*x
    spmv(n, solver->row_ptr, solver->col_idx, solver->val, x_mutable, temp);
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - temp[i];
    }
    
    // Apply preconditioner z = M^-1 * r (Jacobi)
    for (int i = 0; i < n; i++) {
        z[i] = r[i] / diag[i];
    }
    
    // p = z
    memcpy(p, z, n * sizeof(double));
    
    double rz = 0.0;
    for (int i = 0; i < n; i++) {
        rz += r[i] * z[i];
    }
    
    double initial_residual = 0.0;
    for (int i = 0; i < n; i++) {
        initial_residual += r[i] * r[i];
    }
    initial_residual = sqrt(initial_residual);
    
    // CG iteration
    for (int iter = 0; iter < max_iter; iter++) {
        // temp = A * p
        spmv(n, solver->row_ptr, solver->col_idx, solver->val, p, temp);
        
        // alpha = (r,z)/(p,temp)
        double ptemp = 0.0;
        for (int i = 0; i < n; i++) {
            ptemp += p[i] * temp[i];
        }
        double alpha = rz / ptemp;
        
        // x = x + alpha*p
        for (int i = 0; i < n; i++) {
            x_mutable[i] += alpha * p[i];
        }
        
        // r = r - alpha*temp
        for (int i = 0; i < n; i++) {
            r[i] -= alpha * temp[i];
        }
        
        // Check convergence
        double residual = 0.0;
        for (int i = 0; i < n; i++) {
            residual += r[i] * r[i];
        }
        residual = sqrt(residual);
        
        if (residual / initial_residual < tol) {
            break;
        }
        
        // Apply preconditioner z = M^-1 * r
        for (int i = 0; i < n; i++) {
            z[i] = r[i] / diag[i];
        }
        
        // beta = (r_new,z_new)/(r_old,z_old)
        double rz_new = 0.0;
        for (int i = 0; i < n; i++) {
            rz_new += r[i] * z[i];
        }
        double beta = rz_new / rz;
        rz = rz_new;
        
        // p = z + beta*p
        for (int i = 0; i < n; i++) {
            p[i] = z[i] + beta * p[i];
        }
    }
}

//! complex system
void analyse_complex(struct MySolver_complex* solver, const int n, const int* row_ptr, const int* col_idx)
{
    solver->n = n;
    
    // Allocate and copy CSR structure
    solver->row_ptr = (int*)malloc((n + 1) * sizeof(int));
    memcpy(solver->row_ptr, row_ptr, (n + 1) * sizeof(int));
    
    int nnz = row_ptr[n];
    solver->col_idx = (int*)malloc(nnz * sizeof(int));
    memcpy(solver->col_idx, col_idx, nnz * sizeof(int));
    
    // Allocate memory for values
    solver->val_re = (double*)malloc(nnz * sizeof(double));
    solver->val_im = (double*)malloc(nnz * sizeof(double));
    
    // Allocate diagonal array for Jacobi preconditioner
    solver->diag_re = (double*)malloc(n * sizeof(double));
    solver->diag_im = (double*)malloc(n * sizeof(double));
    
    // Allocate memory for BiCGStab method vectors
    solver->r_re = (double*)malloc(n * sizeof(double));
    solver->r_im = (double*)malloc(n * sizeof(double));
    solver->p_re = (double*)malloc(n * sizeof(double));
    solver->p_im = (double*)malloc(n * sizeof(double));
    solver->z_re = (double*)malloc(n * sizeof(double));
    solver->z_im = (double*)malloc(n * sizeof(double));
    solver->temp_re = (double*)malloc(n * sizeof(double));
    solver->temp_im = (double*)malloc(n * sizeof(double));
}

void preprocess_complex(struct MySolver_complex* solver, const int n, const double* val, const double* val_im)
{
    // Copy matrix values
    int nnz = solver->row_ptr[n];
    memcpy(solver->val_re, val, nnz * sizeof(double));
    memcpy(solver->val_im, val_im, nnz * sizeof(double));
    
    // Extract diagonal elements for Jacobi preconditioner
    for (int i = 0; i < n; i++) {
        solver->diag_re[i] = 1.0;
        solver->diag_im[i] = 0.0;
        
        for (int j = solver->row_ptr[i]; j < solver->row_ptr[i + 1]; j++) {
            if (solver->col_idx[j] == i) {
                solver->diag_re[i] = solver->val_re[j];
                solver->diag_im[i] = solver->val_im[j];
                break;
            }
        }
        
        // Ensure diagonal is not zero
        if (fabs(solver->diag_re[i]) < 1e-14 && fabs(solver->diag_im[i]) < 1e-14) {
            solver->diag_re[i] = 1.0;
            solver->diag_im[i] = 0.0;
        }
    }
}

void iterative_solver_complex(struct MySolver_complex* solver, const int n, 
                              const double* x, const double* x_im, 
                              const double* b, const double* b_im)
{
    double* x_re = (double*)x;
    double* x_im = (double*)x_im;
    double* r_re = solver->r_re;
    double* r_im = solver->r_im;
    double* p_re = solver->p_re;
    double* p_im = solver->p_im;
    
    // Initialize to simple diagonal preconditioner solution
    for (int i = 0; i < n; i++) {
        double denom = solver->diag_re[i] * solver->diag_re[i] + 
                       solver->diag_im[i] * solver->diag_im[i];
        
        // (a+bi) / (c+di) = (a+bi)*(c-di)/((c*c)+(d*d))
        x_re[i] = (b[i] * solver->diag_re[i] + b_im[i] * solver->diag_im[i]) / denom;
        x_im[i] = (b_im[i] * solver->diag_re[i] - b[i] * solver->diag_im[i]) / denom;
    }
    
    // Note: A full BiCGStab implementation for complex systems would be more involved
    // This is a simplified placeholder implementation
    
    // Calculate residual r = b - A*x
    spmv_complex(n, solver->row_ptr, solver->col_idx, 
                solver->val_re, solver->val_im, 
                x_re, x_im, 
                solver->temp_re, solver->temp_im);
    
    for (int i = 0; i < n; i++) {
        r_re[i] = b[i] - solver->temp_re[i];
        r_im[i] = b_im[i] - solver->temp_im[i];
    }
}
#endif
