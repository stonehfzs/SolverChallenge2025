#include <stdlib.h>
#include <string.h>
#include <math.h>

struct MySolver {
    int n;
    int* row_ptr;
    int* col_idx;
    double* val; // CSR格式
    
    // For direct solver - ILU factorization
    double* L_val;    // Lower triangular values
    double* U_val;    // Upper triangular values
    int* L_col_idx;   // Lower triangular column indices
    int* U_col_idx;   // Upper triangular column indices
    int* L_row_ptr;   // Lower triangular row pointers
    int* U_row_ptr;   // Upper triangular row pointers
    double* diag;     // Diagonal elements
    
    // For iterative solver - BiCGStab/GMRES
    double* r;        // residual vector
    double* p;        // search direction (BiCGStab)
    double* v;        // auxiliary vector (BiCGStab)
    double* s;        // auxiliary vector (BiCGStab)
    double* t;        // auxiliary vector (BiCGStab)
    double* z;        // preconditioned residual
    double* temp;     // temporary vector
    double* r_star;   // conjugate residual (BiCGStab)
    
    // For GMRES
    double** V;       // Krylov subspace basis (for GMRES)
    double* h;        // Hessenberg matrix (flattened)
    double* c;        // Givens rotation cosines
    double* s_givens; // Givens rotation sines
    double* g;        // RHS of least squares problem
    int restart;      // GMRES restart parameter
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
    double* v_re;    // auxiliary vector (BiCGStab)
    double* v_im;
    double* s_re;    // auxiliary vector (BiCGStab)
    double* s_im;
    double* t_re;    // auxiliary vector (BiCGStab)
    double* t_im;
    double* z_re;    // preconditioned residual
    double* z_im;
    double* temp_re; // temporary vector
    double* temp_im;
    double* r_star_re; // conjugate residual (BiCGStab)
    double* r_star_im;
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
    
    // Allocate memory for ILU factorization (simplified - use same structure as A)
    solver->L_val = (double*)malloc(nnz * sizeof(double));
    solver->U_val = (double*)malloc(nnz * sizeof(double));
    solver->L_col_idx = (int*)malloc(nnz * sizeof(int));
    solver->U_col_idx = (int*)malloc(nnz * sizeof(int));
    solver->L_row_ptr = (int*)malloc((n + 1) * sizeof(int));
    solver->U_row_ptr = (int*)malloc((n + 1) * sizeof(int));
    solver->diag = (double*)malloc(n * sizeof(double));
}

void direct_solver(struct MySolver* solver, const int n, const double* val)
{
    // Copy matrix values
    int nnz = solver->row_ptr[n];
    memcpy(solver->val, val, nnz * sizeof(double));
    
    // Simple ILU(0) with diagonal extraction for now
    // Extract diagonal elements for preconditioning
    for (int i = 0; i < n; i++) {
        solver->diag[i] = 1.0; // Default
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
    double* x_mutable = (double*)x;
    
    // Simple diagonal preconditioner (improved Jacobi) solution
    for (int i = 0; i < n; i++) {
        x_mutable[i] = b[i] / solver->diag[i];
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
        double* x_imag = (double*)x_im + i;
        
        // (a+bi) / (c+di) = (a+bi)*(c-di)/((c*c)+(d*d))
        double denom = solver->diag_re[i] * solver->diag_re[i] + 
                       solver->diag_im[i] * solver->diag_im[i];
        
        *x_re = (b[i] * solver->diag_re[i] + b_im[i] * solver->diag_im[i]) / denom;
        *x_imag = (b_im[i] * solver->diag_re[i] - b[i] * solver->diag_im[i]) / denom;
    }
}
#endif

#ifdef ITERATIVE_SOLVER
/* Note: Using spmv and spmv_complex functions from utlise.h instead of redefining them */

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
    
    // Allocate memory for BiCGStab method vectors
    solver->diag = (double*)malloc(n * sizeof(double));
    solver->r = (double*)malloc(n * sizeof(double));
    solver->p = (double*)malloc(n * sizeof(double));
    solver->v = (double*)malloc(n * sizeof(double));
    solver->s = (double*)malloc(n * sizeof(double));
    solver->t = (double*)malloc(n * sizeof(double));
    solver->z = (double*)malloc(n * sizeof(double));
    solver->temp = (double*)malloc(n * sizeof(double));
    solver->r_star = (double*)malloc(n * sizeof(double));
    
    // Initialize GMRES parameters
    solver->restart = 50; // GMRES restart parameter
    solver->V = (double**)malloc((solver->restart + 1) * sizeof(double*));
    for (int i = 0; i <= solver->restart; i++) {
        solver->V[i] = (double*)malloc(n * sizeof(double));
    }
    solver->h = (double*)malloc((solver->restart + 1) * solver->restart * sizeof(double));
    solver->c = (double*)malloc(solver->restart * sizeof(double));
    solver->s_givens = (double*)malloc(solver->restart * sizeof(double));
    solver->g = (double*)malloc((solver->restart + 1) * sizeof(double));
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
    double* r_star = solver->r_star;
    double* p = solver->p;
    double* v = solver->v;
    double* s = solver->s;
    double* t = solver->t;
    double* z = solver->z;
    double* temp = solver->temp;
    double* diag = solver->diag;
    
    const int max_iter = 2000;
    const double tol = 1e-10;
    
    // Initialize x to zero if not provided
    for (int i = 0; i < n; i++) {
        x_mutable[i] = 0.0;
    }
    
    // r = b - A*x (x starts as zero, so r = b)
    memcpy(r, b, n * sizeof(double));
    
    // Choose r_star = r
    memcpy(r_star, r, n * sizeof(double));
    
    // p = r
    memcpy(p, r, n * sizeof(double));
    
    double rho = 1.0, alpha = 1.0, omega = 1.0;
    
    double initial_residual = 0.0;
    for (int i = 0; i < n; i++) {
        initial_residual += r[i] * r[i];
    }
    initial_residual = sqrt(initial_residual);
    
    if (initial_residual < 1e-16) return; // Already converged
    
    // BiCGStab iteration
    for (int iter = 0; iter < max_iter; iter++) {
        double rho_old = rho;
        
        // rho = (r_star, r)
        rho = 0.0;
        for (int i = 0; i < n; i++) {
            rho += r_star[i] * r[i];
        }
        
        if (fabs(rho) < 1e-16) break; // Method fails
        
        double beta = (rho / rho_old) * (alpha / omega);
        
        // p = r + beta * (p - omega * v)
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        
        // Apply preconditioner: z = M^-1 * p (Jacobi)
        for (int i = 0; i < n; i++) {
            z[i] = p[i] / diag[i];
        }
        
        // v = A * z
        spmv(n, solver->row_ptr, solver->col_idx, solver->val, z, v);
        
        // alpha = rho / (r_star, v)
        double r_star_v = 0.0;
        for (int i = 0; i < n; i++) {
            r_star_v += r_star[i] * v[i];
        }
        
        if (fabs(r_star_v) < 1e-16) break; // Method fails
        
        alpha = rho / r_star_v;
        
        // s = r - alpha * v
        for (int i = 0; i < n; i++) {
            s[i] = r[i] - alpha * v[i];
        }
        
        // Check convergence on s
        double s_norm = 0.0;
        for (int i = 0; i < n; i++) {
            s_norm += s[i] * s[i];
        }
        s_norm = sqrt(s_norm);
        
        if (s_norm / initial_residual < tol) {
            // x = x + alpha * z
            for (int i = 0; i < n; i++) {
                x_mutable[i] += alpha * z[i];
            }
            break;
        }
        
        // Apply preconditioner: temp = M^-1 * s (Jacobi)
        for (int i = 0; i < n; i++) {
            temp[i] = s[i] / diag[i];
        }
        
        // t = A * temp
        spmv(n, solver->row_ptr, solver->col_idx, solver->val, temp, t);
        
        // omega = (t, s) / (t, t)
        double t_s = 0.0, t_t = 0.0;
        for (int i = 0; i < n; i++) {
            t_s += t[i] * s[i];
            t_t += t[i] * t[i];
        }
        
        if (fabs(t_t) < 1e-16) break; // Method fails
        
        omega = t_s / t_t;
        
        // x = x + alpha * z + omega * temp
        for (int i = 0; i < n; i++) {
            x_mutable[i] += alpha * z[i] + omega * temp[i];
        }
        
        // r = s - omega * t
        for (int i = 0; i < n; i++) {
            r[i] = s[i] - omega * t[i];
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
        
        if (fabs(omega) < 1e-16) break; // Method fails
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
    solver->v_re = (double*)malloc(n * sizeof(double));
    solver->v_im = (double*)malloc(n * sizeof(double));
    solver->s_re = (double*)malloc(n * sizeof(double));
    solver->s_im = (double*)malloc(n * sizeof(double));
    solver->t_re = (double*)malloc(n * sizeof(double));
    solver->t_im = (double*)malloc(n * sizeof(double));
    solver->z_re = (double*)malloc(n * sizeof(double));
    solver->z_im = (double*)malloc(n * sizeof(double));
    solver->temp_re = (double*)malloc(n * sizeof(double));
    solver->temp_im = (double*)malloc(n * sizeof(double));
    solver->r_star_re = (double*)malloc(n * sizeof(double));
    solver->r_star_im = (double*)malloc(n * sizeof(double));
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
    double* x_imag = (double*)x_im;
    double* r_re = solver->r_re;
    double* r_im = solver->r_im;
    double* r_star_re = solver->r_star_re;
    double* r_star_im = solver->r_star_im;
    double* p_re = solver->p_re;
    double* p_im = solver->p_im;
    double* v_re = solver->v_re;
    double* v_im = solver->v_im;
    double* s_re = solver->s_re;
    double* s_im = solver->s_im;
    double* t_re = solver->t_re;
    double* t_im = solver->t_im;
    double* z_re = solver->z_re;
    double* z_im = solver->z_im;
    double* temp_re = solver->temp_re;
    double* temp_im = solver->temp_im;
    
    const int max_iter = 2000;
    const double tol = 1e-10;
    
    // Initialize x to zero
    for (int i = 0; i < n; i++) {
        x_re[i] = 0.0;
        x_imag[i] = 0.0;
    }
    
    // r = b - A*x (x starts as zero, so r = b)
    memcpy(r_re, b, n * sizeof(double));
    memcpy(r_im, b_im, n * sizeof(double));
    
    // Choose r_star = r
    memcpy(r_star_re, r_re, n * sizeof(double));
    memcpy(r_star_im, r_im, n * sizeof(double));
    
    // p = r
    memcpy(p_re, r_re, n * sizeof(double));
    memcpy(p_im, r_im, n * sizeof(double));
    
    double rho_re = 1.0, rho_im = 0.0;
    double alpha_re = 1.0, alpha_im = 0.0;
    double omega_re = 1.0, omega_im = 0.0;
    
    double initial_residual = 0.0;
    for (int i = 0; i < n; i++) {
        initial_residual += r_re[i] * r_re[i] + r_im[i] * r_im[i];
    }
    initial_residual = sqrt(initial_residual);
    
    if (initial_residual < 1e-16) return; // Already converged
    
    // BiCGStab iteration for complex systems
    for (int iter = 0; iter < max_iter; iter++) {
        double rho_old_re = rho_re, rho_old_im = rho_im;
        
        // rho = (r_star, r) - complex dot product
        rho_re = 0.0; rho_im = 0.0;
        for (int i = 0; i < n; i++) {
            rho_re += r_star_re[i] * r_re[i] + r_star_im[i] * r_im[i];
            rho_im += r_star_re[i] * r_im[i] - r_star_im[i] * r_re[i];
        }
        
        double rho_norm = rho_re * rho_re + rho_im * rho_im;
        if (rho_norm < 1e-16) break; // Method fails
        
        // beta = (rho / rho_old) * (alpha / omega)
        double rho_ratio_re = (rho_re * rho_old_re + rho_im * rho_old_im) / 
                              (rho_old_re * rho_old_re + rho_old_im * rho_old_im);
        double rho_ratio_im = (rho_im * rho_old_re - rho_re * rho_old_im) / 
                              (rho_old_re * rho_old_re + rho_old_im * rho_old_im);
        
        double alpha_omega_re = (alpha_re * omega_re - alpha_im * omega_im);
        double alpha_omega_im = (alpha_re * omega_im + alpha_im * omega_re);
        double omega_norm = omega_re * omega_re + omega_im * omega_im;
        
        double beta_re = rho_ratio_re * alpha_omega_re / omega_norm - 
                         rho_ratio_im * alpha_omega_im / omega_norm;
        double beta_im = rho_ratio_re * alpha_omega_im / omega_norm + 
                         rho_ratio_im * alpha_omega_re / omega_norm;
        
        // p = r + beta * (p - omega * v)
        for (int i = 0; i < n; i++) {
            double temp1_re = p_re[i] - (omega_re * v_re[i] - omega_im * v_im[i]);
            double temp1_im = p_im[i] - (omega_re * v_im[i] + omega_im * v_re[i]);
            
            p_re[i] = r_re[i] + beta_re * temp1_re - beta_im * temp1_im;
            p_im[i] = r_im[i] + beta_re * temp1_im + beta_im * temp1_re;
        }
        
        // Apply preconditioner: z = M^-1 * p (complex Jacobi)
        for (int i = 0; i < n; i++) {
            double denom = solver->diag_re[i] * solver->diag_re[i] + 
                           solver->diag_im[i] * solver->diag_im[i];
            if (denom > 1e-16) {
                z_re[i] = (p_re[i] * solver->diag_re[i] + p_im[i] * solver->diag_im[i]) / denom;
                z_im[i] = (p_im[i] * solver->diag_re[i] - p_re[i] * solver->diag_im[i]) / denom;
            } else {
                z_re[i] = p_re[i];
                z_im[i] = p_im[i];
            }
        }
        
        // v = A * z
        spmv_complex(n, solver->row_ptr, solver->col_idx, 
                    solver->val_re, solver->val_im, 
                    z_re, z_im, v_re, v_im);
        
        // alpha = rho / (r_star, v)
        double r_star_v_re = 0.0, r_star_v_im = 0.0;
        for (int i = 0; i < n; i++) {
            r_star_v_re += r_star_re[i] * v_re[i] + r_star_im[i] * v_im[i];
            r_star_v_im += r_star_re[i] * v_im[i] - r_star_im[i] * v_re[i];
        }
        
        double r_star_v_norm = r_star_v_re * r_star_v_re + r_star_v_im * r_star_v_im;
        if (r_star_v_norm < 1e-16) break; // Method fails
        
        alpha_re = (rho_re * r_star_v_re + rho_im * r_star_v_im) / r_star_v_norm;
        alpha_im = (rho_im * r_star_v_re - rho_re * r_star_v_im) / r_star_v_norm;
        
        // s = r - alpha * v
        for (int i = 0; i < n; i++) {
            s_re[i] = r_re[i] - (alpha_re * v_re[i] - alpha_im * v_im[i]);
            s_im[i] = r_im[i] - (alpha_re * v_im[i] + alpha_im * v_re[i]);
        }
        
        // Check convergence on s
        double s_norm = 0.0;
        for (int i = 0; i < n; i++) {
            s_norm += s_re[i] * s_re[i] + s_im[i] * s_im[i];
        }
        s_norm = sqrt(s_norm);
        
        if (s_norm / initial_residual < tol) {
            // x = x + alpha * z
            for (int i = 0; i < n; i++) {
                x_re[i] += alpha_re * z_re[i] - alpha_im * z_im[i];
                x_imag[i] += alpha_re * z_im[i] + alpha_im * z_re[i];
            }
            break;
        }
        
        // Apply preconditioner: temp = M^-1 * s
        for (int i = 0; i < n; i++) {
            double denom = solver->diag_re[i] * solver->diag_re[i] + 
                           solver->diag_im[i] * solver->diag_im[i];
            if (denom > 1e-16) {
                temp_re[i] = (s_re[i] * solver->diag_re[i] + s_im[i] * solver->diag_im[i]) / denom;
                temp_im[i] = (s_im[i] * solver->diag_re[i] - s_re[i] * solver->diag_im[i]) / denom;
            } else {
                temp_re[i] = s_re[i];
                temp_im[i] = s_im[i];
            }
        }
        
        // t = A * temp
        spmv_complex(n, solver->row_ptr, solver->col_idx, 
                    solver->val_re, solver->val_im, 
                    temp_re, temp_im, t_re, t_im);
        
        // omega = (t, s) / (t, t)
        double t_s_re = 0.0, t_s_im = 0.0, t_t = 0.0;
        for (int i = 0; i < n; i++) {
            t_s_re += t_re[i] * s_re[i] + t_im[i] * s_im[i];
            t_s_im += t_re[i] * s_im[i] - t_im[i] * s_re[i];
            t_t += t_re[i] * t_re[i] + t_im[i] * t_im[i];
        }
        
        if (t_t < 1e-16) break; // Method fails
        
        omega_re = t_s_re / t_t;
        omega_im = t_s_im / t_t;
        
        // x = x + alpha * z + omega * temp
        for (int i = 0; i < n; i++) {
            x_re[i] += (alpha_re * z_re[i] - alpha_im * z_im[i]) + 
                       (omega_re * temp_re[i] - omega_im * temp_im[i]);
            x_imag[i] += (alpha_re * z_im[i] + alpha_im * z_re[i]) + 
                         (omega_re * temp_im[i] + omega_im * temp_re[i]);
        }
        
        // r = s - omega * t
        for (int i = 0; i < n; i++) {
            r_re[i] = s_re[i] - (omega_re * t_re[i] - omega_im * t_im[i]);
            r_im[i] = s_im[i] - (omega_re * t_im[i] + omega_im * t_re[i]);
        }
        
        // Check convergence
        double residual = 0.0;
        for (int i = 0; i < n; i++) {
            residual += r_re[i] * r_re[i] + r_im[i] * r_im[i];
        }
        residual = sqrt(residual);
        
        if (residual / initial_residual < tol) {
            break;
        }
        
        double omega_norm_check = omega_re * omega_re + omega_im * omega_im;
        if (omega_norm_check < 1e-16) break; // Method fails
    }
}
#endif
