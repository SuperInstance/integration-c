/*
 * Common spectral gap computation for integration tests.
 * Uses Jacobi eigenvalue algorithm for small symmetric matrices.
 */

/* Compute all eigenvalues of a symmetric matrix using Jacobi rotations.
   eigs: output array of length n (will be sorted ascending).
   A: input symmetric matrix, n×n (will be modified).
   Returns 0 on success. */
static int jacobi_eigenvalues(double *A, int n, double *eigs) {
    /* Copy to work array */
    double *W = malloc(n*n*sizeof(double));
    memcpy(W, A, n*n*sizeof(double));
    
    for (int sweep=0; sweep<100*n*n; sweep++) {
        /* Find largest off-diagonal element */
        double max_off = 0;
        int p=0, q=1;
        for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
            if (fabs(W[i*n+j])>max_off) { max_off=fabs(W[i*n+j]); p=i; q=j; }
        }
        if (max_off < 1e-14) break;
        
        /* Compute rotation angle */
        double app = W[p*n+p], aqq = W[q*n+q], apq = W[p*n+q];
        double theta;
        if (fabs(app-aqq) < 1e-30) {
            theta = M_PI/4;
        } else {
            theta = 0.5*atan2(2*apq, app-aqq);
        }
        double c = cos(theta), s = sin(theta);
        
        /* Apply rotation: W' = G^T W G */
        /* Only need to update rows/cols p and q */
        for (int i=0;i<n;i++) {
            if (i==p || i==q) continue;
            double wip = W[i*n+p], wiq = W[i*n+q];
            W[i*n+p] = c*wip + s*wiq;
            W[p*n+i] = W[i*n+p];
            W[i*n+q] = -s*wip + c*wiq;
            W[q*n+i] = W[i*n+q];
        }
        double new_pp = c*c*app + 2*s*c*apq + s*s*aqq;
        double new_qq = s*s*app - 2*s*c*apq + c*c*aqq;
        double new_pq = 0; /* By construction */
        W[p*n+p] = new_pp;
        W[q*n+q] = new_qq;
        W[p*n+q] = new_pq;
        W[q*n+p] = new_pq;
    }
    
    /* Extract diagonal as eigenvalues */
    for (int i=0;i<n;i++) eigs[i] = W[i*n+i];
    
    /* Sort ascending */
    for (int i=0;i<n-1;i++) for (int j=i+1;j<n;j++) {
        if (eigs[j] < eigs[i]) { double t=eigs[i]; eigs[i]=eigs[j]; eigs[j]=t; }
    }
    
    free(W);
    return 0;
}

/* Spectral gap = second smallest eigenvalue of Laplacian */
static double compute_spectral_gap(double *L, int n) {
    if (n<=1) return 0.0;
    double *eigs = malloc(n*sizeof(double));
    double *W = malloc(n*n*sizeof(double));
    memcpy(W, L, n*n*sizeof(double));
    /* Symmetrize (should already be symmetric) */
    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
        double avg = (W[i*n+j]+W[j*n+i])/2.0;
        W[i*n+j]=W[j*n+i]=avg;
    }
    jacobi_eigenvalues(W, n, eigs);
    double gap = eigs[1]; /* second smallest */
    free(eigs); free(W);
    return gap;
}
