#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include <cblas.h>
#include <lapacke.h>

// Macro for error handling
#define CHECK_STATUS(status) if (status) { fits_report_error(stderr, status); exit(status); }

typedef struct {
    double *data;     // Flattened data (n_frames x n_pixels), row-major
    double *time;     // Time stamps
    int n_frames;
    int n_pixels;     // width * height
    int width;
    int height;
} DataSet;

// Load Mask Data (2D FITS)
double* load_mask(const char *filename, int width, int height) {
    fitsfile *fptr;
    int status = 0;
    int naxis;
    long naxes[2];

    printf("Loading mask %s...\n", filename);

    if (fits_open_file(&fptr, filename, READONLY, &status)) {
        fits_report_error(stderr, status);
        exit(status);
    }

    if (fits_get_img_dim(fptr, &naxis, &status)) CHECK_STATUS(status);
    if (naxis != 2) {
        fprintf(stderr, "Error: Mask FITS file must have 2 dimensions (X, Y)\n");
        exit(1);
    }

    if (fits_get_img_size(fptr, 2, naxes, &status)) CHECK_STATUS(status);

    if (naxes[0] != width || naxes[1] != height) {
        fprintf(stderr, "Error: Mask dimensions (%ld x %ld) do not match dataset dimensions (%d x %d)\n",
                naxes[0], naxes[1], width, height);
        exit(1);
    }

    size_t data_size = (size_t)width * (size_t)height;
    double *mask = (double*) malloc(data_size * sizeof(double));
    if (!mask) {
        fprintf(stderr, "Error: Memory allocation failed for mask data\n");
        exit(1);
    }

    long fpixel[2] = {1, 1};
    if (fits_read_pix(fptr, TDOUBLE, fpixel, data_size, NULL, mask, NULL, &status)) CHECK_STATUS(status);

    fits_close_file(fptr, &status); CHECK_STATUS(status);

    return mask;
}

// Load FITS and timing data
DataSet* load_dataset(const char *fits_file, const char *txt_file) {
    fitsfile *fptr;
    int status = 0;
    int naxis;
    long naxes[3];

    printf("Loading %s and %s...\n", fits_file, txt_file);

    if (fits_open_file(&fptr, fits_file, READONLY, &status)) {
        fits_report_error(stderr, status);
        exit(status);
    }

    if (fits_get_img_dim(fptr, &naxis, &status)) CHECK_STATUS(status);
    if (naxis != 3) {
        fprintf(stderr, "Error: FITS file must have 3 dimensions (X, Y, Time)\n");
        exit(1);
    }

    if (fits_get_img_size(fptr, 3, naxes, &status)) CHECK_STATUS(status);

    DataSet *ds = (DataSet*) malloc(sizeof(DataSet));
    ds->width = naxes[0];
    ds->height = naxes[1];
    ds->n_frames = naxes[2];
    ds->n_pixels = ds->width * ds->height;

    printf("Dimensions: %d x %d x %d\n", ds->width, ds->height, ds->n_frames);

    // Allocate memory for data
    // Storing as (Time x Pixels)
    size_t data_size = (size_t)ds->n_frames * (size_t)ds->n_pixels;
    ds->data = (double*) malloc(data_size * sizeof(double));
    if (!ds->data) {
        fprintf(stderr, "Error: Memory allocation failed for image data\n");
        exit(1);
    }

    // Read data
    // FITS stores as Time blocks of (Y, X). We want to read all.
    // fits_read_img reads data into the array.
    // If the file is float or double, we can read directly.
    // We assume we want doubles.
    long fpixel[3] = {1, 1, 1};
    if (fits_read_pix(fptr, TDOUBLE, fpixel, data_size, NULL, ds->data, NULL, &status)) CHECK_STATUS(status);

    fits_close_file(fptr, &status); CHECK_STATUS(status);

    // Read Timing Data
    ds->time = (double*) malloc(ds->n_frames * sizeof(double));
    FILE *f = fopen(txt_file, "r");
    if (!f) {
        perror("Error opening timing file");
        exit(1);
    }

    for (int i = 0; i < ds->n_frames; i++) {
        if (fscanf(f, "%lf", &ds->time[i]) != 1) {
            fprintf(stderr, "Error reading timing file at line %d\n", i+1);
            exit(1);
        }
    }
    fclose(f);

    return ds;
}

void free_dataset(DataSet *ds) {
    if (ds) {
        if (ds->data) free(ds->data);
        if (ds->time) free(ds->time);
        free(ds);
    }
}

// Perform PCA to reduce dimensionality
// Input: ds (N x P)
// Output: reduced_data (N x k) allocated inside
// Algorithm:
// 1. Center the data (remove mean of each pixel over time? Or mean image?)
//    Usually for PCA on images, we subtract the mean image (average over frames).
// 2. Compute Gram matrix G = X * X^T (N x N)
//    (Since P >> N usually)
// 3. Eigen decomposition of G -> U, Lambda
// 4. Scores = U * sqrt(Lambda)
//    We want top k components.
void perform_pca(DataSet *ds, int k, double **reduced_data) {
    int N = ds->n_frames;
    int P = ds->n_pixels;

    printf("Performing PCA... N=%d, P=%d, k=%d\n", N, P, k);

    // 1. Compute Mean Image
    double *mean_image = (double*) calloc(P, sizeof(double));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < P; j++) {
            mean_image[j] += ds->data[i * P + j];
        }
    }
    for (int j = 0; j < P; j++) {
        mean_image[j] /= N;
    }

    // 2. Center Data (in place? No, let's create a temporary centered copy or do it on the fly?)
    // Creating a copy is memory intensive.
    // X * X^T where X_ij = Data_ij - Mean_j
    // (X X^T)_ab = sum_j (D_aj - M_j)(D_bj - M_j)
    //            = sum_j (D_aj D_bj - D_aj M_j - D_bj M_j + M_j^2)
    //            = (D D^T)_ab - N * (row_mean_a * row_mean_b)? No.
    // Let's just create centered data if memory allows.
    // If P is large, we should be careful.
    // But since we used `malloc` for `ds->data`, we have enough RAM for one copy?
    // Let's modify ds->data in place since we don't need original values?
    // Yes, for this task, centered data is fine.

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < P; j++) {
            ds->data[i * P + j] -= mean_image[j];
        }
    }
    free(mean_image);

    // 3. Compute Gram Matrix G = X * X^T
    // Use CBLAS. X is N x P (Row major).
    // G = X * X^T. Result is N x N.
    // cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)
    // C = A * B.
    // We want C = X * X^T.
    // A = X (N x P), B = X^T (P x N).
    // M=N, N=N, K=P.
    // Order = CblasRowMajor

    double *G = (double*) malloc(N * N * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                N, N, P,
                1.0, ds->data, P,
                ds->data, P,
                0.0, G, N);

    // 4. Eigen Decomposition
    // LAPACKE_dsyev(matrix_layout, jobz, uplo, n, a, lda, w)
    // jobz = 'V' (compute eigenvectors)
    // uplo = 'U' (upper triangle)
    // a = G (will be overwritten by eigenvectors)
    // w = eigenvalues (array of size N)

    double *w = (double*) malloc(N * sizeof(double));
    int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', N, G, N, w);
    if (info > 0) {
        fprintf(stderr, "The algorithm failed to compute eigenvalues.\n");
        exit(1);
    }

    // w is in ascending order. We want top k (last k).
    // Eigenvectors are in columns of G (if using FORTRAN/ColumnMajor) or Rows (if RowMajor)?
    // LAPACK_ROW_MAJOR: G contains eigenvectors as rows?
    // Documentation says: "If JOBZ = 'V', then if INFO = 0, A contains the orthonormal eigenvectors of the matrix A."
    // For Row Major, the eigenvectors are stored row-wise?
    // Wait, usually eigenvectors are columns.
    // In LAPACKE with Row Major, standard mapping applies.
    // If A is symmetric, U^T A U = Lambda.
    // The columns of U are eigenvectors.
    // But in Row Major storage, a column is not contiguous.
    // Let's assume standard LAPACK convention:
    // With LAPACK_ROW_MAJOR, A[i*lda + j] corresponds to matrix element (i,j).
    // The eigenvectors are returned in the same storage.
    // So the j-th eigenvector is the j-th column.

    // We want the last k eigenvectors corresponding to the largest k eigenvalues.
    // Eigenvalues w[N-1], w[N-2], ..., w[N-k].
    // Corresponding eigenvectors are in columns N-1, N-2, ...

    *reduced_data = (double*) malloc(N * k * sizeof(double));
    double *scores = *reduced_data;

    // Scores = U * S.
    // U are the eigenvectors of X X^T.
    // Wait, let X = U_x S_x V_x^T.
    // X X^T = U_x S_x^2 U_x^T.
    // So eigenvectors of G are U_x. Square root of eigenvalues are S_x.
    // The reduced representation (scores) we want is U_x S_x.
    // (Or just U_x if we want normalized temporal modes).
    // Let's use U_x * S_x.

    for (int i = 0; i < k; i++) {
        int eig_idx = N - 1 - i; // Index of eigenvalue (largest first)
        double lambda = w[eig_idx];
        double sigma = sqrt(lambda > 0 ? lambda : 0);

        // Copy eigenvector column `eig_idx` from G to column `i` of scores
        // G is N x N.
        for (int r = 0; r < N; r++) {
            // scores[r][i] = G[r][eig_idx] * sigma
            scores[r * k + i] = G[r * N + eig_idx] * sigma;
        }
    }

    free(G);
    free(w);
}

// Compute Canonical Correlation
// Input: DataA (N x k), DataB (N x k)
// Returns the first canonical correlation coefficient
double compute_cca_correlation(double *dataA, double *dataB, int n, int k) {
    if (n <= k) return 0.0; // Not enough samples

    // 1. Center the data (columns)
    // Create copies to work on
    double *A = (double*) malloc(n * k * sizeof(double));
    double *B = (double*) malloc(n * k * sizeof(double));
    memcpy(A, dataA, n * k * sizeof(double));
    memcpy(B, dataB, n * k * sizeof(double));

    for (int j = 0; j < k; j++) {
        double meanA = 0, meanB = 0;
        for (int i = 0; i < n; i++) {
            meanA += A[i * k + j];
            meanB += B[i * k + j];
        }
        meanA /= n;
        meanB /= n;
        for (int i = 0; i < n; i++) {
            A[i * k + j] -= meanA;
            B[i * k + j] -= meanB;
        }
    }

    // 2. QR Decomposition
    // We want Q_A and Q_B from QR decomposition of A and B.
    // A = Q_A R_A.
    // LAPACKE_dgeqrf computes QR.
    // LAPACKE_dorgqr generates Q.

    double *tauA = (double*) malloc(k * sizeof(double));
    double *tauB = (double*) malloc(k * sizeof(double));

    // dgeqrf(layout, m, n, a, lda, tau)
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, n, k, A, k, tauA);
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, n, k, B, k, tauB);

    // Generate Q. Q replaces A.
    // dorgqr(layout, m, n, k, a, lda, tau)
    // We only need the first k columns of Q (since A is n x k).
    // Actually, A is replaced by Q.
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, n, k, k, A, k, tauA);
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, n, k, k, B, k, tauB);

    // 3. Compute C = Q_A^T * Q_B
    // Q_A is n x k. Q_B is n x k.
    // C is k x k.
    double *C = (double*) malloc(k * k * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                k, k, n,
                1.0, A, k,
                B, k,
                0.0, C, k);

    // 4. SVD of C
    // The singular values of C are the canonical correlations.
    // We want the largest one.
    double *s = (double*) malloc(k * sizeof(double));
    double *superb = (double*) malloc((k - 1) * sizeof(double)); // For dgesvd if needed, but lapacke handles workspace?
    // Using dgesvd
    // dgesvd(layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb)
    // We don't need U and VT.

    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'N', 'N', k, k, C, k, s, NULL, 1, NULL, 1, superb);

    double max_corr = s[0]; // Singular values are sorted descending

    free(A); free(B);
    free(tauA); free(tauB);
    free(C); free(s); free(superb);

    return max_corr;
}

// Function to perform interpolation and correlation for a given offset
double evaluate_offset(DataSet *dsA, double *reducedA, DataSet *dsB, double *reducedB, int k, double offset) {
    // We shift B by offset: t_B_shifted = t_B + offset.
    // We want to compare A(t) with B(t + offset).
    // Or A(t) with B_shifted(t).
    // If t in A is t_A, we want B at t_A.
    // B's original times are t_B.
    // We interpolate B at t_A - offset?
    // Wait. If B is delayed by offset, B_shifted(t) = B(t - offset).
    // No, usually offset means time difference.
    // Let's say we want corr(A(t), B(t+offset)).
    // So for each t_A in A, we need value of B at (t_A + offset).

    int n_A = dsA->n_frames;
    int n_B = dsB->n_frames;
    double *t_A = dsA->time;
    double *t_B = dsB->time;

    // Find overlap range
    // We need t_A + offset to be within [t_B[0], t_B[n_B-1]]

    double t_B_min = t_B[0];
    double t_B_max = t_B[n_B-1];

    // Collect valid indices from A
    int *valid_indices = (int*) malloc(n_A * sizeof(int));
    int valid_count = 0;

    for (int i = 0; i < n_A; i++) {
        double target_t = t_A[i] + offset;
        if (target_t >= t_B_min && target_t <= t_B_max) {
            valid_indices[valid_count++] = i;
        }
    }

    if (valid_count < 3 * k || valid_count < 20) { // Minimum samples (prevent overfitting on small overlaps)
        free(valid_indices);
        return -1.0;
    }

    // Construct matrices for CCA
    double *subA = (double*) malloc(valid_count * k * sizeof(double));
    double *subB = (double*) malloc(valid_count * k * sizeof(double));

    for (int i = 0; i < valid_count; i++) {
        int idx_A = valid_indices[i];
        double target_t = t_A[idx_A] + offset;

        // Copy A row
        for (int j = 0; j < k; j++) {
            subA[i * k + j] = reducedA[idx_A * k + j];
        }

        // Interpolate B row
        // Binary search for target_t in t_B
        int low = 0, high = n_B - 2;
        int idx_B = 0;
        while (low <= high) {
            int mid = (low + high) / 2;
            if (t_B[mid] <= target_t && t_B[mid+1] >= target_t) {
                idx_B = mid;
                break;
            } else if (t_B[mid] < target_t) {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        // Fallback for edge cases (floating point precision)
        if (t_B[idx_B] > target_t || t_B[idx_B+1] < target_t) {
             // Linear scan if binary search failed slightly
             for(int ii=0; ii<n_B-1; ii++) {
                 if (t_B[ii] <= target_t && t_B[ii+1] >= target_t) {
                     idx_B = ii;
                     break;
                 }
             }
        }

        double t0 = t_B[idx_B];
        double t1 = t_B[idx_B+1];
        double alpha = (target_t - t0) / (t1 - t0);

        for (int j = 0; j < k; j++) {
            double v0 = reducedB[idx_B * k + j];
            double v1 = reducedB[(idx_B + 1) * k + j];
            subB[i * k + j] = v0 + alpha * (v1 - v0);
        }
    }

    double corr = compute_cca_correlation(subA, subB, valid_count, k);

    free(valid_indices);
    free(subA);
    free(subB);
    return corr;
}


int main(int argc, char *argv[]) {
    // Parameters
    const char *fileA_fits = "dataA.fits";
    const char *fileA_txt = "dataA.txt";
    const char *fileB_fits = "dataB.fits";
    const char *fileB_txt = "dataB.txt";
    const char *maskA_file = NULL;
    const char *maskB_file = NULL;

    int K = 10; // Number of PCA components

    // Parse command-line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--maskA") == 0) {
            if (i + 1 < argc) {
                maskA_file = argv[++i];
            } else {
                fprintf(stderr, "Error: --maskA requires a filename argument.\n");
                exit(1);
            }
        } else if (strcmp(argv[i], "--maskB") == 0) {
             if (i + 1 < argc) {
                maskB_file = argv[++i];
            } else {
                fprintf(stderr, "Error: --maskB requires a filename argument.\n");
                exit(1);
            }
        }
    }

    // Load Data A
    DataSet *dsA = load_dataset(fileA_fits, fileA_txt);

    // Apply Mask A if specified
    if (maskA_file) {
        double *maskA = load_mask(maskA_file, dsA->width, dsA->height);
        // Multiply
        for (int i = 0; i < dsA->n_frames; i++) {
            for (int j = 0; j < dsA->n_pixels; j++) {
                dsA->data[i * dsA->n_pixels + j] *= maskA[j];
            }
        }
        free(maskA);
        printf("Applied mask %s to Dataset A\n", maskA_file);
    }

    // Load Data B
    DataSet *dsB = load_dataset(fileB_fits, fileB_txt);

    // Apply Mask B if specified
    if (maskB_file) {
        double *maskB = load_mask(maskB_file, dsB->width, dsB->height);
        // Multiply
        for (int i = 0; i < dsB->n_frames; i++) {
            for (int j = 0; j < dsB->n_pixels; j++) {
                dsB->data[i * dsB->n_pixels + j] *= maskB[j];
            }
        }
        free(maskB);
        printf("Applied mask %s to Dataset B\n", maskB_file);
    }

    // Perform PCA
    double *reducedA, *reducedB;
    perform_pca(dsA, K, &reducedA);
    perform_pca(dsB, K, &reducedB);

    // Determine Scan Range
    // We scan offset such that A(t) ~ B(t + offset)
    // Offset range: min(t_B) - max(t_A) to max(t_B) - min(t_A)
    double min_offset = dsB->time[0] - dsA->time[dsA->n_frames-1];
    double max_offset = dsB->time[dsB->n_frames-1] - dsA->time[0];

    // Estimate step size (e.g., 1/10th of min frame duration)
    double dtA = (dsA->time[dsA->n_frames-1] - dsA->time[0]) / (dsA->n_frames - 1);
    double dtB = (dsB->time[dsB->n_frames-1] - dsB->time[0]) / (dsB->n_frames - 1);
    double step = (dtA < dtB ? dtA : dtB) / 5.0;

    printf("Scanning offsets from %.2f to %.2f with step %.4f\n", min_offset, max_offset, step);

    double best_offset = 0;
    double max_corr = -1.0;

    // To speed up, we can trim the search range to meaningful overlaps
    // But let's stick to the full range for now.

    for (double offset = min_offset; offset <= max_offset; offset += step) {
        double corr = evaluate_offset(dsA, reducedA, dsB, reducedB, K, offset);
        if (corr > max_corr) {
            max_corr = corr;
            best_offset = offset;
            // printf("Offset: %.4f, Corr: %.4f\n", offset, corr);
        }
    }

    // Refine search around best_offset
    double refined_step = step / 20.0;
    double search_width = step * 2.0;

    printf("Refining search around %.4f...\n", best_offset);

    for (double offset = best_offset - search_width; offset <= best_offset + search_width; offset += refined_step) {
         double corr = evaluate_offset(dsA, reducedA, dsB, reducedB, K, offset);
         if (corr > max_corr) {
            max_corr = corr;
            best_offset = offset;
         }
    }

    printf("\nBest Time Offset: %.6f\n", best_offset);
    printf("Max Correlation: %.6f\n", max_corr);

    // Cleanup
    free(reducedA);
    free(reducedB);
    free_dataset(dsA);
    free_dataset(dsB);

    return 0;
}
