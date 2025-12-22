#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <ctype.h>
#include <math.h>
#include <fitsio.h>
#include <unistd.h>

// Macro for error handling
#define CHECK_STATUS(status) if (status) { fits_report_error(stderr, status); exit(status); }

#define MAX_LINE_LENGTH 1024
#define MAX_FILENAME 512

typedef struct {
    char filename[MAX_FILENAME];
} FileInfo;

typedef struct {
    long global_idx;
    double timestamp;
    char fits_filename[MAX_FILENAME];
    long local_cube_idx;
} FrameInfo;

typedef struct {
    double *times;
    double *times_end;
    double *diffs; // SumSq diffs between pairs
    int n_pairs;   // Number of pairs
    int n_frames;  // Number of original frames
} DiffResult;

int compare_filenames(const void *a, const void *b) {
    const FileInfo *fa = (const FileInfo *)a;
    const FileInfo *fb = (const FileInfo *)b;
    return strcmp(fa->filename, fb->filename);
}

void get_fits_filename(const char *txt_filename, char *fits_filename) {
    strcpy(fits_filename, txt_filename);
    char *dot = strrchr(fits_filename, '.');
    if (dot) {
        strcpy(dot, ".fits");
    } else {
        strcat(fits_filename, ".fits");
    }
}

double* load_mask(const char *stream_name, int *width, int *height) {
    char mask_filename[MAX_FILENAME];
    snprintf(mask_filename, sizeof(mask_filename), "%s.mask.fits", stream_name);

    if (access(mask_filename, F_OK) != 0) {
        return NULL;
    }

    fitsfile *fptr;
    int status = 0;
    int naxis;
    long naxes[2];

    printf("Loading mask %s...\n", mask_filename);
    if (fits_open_file(&fptr, mask_filename, READONLY, &status)) {
        fits_report_error(stderr, status);
        return NULL;
    }

    if (fits_get_img_dim(fptr, &naxis, &status)) CHECK_STATUS(status);
    if (naxis != 2) {
        fprintf(stderr, "Error: Mask FITS file must have 2 dimensions (X, Y)\n");
        exit(1);
    }
    if (fits_get_img_size(fptr, 2, naxes, &status)) CHECK_STATUS(status);

    *width = naxes[0];
    *height = naxes[1];

    size_t data_size = (*width) * (*height);
    double *mask = (double*) malloc(data_size * sizeof(double));
    long fpixel[2] = {1, 1};
    if (fits_read_pix(fptr, TDOUBLE, fpixel, data_size, NULL, mask, NULL, &status)) CHECK_STATUS(status);

    fits_close_file(fptr, &status); CHECK_STATUS(status);
    return mask;
}

void read_frame_data(const char *stream_name, FrameInfo *frame, int expected_pixels, double *buffer, int *width_chk, int *height_chk) {
    char filepath[MAX_FILENAME * 2];
    snprintf(filepath, sizeof(filepath), "%s/%s", stream_name, frame->fits_filename);

    fitsfile *fptr;
    int status = 0;

    if (fits_open_file(&fptr, filepath, READONLY, &status)) {
        fits_report_error(stderr, status);
        memset(buffer, 0, expected_pixels * sizeof(double));
        return;
    }

    if (width_chk && height_chk) {
        int naxis;
        long naxes[3];
        if (fits_get_img_dim(fptr, &naxis, &status) == 0 && fits_get_img_size(fptr, 3, naxes, &status) == 0) {
             if (naxes[0] != *width_chk || naxes[1] != *height_chk) {
                 fprintf(stderr, "Error: Dimension mismatch in %s\n", filepath);
                 exit(1);
             }
        }
    }

    long fpixel[3] = {1, 1, frame->local_cube_idx + 1};
    if (fits_read_pix(fptr, TDOUBLE, fpixel, expected_pixels, NULL, buffer, NULL, &status)) {
        fits_report_error(stderr, status);
    }
    fits_close_file(fptr, &status);
}

void get_timestamps(const char *stream_name, double start_time, double end_time, double **timestamps, int *n_timestamps) {
    DIR *dir = opendir(stream_name);
    if (!dir) {
        perror("Error opening stream directory for timestamp scan");
        exit(1);
    }

    struct dirent *ent;
    FileInfo *files = NULL;
    int file_count = 0;
    int capacity = 100;

    files = (FileInfo *)malloc(capacity * sizeof(FileInfo));

    while ((ent = readdir(dir)) != NULL) {
        if (strncmp(ent->d_name, stream_name, strlen(stream_name)) == 0 &&
            strstr(ent->d_name, ".txt") != NULL &&
            ent->d_name[strlen(ent->d_name) - 4] == '.') {
            if (file_count >= capacity) {
                capacity *= 2;
                files = (FileInfo *)realloc(files, capacity * sizeof(FileInfo));
            }
            strncpy(files[file_count].filename, ent->d_name, MAX_FILENAME - 1);
            files[file_count].filename[MAX_FILENAME - 1] = '\0';
            file_count++;
        }
    }
    closedir(dir);

    qsort(files, file_count, sizeof(FileInfo), compare_filenames);

    int count = 0;
    int cap = 1000;
    double *ts_list = (double*)malloc(cap * sizeof(double));

    for (int i = 0; i < file_count; i++) {
        char filepath[MAX_FILENAME * 2];
        snprintf(filepath, sizeof(filepath), "%s/%s", stream_name, files[i].filename);
        FILE *fin = fopen(filepath, "r");
        if (!fin) continue;

        char line[MAX_LINE_LENGTH];
        while (fgets(line, sizeof(line), fin)) {
            char *ptr = line;
            while (isspace(*ptr)) ptr++;
            if (*ptr == '#' || *ptr == '\0') continue;

            long col1, col2;
            double col3, col4, col5;
            if (sscanf(line, "%ld %ld %lf %lf %lf", &col1, &col2, &col3, &col4, &col5) >= 5) {
                if (col5 >= start_time && col5 <= end_time) {
                    if (count >= cap) {
                        cap *= 2;
                        ts_list = (double*)realloc(ts_list, cap * sizeof(double));
                    }
                    ts_list[count++] = col5;
                }
            }
        }
        fclose(fin);
    }
    free(files);

    *timestamps = ts_list;
    *n_timestamps = count;
    printf("Collected %d timestamps from %s\n", count, stream_name);
}

// Helper: Calculate diffs for all pairs (i, j) with j > i
void compute_all_diffs(double *data, int n_frames, int n_valid, double *timestamps, DiffResult *result, const char *stream_name) {
    if (n_frames < 2) return;

    // N*(N-1)/2 pairs
    long n_pairs = (long)n_frames * (n_frames - 1) / 2;
    result->n_frames = n_frames;
    result->n_pairs = n_pairs;
    result->times = (double*) malloc(n_pairs * sizeof(double));
    result->times_end = (double*) malloc(n_pairs * sizeof(double));
    result->diffs = (double*) malloc(n_pairs * sizeof(double));

    long pair_idx = 0;
    char diff_filename[MAX_FILENAME];
    snprintf(diff_filename, sizeof(diff_filename), "%s.diff.txt", stream_name);
    FILE *fdiff = fopen(diff_filename, "w");

    for (int i = 0; i < n_frames; i++) {
        for (int j = i + 1; j < n_frames; j++) {
            double sum_sq = 0.0;
            // Vector i start: i * n_valid
            // Vector j start: j * n_valid
            for (int p = 0; p < n_valid; p++) {
                double val_i = data[(size_t)i * n_valid + p];
                double val_j = data[(size_t)j * n_valid + p];
                double d = val_j - val_i;
                sum_sq += d * d;
            }

            result->times[pair_idx] = timestamps[i];
            result->times_end[pair_idx] = timestamps[j];
            result->diffs[pair_idx] = sum_sq;

            fprintf(fdiff, "%.6f %.6f %.9f\n", timestamps[i], timestamps[j], sum_sq);
            pair_idx++;
        }
    }
    fclose(fdiff);
    printf("Written diff file to %s (%ld pairs)\n", diff_filename, n_pairs);
}

// Returns DiffResult if do_diff is true, else NULL
DiffResult* process_stream(const char *stream_name, double start_time, double end_time, const char *output_txt_file, int do_mask_processing, double dt, double *target_ts, int n_target, int do_diff) {
    printf("Processing stream %s (Time: %.4f - %.4f, Resampling: %s, Mode: %s)\n",
            stream_name, start_time, end_time,
            (target_ts ? "Target Grid" : (dt > 0 ? "Fixed dt" : "None")),
            (do_diff ? "Diff" : "Standard"));

    DIR *dir = opendir(stream_name);
    if (!dir) {
        perror("Error opening stream directory");
        exit(1);
    }

    struct dirent *ent;
    FileInfo *files = NULL;
    int file_count = 0;
    int capacity = 100;

    files = (FileInfo *)malloc(capacity * sizeof(FileInfo));

    while ((ent = readdir(dir)) != NULL) {
        if (strncmp(ent->d_name, stream_name, strlen(stream_name)) == 0 &&
            strstr(ent->d_name, ".txt") != NULL &&
            ent->d_name[strlen(ent->d_name) - 4] == '.') {

            if (file_count >= capacity) {
                capacity *= 2;
                files = (FileInfo *)realloc(files, capacity * sizeof(FileInfo));
            }
            strncpy(files[file_count].filename, ent->d_name, MAX_FILENAME - 1);
            files[file_count].filename[MAX_FILENAME - 1] = '\0';
            file_count++;
        }
    }
    closedir(dir);

    qsort(files, file_count, sizeof(FileInfo), compare_filenames);

    FrameInfo *frames = NULL;
    int frame_count = 0;
    int frame_capacity = 1000;
    frames = (FrameInfo *)malloc(frame_capacity * sizeof(FrameInfo));

    for (int i = 0; i < file_count; i++) {
        char filepath[MAX_FILENAME * 2];
        snprintf(filepath, sizeof(filepath), "%s/%s", stream_name, files[i].filename);

        FILE *fin = fopen(filepath, "r");
        if (!fin) continue;

        char line[MAX_LINE_LENGTH];
        while (fgets(line, sizeof(line), fin)) {
            char *ptr = line;
            while (isspace(*ptr)) ptr++;
            if (*ptr == '#' || *ptr == '\0') continue;

            long col1, col2;
            double col3, col4, col5;
            int items = sscanf(line, "%ld %ld %lf %lf %lf", &col1, &col2, &col3, &col4, &col5);

            if (items >= 5) {
                double margin = (dt > 0 || target_ts) ? 0.1 : 0.0;
                if (dt > 0) margin = dt * 2.0;

                if (col5 >= start_time - margin && col5 <= end_time + margin) {
                    if (frame_count >= frame_capacity) {
                        frame_capacity *= 2;
                        frames = (FrameInfo *)realloc(frames, frame_capacity * sizeof(FrameInfo));
                    }
                    frames[frame_count].global_idx = col2;
                    frames[frame_count].timestamp = col5;
                    frames[frame_count].local_cube_idx = col1;
                    get_fits_filename(files[i].filename, frames[frame_count].fits_filename);
                    frame_count++;
                }
            }
        }
        fclose(fin);
    }
    free(files);

    if (frame_count == 0) {
        printf("No frames found in window.\n");
        free(frames);
        return NULL;
    }

    int width = 0, height = 0;
    double *mask = NULL;

    // Always load mask if diffing (to normalize summation) or if requested
    if (do_mask_processing || dt > 0 || target_ts || do_diff) {
        mask = load_mask(stream_name, &width, &height);

        if (!mask) {
            char first_fits[MAX_FILENAME * 2];
            snprintf(first_fits, sizeof(first_fits), "%s/%s", stream_name, frames[0].fits_filename);

            fitsfile *fptr;
            int status = 0;
            int naxis;
            long naxes[3];

            if (fits_open_file(&fptr, first_fits, READONLY, &status)) CHECK_STATUS(status);
            if (fits_get_img_dim(fptr, &naxis, &status)) CHECK_STATUS(status);
            if (fits_get_img_size(fptr, 3, naxes, &status)) CHECK_STATUS(status);

            width = naxes[0];
            height = naxes[1];
            fits_close_file(fptr, &status); CHECK_STATUS(status);

            printf("Created default mask %d x %d\n", width, height);
            mask = (double*) malloc(width * height * sizeof(double));
            for(int k=0; k<width*height; k++) mask[k] = 1.0;
        }
    }

    int n_pixels = width * height;
    int *valid_indices = NULL;
    int n_valid = 0;

    if (mask) {
        valid_indices = (int*) malloc(n_pixels * sizeof(int));
        for (int i = 0; i < n_pixels; i++) {
            if (mask[i] > 0.001) {
                valid_indices[n_valid++] = i;
            }
        }
        printf("Valid pixels in mask: %d / %d\n", n_valid, n_pixels);
    }

    DiffResult *result = NULL;
    if (do_diff) {
        result = (DiffResult*) malloc(sizeof(DiffResult));
        result->n_frames = 0;
        result->n_pairs = 0;
    }

    // --- Processing Logic ---

    // RESAMPLING MODE
    if (dt > 0 || target_ts) {
        int n_resampled = 0;
        if (target_ts) {
            n_resampled = n_target;
        } else {
            n_resampled = (int)ceil((end_time - start_time) / dt);
            if (n_resampled < 1) n_resampled = 1;
        }

        double *out_data = (double*) malloc((size_t)n_valid * (size_t)n_resampled * sizeof(double));
        double *timestamps = (double*) malloc(n_resampled * sizeof(double));

        FILE *fout = fopen(output_txt_file, "w");
        double *buf0 = (double*) malloc(n_pixels * sizeof(double));
        double *buf1 = (double*) malloc(n_pixels * sizeof(double));
        int idx0 = -1, idx1 = -1;

        for (int k = 0; k < n_resampled; k++) {
            double T;
            if (target_ts) T = target_ts[k];
            else { T = start_time + k * dt; if (T > end_time) break; }
            timestamps[k] = T;

            int i_left = -1, i_right = -1;
            int right_idx = -1;
            for (int i = 0; i < frame_count; i++) { if (frames[i].timestamp > T) { right_idx = i; break; } }
            if (right_idx == -1) { i_left = frame_count - 1; i_right = frame_count - 1; }
            else if (right_idx == 0) { i_left = 0; i_right = 0; }
            else { i_left = right_idx - 1; i_right = right_idx; }

            if (idx0 != i_left) { read_frame_data(stream_name, &frames[i_left], n_pixels, buf0, &width, &height); idx0 = i_left; }
            if (idx1 != i_right) { read_frame_data(stream_name, &frames[i_right], n_pixels, buf1, &width, &height); idx1 = i_right; }

            double t0 = frames[i_left].timestamp;
            double t1 = frames[i_right].timestamp;
            double alpha = 0.0;
            if (i_left != i_right && (t1 - t0) > 1e-9) alpha = (T - t0) / (t1 - t0);

            // Extract & Interpolate
            for (int p = 0; p < n_valid; p++) {
                int pix = valid_indices[p];
                double res = (buf0[pix] * mask[pix]) + alpha * ((buf1[pix] * mask[pix]) - (buf0[pix] * mask[pix]));
                out_data[(size_t)k * n_valid + p] = res;
            }

            fprintf(fout, "%d %.6f interpolated 0\n", k, T);
        }
        fclose(fout);
        free(buf0); free(buf1);

        if (do_diff) {
            compute_all_diffs(out_data, n_resampled, n_valid, timestamps, result, stream_name);
        }
        free(timestamps);

        fitsfile *fptr; int status = 0;
        char out_fits_name[MAX_FILENAME];
        snprintf(out_fits_name, sizeof(out_fits_name), "!%s_processed.fits", stream_name);
        fits_create_file(&fptr, out_fits_name, &status);
        long out_naxes[2] = {n_valid, n_resampled};
        fits_create_img(fptr, DOUBLE_IMG, 2, out_naxes, &status);
        long firstpix[2] = {1, 1};
        fits_write_pix(fptr, TDOUBLE, firstpix, (size_t)n_valid * (size_t)n_resampled, out_data, &status);
        fits_close_file(fptr, &status);
        free(out_data);

    } else {
        // RAW MODE
        // Original non-resampled logic
        FILE *fout = fopen(output_txt_file, "w");
        for (int i = 0; i < frame_count; i++) {
            if (frames[i].timestamp >= start_time && frames[i].timestamp <= end_time) {
                fprintf(fout, "%ld %.6f %s %ld\n", frames[i].global_idx, frames[i].timestamp, frames[i].fits_filename, frames[i].local_cube_idx);
            }
        }
        fclose(fout);

        // Always do mask processing if do_mask or do_diff (to get vectors)
        // If do_diff is on, we need vectors.
        if (do_mask_processing || do_diff) {
             double *out_data = (double*) malloc((size_t)n_valid * (size_t)frame_count * sizeof(double));
             double *timestamps = (double*) malloc(frame_count * sizeof(double));
             double *frame_buf = (double*) malloc(n_pixels * sizeof(double));
             char current_file[MAX_FILENAME] = "";
             fitsfile *fptr = NULL; int status = 0;

             for (int i = 0; i < frame_count; i++) {
                 timestamps[i] = frames[i].timestamp;
                 char filepath[MAX_FILENAME * 2];
                 snprintf(filepath, sizeof(filepath), "%s/%s", stream_name, frames[i].fits_filename);
                 if (strcmp(filepath, current_file) != 0) {
                     if (fptr) fits_close_file(fptr, &status);
                     fits_open_file(&fptr, filepath, READONLY, &status);
                     strcpy(current_file, filepath);
                 }
                 long fpixel[3] = {1, 1, frames[i].local_cube_idx + 1};
                 fits_read_pix(fptr, TDOUBLE, fpixel, n_pixels, NULL, frame_buf, NULL, &status);
                 for(int k=0; k<n_valid; k++) out_data[(size_t)i * n_valid + k] = frame_buf[valid_indices[k]] * mask[valid_indices[k]];
             }
             if (fptr) fits_close_file(fptr, &status);
             free(frame_buf);

             if (do_mask_processing) {
                 fitsfile *fout_fits;
                 char out_fits_name[MAX_FILENAME];
                 snprintf(out_fits_name, sizeof(out_fits_name), "!%s_processed.fits", stream_name);
                 fits_create_file(&fout_fits, out_fits_name, &status);
                 long out_naxes[2] = {n_valid, frame_count};
                 fits_create_img(fout_fits, DOUBLE_IMG, 2, out_naxes, &status);
                 long firstpix[2] = {1, 1};
                 fits_write_pix(fout_fits, TDOUBLE, firstpix, (size_t)n_valid * (size_t)frame_count, out_data, &status);
                 fits_close_file(fout_fits, &status);
             }

             if (do_diff) {
                 compute_all_diffs(out_data, frame_count, n_valid, timestamps, result, stream_name);
             }

             free(timestamps);
             free(out_data);
        }
    }

    if (mask) free(mask);
    if (valid_indices) free(valid_indices);
    free(frames);

    return result;
}

int main(int argc, char *argv[]) {
    int do_mask = 0;
    double dt = 0.0;
    int resample_A = 0;
    int resample_B = 0;
    int ndtmax = -1;
    int do_diff = 0;
    int arg_idx = 1;

    while (arg_idx < argc && argv[arg_idx][0] == '-') {
        if (strcmp(argv[arg_idx], "-m") == 0) {
            do_mask = 1;
            arg_idx++;
        } else if (strcmp(argv[arg_idx], "-dt") == 0) {
            dt = atof(argv[++arg_idx]);
            arg_idx++;
        } else if (strcmp(argv[arg_idx], "-A") == 0) {
            resample_A = 1;
            arg_idx++;
        } else if (strcmp(argv[arg_idx], "-B") == 0) {
            resample_B = 1;
            arg_idx++;
        } else if (strcmp(argv[arg_idx], "-odiff") == 0) {
            ndtmax = atoi(argv[++arg_idx]);
            do_diff = 1;
            arg_idx++;
        } else {
            fprintf(stderr, "Unknown option %s\n", argv[arg_idx]);
            return 1;
        }
    }

    if (resample_A && resample_B) { fprintf(stderr, "Error: Cannot use both -A and -B.\n"); return 1; }

    if (argc - arg_idx != 4) {
        fprintf(stderr, "Usage: %s [-m] [-dt val | -A | -B] [-odiff ndtmax] <streamA> <streamB> <timestart> <timeend>\n", argv[0]);
        return 1;
    }

    const char *streamA = argv[arg_idx];
    const char *streamB = argv[arg_idx+1];
    double timestart = atof(argv[arg_idx+2]);
    double timeend = atof(argv[arg_idx+3]);

    char outA[256], outB[256];
    snprintf(outA, sizeof(outA), "%s.flist.txt", streamA);
    snprintf(outB, sizeof(outB), "%s.flist.txt", streamB);

    double *timestamps = NULL;
    int n_timestamps = 0;

    if (resample_A) get_timestamps(streamA, timestart, timeend, &timestamps, &n_timestamps);
    if (resample_B) get_timestamps(streamB, timestart, timeend, &timestamps, &n_timestamps);

    // Call process_stream for both, passing all flags
    DiffResult *resA = process_stream(streamA, timestart, timeend, outA, do_mask, dt, timestamps, n_timestamps, do_diff);
    DiffResult *resB = process_stream(streamB, timestart, timeend, outB, do_mask, dt, timestamps, n_timestamps, do_diff);

    if (timestamps) free(timestamps);

    if (do_diff && resA && resB) {
        for (int ndt = -ndtmax; ndt <= ndtmax; ndt++) {
            char outfile[MAX_FILENAME];
            snprintf(outfile, sizeof(outfile), "diff_ndt_%d.txt", ndt);
            FILE *f = fopen(outfile, "w");

            // Loop over pairs in A
            // resA->times has N_pairs elements.
            // pair_idx = i*(2N - 1 - i)/2 + (j - i - 1)
            // But we just iterated flatly.
            // How do we match "Line N" and "Line N+ndt"?
            // We assume the user implies "Line N of the FILE", which matches "pair N of the stored list".
            // Since both files are generated by the same logic (all unique pairs, sorted by i then j),
            // matching Line N to Line N+ndt matches the N-th pair of A to the (N+ndt)-th pair of B.

            long N_pairs_A = resA->n_pairs;
            long N_pairs_B = resB->n_pairs;

            for (long N = 0; N < N_pairs_A; N++) {
                long K = N + ndt;
                if (K >= 0 && K < N_pairs_B) {
                    fprintf(f, "%.6f %.6f %.9f %.6f %.6f %.9f\n",
                            resA->times[N], resA->times_end[N], resA->diffs[N],
                            resB->times[K], resB->times_end[K], resB->diffs[K]);
                }
            }
            fclose(f);
        }
        if (resA) { free(resA->times); free(resA->times_end); free(resA->diffs); free(resA); }
        if (resB) { free(resB->times); free(resB->times_end); free(resB->diffs); free(resB); }
    }

    return 0;
}
