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

// Load mask or return NULL if not found. Sets width/height.
double* load_mask(const char *stream_name, int *width, int *height) {
    char mask_filename[MAX_FILENAME];
    snprintf(mask_filename, sizeof(mask_filename), "%s.mask.fits", stream_name);

    if (access(mask_filename, F_OK) != 0) {
        // printf("No mask file found (%s). Will use default all-ones mask.\n", mask_filename);
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

// Helper to read a frame into buffer
void read_frame_data(const char *stream_name, FrameInfo *frame, int expected_pixels, double *buffer, int *width_chk, int *height_chk) {
    char filepath[MAX_FILENAME * 2];
    snprintf(filepath, sizeof(filepath), "%s/%s", stream_name, frame->fits_filename);

    fitsfile *fptr;
    int status = 0;

    if (fits_open_file(&fptr, filepath, READONLY, &status)) {
        // Warning?
        fits_report_error(stderr, status);
        memset(buffer, 0, expected_pixels * sizeof(double));
        return;
    }

    // Check dimensions if pointers provided
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

// Collect timestamps from a stream
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

    // Scan files for timestamps
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

void process_stream(const char *stream_name, double start_time, double end_time, const char *output_txt_file, int do_mask_processing, double dt, double *target_ts, int n_target) {
    printf("Processing stream %s (Time: %.4f - %.4f, Resampling: %s)\n",
            stream_name, start_time, end_time,
            (target_ts ? "Target Grid" : (dt > 0 ? "Fixed dt" : "None")));

    DIR *dir = opendir(stream_name);
    if (!dir) {
        perror("Error opening stream directory");
        exit(1);
    }

    // Collect files
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
                // If resampling (dt > 0 or target_ts), widen window slightly
                double margin = (dt > 0 || target_ts) ? 0.1 : 0.0; // Arbitrary small margin if dt is unknown, or rely on logic
                if (dt > 0) margin = dt * 2.0;

                // If target_ts provided, we might need frames covering the full range of target_ts
                // The passed start_time/end_time should cover it.

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
        return;
    }

    // Determine dimensions from first frame/mask
    int width = 0, height = 0;
    double *mask = NULL;

    if (do_mask_processing || dt > 0 || target_ts) {
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

    // Processing Logic

    if (dt > 0 || target_ts) {
        // RESAMPLING MODE
        int n_resampled = 0;

        // Use target_ts if provided, otherwise generate grid from dt
        if (target_ts) {
            n_resampled = n_target;
        } else {
            n_resampled = (int)ceil((end_time - start_time) / dt);
            if (n_resampled < 1) n_resampled = 1;
        }

        double *out_data = (double*) malloc((size_t)n_valid * (size_t)n_resampled * sizeof(double));

        FILE *fout = fopen(output_txt_file, "w");
        if (!fout) { perror("Error opening output file"); exit(1); }

        double *buf0 = (double*) malloc(n_pixels * sizeof(double));
        double *buf1 = (double*) malloc(n_pixels * sizeof(double));
        int idx0 = -1;
        int idx1 = -1;

        for (int k = 0; k < n_resampled; k++) {
            double T;
            if (target_ts) {
                T = target_ts[k];
            } else {
                T = start_time + k * dt;
                if (T > end_time) break;
            }

            int i_left = -1;
            int i_right = -1;

            // Find bracket for T
            int right_idx = -1;
            for (int i = 0; i < frame_count; i++) {
                if (frames[i].timestamp > T) {
                    right_idx = i;
                    break;
                }
            }

            if (right_idx == -1) {
                i_left = frame_count - 1;
                i_right = frame_count - 1;
            } else if (right_idx == 0) {
                i_left = 0;
                i_right = 0;
            } else {
                i_left = right_idx - 1;
                i_right = right_idx;
            }

            if (idx0 != i_left) {
                read_frame_data(stream_name, &frames[i_left], n_pixels, buf0, &width, &height);
                idx0 = i_left;
            }
            if (idx1 != i_right) {
                read_frame_data(stream_name, &frames[i_right], n_pixels, buf1, &width, &height);
                idx1 = i_right;
            }

            double t0 = frames[i_left].timestamp;
            double t1 = frames[i_right].timestamp;
            double alpha = 0.0;

            if (i_left != i_right && (t1 - t0) > 1e-9) {
                alpha = (T - t0) / (t1 - t0);
            }

            for (int p = 0; p < n_valid; p++) {
                int pix = valid_indices[p];
                double val0 = buf0[pix] * mask[pix];
                double val1 = buf1[pix] * mask[pix];
                double res = val0 + alpha * (val1 - val0);
                out_data[(size_t)k * n_valid + p] = res;
            }

            fprintf(fout, "%d %.6f interpolated 0\n", k, T);
        }

        fclose(fout);
        free(buf0);
        free(buf1);

        fitsfile *fptr;
        int status = 0;
        char out_fits_name[MAX_FILENAME];
        snprintf(out_fits_name, sizeof(out_fits_name), "!%s_processed.fits", stream_name);
        if (fits_create_file(&fptr, out_fits_name, &status)) CHECK_STATUS(status);
        long out_naxes[2] = {n_valid, n_resampled};
        if (fits_create_img(fptr, DOUBLE_IMG, 2, out_naxes, &status)) CHECK_STATUS(status);
        long firstpix[2] = {1, 1};
        if (fits_write_pix(fptr, TDOUBLE, firstpix, (size_t)n_valid * (size_t)n_resampled, out_data, &status)) CHECK_STATUS(status);
        fits_close_file(fptr, &status); CHECK_STATUS(status);

        printf("Written resampled data to %s (%ld x %d)\n", out_fits_name, out_naxes[0], (int)out_naxes[1]);
        free(out_data);

    } else {
        // ORIGINAL LOGIC (No Resampling)
        FILE *fout = fopen(output_txt_file, "w");
        if (!fout) { perror("Error opening output file"); exit(1); }
        for (int i = 0; i < frame_count; i++) {
            if (frames[i].timestamp >= start_time && frames[i].timestamp <= end_time) {
                fprintf(fout, "%ld %.6f %s %ld\n",
                        frames[i].global_idx,
                        frames[i].timestamp,
                        frames[i].fits_filename,
                        frames[i].local_cube_idx);
            }
        }
        fclose(fout);
        printf("Written %d frames to %s\n", frame_count, output_txt_file);

        if (do_mask_processing && frame_count > 0) {
            double *out_data = (double*) malloc((size_t)n_valid * (size_t)frame_count * sizeof(double));
            double *frame_buf = (double*) malloc(n_pixels * sizeof(double));
            char current_file[MAX_FILENAME] = "";
            fitsfile *fptr = NULL;
            int status = 0;

            for (int i = 0; i < frame_count; i++) {
                if (frames[i].timestamp < start_time || frames[i].timestamp > end_time) continue;

                char filepath[MAX_FILENAME * 2];
                snprintf(filepath, sizeof(filepath), "%s/%s", stream_name, frames[i].fits_filename);

                if (strcmp(filepath, current_file) != 0) {
                    if (fptr) { fits_close_file(fptr, &status); status = 0; }
                    if (fits_open_file(&fptr, filepath, READONLY, &status)) {
                         fits_report_error(stderr, status); status = 0; continue;
                    }
                    int naxis; long naxes[3];
                    if (fits_get_img_dim(fptr, &naxis, &status) == 0 && fits_get_img_size(fptr, 3, naxes, &status) == 0) {
                        if (naxes[0] != width || naxes[1] != height) {
                            fprintf(stderr, "Error: Dimensions mismatch in %s\n", filepath); exit(1);
                        }
                    }
                    strcpy(current_file, filepath);
                }
                long fpixel[3] = {1, 1, frames[i].local_cube_idx + 1};
                if (fits_read_pix(fptr, TDOUBLE, fpixel, n_pixels, NULL, frame_buf, NULL, &status)) {
                     fits_report_error(stderr, status); status = 0; continue;
                }

                for (int k = 0; k < n_valid; k++) {
                    int pix_idx = valid_indices[k];
                    double val = frame_buf[pix_idx] * mask[pix_idx];
                    out_data[(size_t)i * n_valid + k] = val;
                }
            }
            if (fptr) fits_close_file(fptr, &status);

            char out_fits_name[MAX_FILENAME];
            snprintf(out_fits_name, sizeof(out_fits_name), "!%s_processed.fits", stream_name);
            if (fits_create_file(&fptr, out_fits_name, &status)) CHECK_STATUS(status);
            long out_naxes[2] = {n_valid, frame_count};
            if (fits_create_img(fptr, DOUBLE_IMG, 2, out_naxes, &status)) CHECK_STATUS(status);
            long firstpix[2] = {1, 1};
            if (fits_write_pix(fptr, TDOUBLE, firstpix, (size_t)n_valid * (size_t)frame_count, out_data, &status)) CHECK_STATUS(status);
            fits_close_file(fptr, &status); CHECK_STATUS(status);

            printf("Written processed data to %s (%ld x %d)\n", out_fits_name, out_naxes[0], (int)out_naxes[1]);
            free(out_data);
            free(frame_buf);
        }
    }

    if (mask) free(mask);
    if (valid_indices) free(valid_indices);
    free(frames);
}

int main(int argc, char *argv[]) {
    int do_mask = 0;
    double dt = 0.0;
    int resample_A = 0;
    int resample_B = 0;
    int arg_idx = 1;

    while (arg_idx < argc && argv[arg_idx][0] == '-') {
        if (strcmp(argv[arg_idx], "-m") == 0) {
            do_mask = 1;
            arg_idx++;
        } else if (strcmp(argv[arg_idx], "-dt") == 0) {
            if (arg_idx + 1 < argc) {
                dt = atof(argv[++arg_idx]);
                arg_idx++;
            } else {
                fprintf(stderr, "Error: -dt requires a value.\n");
                return 1;
            }
        } else if (strcmp(argv[arg_idx], "-A") == 0) {
            resample_A = 1;
            arg_idx++;
        } else if (strcmp(argv[arg_idx], "-B") == 0) {
            resample_B = 1;
            arg_idx++;
        } else {
            fprintf(stderr, "Unknown option %s\n", argv[arg_idx]);
            return 1;
        }
    }

    if (resample_A && resample_B) {
        fprintf(stderr, "Error: Cannot use both -A and -B.\n");
        return 1;
    }

    if (argc - arg_idx != 4) {
        fprintf(stderr, "Usage: %s [-m] [-dt val | -A | -B] <streamA> <streamB> <timestart> <timeend>\n", argv[0]);
        return 1;
    }

    const char *streamA = argv[arg_idx];
    const char *streamB = argv[arg_idx+1];
    double timestart = atof(argv[arg_idx+2]);
    double timeend = atof(argv[arg_idx+3]);

    char outA[256], outB[256];
    snprintf(outA, sizeof(outA), "%s.flist.txt", streamA);
    snprintf(outB, sizeof(outB), "%s.flist.txt", streamB);

    // Logic for -A or -B
    double *timestamps = NULL;
    int n_timestamps = 0;

    if (resample_A) {
        get_timestamps(streamA, timestart, timeend, &timestamps, &n_timestamps);
        // Process A naturally (or with target=timestamps, effectively identity but safer to ensure exact match)
        // Actually, if we pass target timestamps matching A's original, interpolation returns A's values (mostly).
        // Let's pass the timestamps to A too, to ensure "resampled to A" (which is identity for A).
        process_stream(streamA, timestart, timeend, outA, do_mask, 0.0, timestamps, n_timestamps);
        process_stream(streamB, timestart, timeend, outB, do_mask, 0.0, timestamps, n_timestamps);
    } else if (resample_B) {
        get_timestamps(streamB, timestart, timeend, &timestamps, &n_timestamps);
        process_stream(streamA, timestart, timeend, outA, do_mask, 0.0, timestamps, n_timestamps);
        process_stream(streamB, timestart, timeend, outB, do_mask, 0.0, timestamps, n_timestamps);
    } else {
        process_stream(streamA, timestart, timeend, outA, do_mask, dt, NULL, 0);
        process_stream(streamB, timestart, timeend, outB, do_mask, dt, NULL, 0);
    }

    if (timestamps) free(timestamps);

    return 0;
}
