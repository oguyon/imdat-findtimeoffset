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
        printf("No mask file found (%s). Will use default all-ones mask.\n", mask_filename);
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

void process_stream(const char *stream_name, double start_time, double end_time, const char *output_txt_file, int do_mask_processing) {
    printf("Processing stream %s (Time: %.4f - %.4f)\n", stream_name, start_time, end_time);

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

    // Collect valid frames
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
            // Scan 5 items
            int items = sscanf(line, "%ld %ld %lf %lf %lf", &col1, &col2, &col3, &col4, &col5);

            if (items >= 5) {
                if (col5 >= start_time && col5 <= end_time) {
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

    // Write Text Output
    FILE *fout = fopen(output_txt_file, "w");
    if (!fout) {
        perror("Error opening output file");
        exit(1);
    }
    for (int i = 0; i < frame_count; i++) {
        fprintf(fout, "%ld %.6f %s %ld\n",
                frames[i].global_idx,
                frames[i].timestamp,
                frames[i].fits_filename,
                frames[i].local_cube_idx);
    }
    fclose(fout);
    printf("Written %d frames to %s\n", frame_count, output_txt_file);

    // Data Processing
    if (do_mask_processing && frame_count > 0) {
        int width = 0, height = 0;
        double *mask = load_mask(stream_name, &width, &height);

        // If mask not found, we need to get dimensions from first file
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

        // Identify valid pixels (mask > 0.001)
        int n_pixels = width * height;
        int *valid_indices = (int*) malloc(n_pixels * sizeof(int));
        int n_valid = 0;
        for (int i = 0; i < n_pixels; i++) {
            if (mask[i] > 0.001) {
                valid_indices[n_valid++] = i;
            }
        }

        printf("Valid pixels in mask: %d / %d\n", n_valid, n_pixels);

        // Allocate Output Data: N_valid x N_frames
        double *out_data = (double*) malloc((size_t)n_valid * (size_t)frame_count * sizeof(double));
        if (!out_data) {
             fprintf(stderr, "Memory allocation failed for output data\n");
             exit(1);
        }

        char current_file[MAX_FILENAME] = "";
        fitsfile *fptr = NULL;
        int status = 0;
        double *frame_buf = (double*) malloc(n_pixels * sizeof(double));

        for (int i = 0; i < frame_count; i++) {
            char filepath[MAX_FILENAME * 2];
            snprintf(filepath, sizeof(filepath), "%s/%s", stream_name, frames[i].fits_filename);

            if (strcmp(filepath, current_file) != 0) {
                if (fptr) {
                    fits_close_file(fptr, &status);
                    status = 0;
                }
                if (fits_open_file(&fptr, filepath, READONLY, &status)) {
                    fits_report_error(stderr, status);
                    status = 0;
                    continue;
                }

                // Verify Dimensions
                int file_naxis;
                long file_naxes[3];
                if (fits_get_img_dim(fptr, &file_naxis, &status)) CHECK_STATUS(status);
                if (fits_get_img_size(fptr, 3, file_naxes, &status)) CHECK_STATUS(status);

                if (file_naxes[0] != width || file_naxes[1] != height) {
                    fprintf(stderr, "Error: Dimensions of %s (%ldx%ld) do not match mask dimensions (%dx%d)\n",
                            filepath, file_naxes[0], file_naxes[1], width, height);
                    exit(1);
                }

                strcpy(current_file, filepath);
            }

            // Read Frame
            long fpixel[3] = {1, 1, frames[i].local_cube_idx + 1};
            if (fits_read_pix(fptr, TDOUBLE, fpixel, n_pixels, NULL, frame_buf, NULL, &status)) {
                 fits_report_error(stderr, status);
                 status = 0;
                 continue;
            }

            // Extract and Mask
            for (int k = 0; k < n_valid; k++) {
                int pix_idx = valid_indices[k];
                double val = frame_buf[pix_idx] * mask[pix_idx];
                out_data[(size_t)i * n_valid + k] = val;
            }
        }

        if (fptr) fits_close_file(fptr, &status);
        free(frame_buf);
        free(mask);
        free(valid_indices);

        // Write Output FITS
        char out_fits_name[MAX_FILENAME];
        snprintf(out_fits_name, sizeof(out_fits_name), "%s_processed.fits", stream_name);

        if (fits_create_file(&fptr, out_fits_name, &status)) CHECK_STATUS(status);

        // Create Image
        // NAXIS1 = N_valid, NAXIS2 = frame_count
        long out_naxes[2] = {n_valid, frame_count};
        if (fits_create_img(fptr, DOUBLE_IMG, 2, out_naxes, &status)) CHECK_STATUS(status);

        long firstpix[2] = {1, 1};
        size_t total_elements = (size_t)n_valid * (size_t)frame_count;
        if (fits_write_pix(fptr, TDOUBLE, firstpix, total_elements, out_data, &status)) CHECK_STATUS(status);

        fits_close_file(fptr, &status); CHECK_STATUS(status);

        printf("Written processed data to %s (%ld x %d)\n", out_fits_name, out_naxes[0], (int)out_naxes[1]);
        free(out_data);
    }

    free(frames);
}

int main(int argc, char *argv[]) {
    int do_mask = 0;
    int arg_idx = 1;

    if (argc > 1 && strcmp(argv[1], "-m") == 0) {
        do_mask = 1;
        arg_idx++;
    }

    if (argc - arg_idx != 4) {
        fprintf(stderr, "Usage: %s [-m] <streamA> <streamB> <timestart> <timeend>\n", argv[0]);
        return 1;
    }

    const char *streamA = argv[arg_idx];
    const char *streamB = argv[arg_idx+1];
    double timestart = atof(argv[arg_idx+2]);
    double timeend = atof(argv[arg_idx+3]);

    char outA[256], outB[256];
    snprintf(outA, sizeof(outA), "%s.flist.txt", streamA);
    snprintf(outB, sizeof(outB), "%s.flist.txt", streamB);

    process_stream(streamA, timestart, timeend, outA, do_mask);
    process_stream(streamB, timestart, timeend, outB, do_mask);

    return 0;
}
