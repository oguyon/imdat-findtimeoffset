#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <ctype.h>

// Max line length for reading input files
#define MAX_LINE_LENGTH 1024
#define MAX_FILENAME 512

// Structure to hold file info for sorting
typedef struct {
    char filename[MAX_FILENAME];
    // Potentially extract time from filename if needed for sorting,
    // but alphabetical sorting of filenames "streamA_HH:MM:SS.SSSSSS"
    // works for sorting by time within a day (and we assume single day or standard sorting).
} FileInfo;

int compare_filenames(const void *a, const void *b) {
    const FileInfo *fa = (const FileInfo *)a;
    const FileInfo *fb = (const FileInfo *)b;
    return strcmp(fa->filename, fb->filename);
}

// Function to replace .txt with .fits
void get_fits_filename(const char *txt_filename, char *fits_filename) {
    strcpy(fits_filename, txt_filename);
    char *dot = strrchr(fits_filename, '.');
    if (dot) {
        strcpy(dot, ".fits");
    } else {
        strcat(fits_filename, ".fits");
    }
}

void process_stream(const char *stream_name, double start_time, double end_time, const char *output_filename) {
    printf("Processing stream %s into %s (Time: %.4f - %.4f)\n", stream_name, output_filename, start_time, end_time);

    DIR *dir = opendir(stream_name);
    if (!dir) {
        perror("Error opening stream directory");
        return; // Or exit
    }

    // Collect files
    struct dirent *ent;
    FileInfo *files = NULL;
    int file_count = 0;
    int capacity = 100;

    files = (FileInfo *)malloc(capacity * sizeof(FileInfo));
    if (!files) {
        perror("Memory allocation failed");
        exit(1);
    }

    char prefix[256];
    snprintf(prefix, sizeof(prefix), "%s", stream_name);

    while ((ent = readdir(dir)) != NULL) {
        // Filter: starts with stream_name and ends with .txt
        if (strncmp(ent->d_name, stream_name, strlen(stream_name)) == 0 &&
            strstr(ent->d_name, ".txt") != NULL &&
            ent->d_name[strlen(ent->d_name) - 4] == '.') {

            if (file_count >= capacity) {
                capacity *= 2;
                files = (FileInfo *)realloc(files, capacity * sizeof(FileInfo));
                if (!files) {
                    perror("Memory reallocation failed");
                    exit(1);
                }
            }
            strncpy(files[file_count].filename, ent->d_name, MAX_FILENAME - 1);
            files[file_count].filename[MAX_FILENAME - 1] = '\0';
            file_count++;
        }
    }
    closedir(dir);

    // Sort files
    qsort(files, file_count, sizeof(FileInfo), compare_filenames);

    FILE *fout = fopen(output_filename, "w");
    if (!fout) {
        perror("Error opening output file");
        free(files);
        exit(1);
    }

    for (int i = 0; i < file_count; i++) {
        char filepath[MAX_FILENAME * 2];
        snprintf(filepath, sizeof(filepath), "%s/%s", stream_name, files[i].filename);

        FILE *fin = fopen(filepath, "r");
        if (!fin) {
            fprintf(stderr, "Warning: Could not open %s\n", filepath);
            continue;
        }

        char line[MAX_LINE_LENGTH];
        while (fgets(line, sizeof(line), fin)) {
            // Skip comments and empty lines
            char *ptr = line;
            while (isspace(*ptr)) ptr++;
            if (*ptr == '#' || *ptr == '\0') continue;

            // Parse columns
            // col1: frame index (within file)
            // col2: Main index
            // col3, col4...
            // col5: Absolute time (acquisition)

            long col1_cube_idx;
            long col2_main_idx;
            double col3, col4, col5_time;

            // Format: %ld %ld %lf %lf %lf ...
            // We need to parse up to col 5.
            // Using sscanf
            // Note: Use temporary variables for ignored columns if needed, but we can just scan them.

            int items = sscanf(line, "%ld %ld %lf %lf %lf", &col1_cube_idx, &col2_main_idx, &col3, &col4, &col5_time);

            if (items >= 5) {
                if (col5_time >= start_time && col5_time <= end_time) {
                    // Match!
                    char fits_name[MAX_FILENAME];
                    get_fits_filename(files[i].filename, fits_name);

                    // Output format:
                    // col 1: frame index (extracted from input col 2) -> col2_main_idx
                    // col 2: acquisition time (input col 5) -> col5_time
                    // col 3: FITS file name -> fits_name
                    // col 4: frame index within FITS file (input col 1) -> col1_cube_idx

                    fprintf(fout, "%ld %.6f %s %ld\n", col2_main_idx, col5_time, fits_name, col1_cube_idx);
                }
            }
        }
        fclose(fin);
    }

    fclose(fout);
    free(files);
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <streamA> <streamB> <timestart> <timeend>\n", argv[0]);
        return 1;
    }

    const char *streamA = argv[1];
    const char *streamB = argv[2];
    double timestart = atof(argv[3]);
    double timeend = atof(argv[4]);

    char outA[256], outB[256];
    snprintf(outA, sizeof(outA), "%s.flist.txt", streamA);
    snprintf(outB, sizeof(outB), "%s.flist.txt", streamB);

    process_stream(streamA, timestart, timeend, outA);
    process_stream(streamB, timestart, timeend, outB);

    return 0;
}
