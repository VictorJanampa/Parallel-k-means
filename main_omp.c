#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

#include "image_io.h"
#include "segmentation.h"

#define DEFAULT_N_CLUSTS 4
#define DEFAULT_MAX_ITERS 150
#define DEFAULT_N_THREADS 2
#define DEFAULT_OUT_PATH "result.jpg"

double get_time();

int main(int argc, char **argv)
{
    char *in_path = NULL;
    char *out_path = DEFAULT_OUT_PATH;
    byte_t *data;
    int width, height, n_ch;
    int n_clus = DEFAULT_N_CLUSTS;
    int n_iters = DEFAULT_MAX_ITERS;
    int n_threads = DEFAULT_N_THREADS;
    int seed = time(NULL);
    double sse, start_time, exec_time;



    char optchar;
    while ((optchar = getopt(argc, argv, "k:m:o:s:t:h")) != -1) {
        switch (optchar) {
            case 'k':
                n_clus = strtol(optarg, NULL, 10);
                break;
            case 'm':
                n_iters = strtol(optarg, NULL, 10);
                break;
            case 'o':
                out_path = optarg;
                break;
            case 's':
                seed = strtol(optarg, NULL, 10);
                break;
            case 't':
                n_threads = strtol(optarg, NULL, 10);
                break;
            case 'h':
            default:
                exit(EXIT_FAILURE);
                break;
        }
    }

    in_path = argv[optind];

    if (in_path == NULL) {
        exit(EXIT_FAILURE);
    }

    if (n_clus < 2) {
        fprintf(stderr, "INPUT ERROR: << Invalid number of clusters >> \n");
        exit(EXIT_FAILURE);
    }

    if (n_iters < 1) {
        fprintf(stderr, "INPUT ERROR: << Invalid maximum number of iterations >> \n");
        exit(EXIT_FAILURE);
    }

    if (n_threads < 2) {
        fprintf(stderr, "INPUT ERROR: << Invalid number of threads >> \n");
        exit(EXIT_FAILURE);
    }

    srand(seed);

    data = img_load(in_path, &width, &height, &n_ch);

    start_time = get_time();
    kmeans_segm_omp(data, width, height, n_ch, n_clus, &n_iters, &sse, n_threads);
    exec_time = get_time() - start_time;
    fprintf(stderr, "TIME: %f\n",exec_time);
    img_save(out_path, data, width, height, n_ch);

    free(data);

    return EXIT_SUCCESS;
}

double get_time()
{
    struct timeval timecheck;

    gettimeofday(&timecheck, NULL);

    return timecheck.tv_sec + timecheck.tv_usec / 1000000.0;
}
