#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

void print_usage(char **argv) {
    fprintf(stderr, "Usage: %s [-d] [-c chunk_size] one.fastq.gz two.fastq.gz <combined.fastq>\n", argv[0]);
}


int main(int argc, char **argv) {
    bool debug_mode = false;
    int chunk_size = 0;
    char *file_out_orig = NULL;
    int option = 0;

    opterr = 0;

    while ((option = getopt(argc, argv, "c:do:")) != -1) {
        switch (option) {
            case 'c':
                chunk_size = atoi(optarg);
                break;
            case 'd':
                debug_mode = true;
                break;
            case 'o':
                asprintf(&file_out_orig, "%s", optarg);
                break;
            default:
                print_usage(argv);
                exit(EXIT_FAILURE);
            }
    }

    char* file_1 = NULL;
    char* file_2 = NULL;
    char* file_out = NULL;

    if ((argc - optind) == 2) {
        file_1 = argv[optind];
        file_2 = argv[optind+1];
    } else {
        print_usage(argv);
        exit(EXIT_FAILURE);
    }


    if ((chunk_size != 0) && (!file_out_orig)) {
        fprintf(stderr, "-c cannot be specified when using stdout");
        exit(EXIT_FAILURE);
    }

    if (debug_mode) {
        fprintf(stderr, "R1: %s\n", file_1);
        fprintf(stderr, "R2: %s\n", file_2);
        if (file_out_orig) {
            fprintf(stderr, "Output: %s\n", file_out_orig);
        }
        fprintf(stderr, "Chunk size: %d\n", chunk_size);
    }

    // timers for debugging
    clock_t t_start;
    clock_t t_chunk;
    clock_t t_temp;

    t_start = clock();
    t_chunk = clock();

    /* setup R1 */
    gzFile fp_1;
    kseq_t *seq_1;
    FILE *f_1 = fopen(file_1, "rb");
    if (f_1 == NULL) {
        fprintf(stderr, "Can't open: %s\n", file_1);
        exit(EXIT_FAILURE);
    }
    fp_1 = gzdopen(fileno(f_1), "r");
    seq_1 = kseq_init(fp_1);

    /* setup R2 */
    gzFile fp_2;
    kseq_t *seq_2;
    FILE *f_2 = fopen(file_2, "rb");
    if (f_2 == NULL) {
        fprintf(stderr, "Can't open: %s\n", file_2);
        exit(EXIT_FAILURE);
    }
    fp_2 = gzdopen(fileno(f_2), "r");
    seq_2 = kseq_init(fp_2);

    FILE *f_out = stdout;
    if (file_out_orig) {
        if (chunk_size > 0) {
            asprintf(&file_out, "%s_0.fastq", file_out_orig);
        } else {
            asprintf(&file_out, "%s", file_out_orig);
        }
        f_out = fopen(file_out, "w");
    }

    int debug_chunk = 1000000;

    if (chunk_size > 0) {
        debug_chunk = chunk_size;
    }

    /*
    CR = R1.sequence[:16]
    CY = R1.quality[:16]
    UR = R1.sequence[16:]
    UY = R1.quality[16:]
    */

    char CR[17];
    char CY[17];
    char UR[11];
    char UY[11];

    //name, comment, seq, qual; )
    int counter = 0;
    int file_counter = 0;
    while (kseq_read(seq_1) >= 0) {
        kseq_read(seq_2);

        memcpy(CR, &seq_1->seq.s[0], 16);
        CR[16] = '\0';

        memcpy(CY, &seq_1->qual.s[0], 16);
        CY[16] = '\0';

        memcpy(UR, &seq_1->seq.s[16], 10);
        UR[10] = '\0';

        memcpy(UY, &seq_1->qual.s[16], 10);
        UY[10] = '\0';

        fprintf(f_out, "@%s|||BC||||||QT||||||CR|||%s|||CY|||%s|||UR|||%s|||UY|||%s %s\n", seq_1->name.s, CR, CY, UR, UY, seq_2->comment.s);
        fprintf(f_out, "%s\n+\n%s\n", seq_2->seq.s, seq_2->qual.s);

        counter = counter + 1;

        if ((debug_mode) && ((counter % debug_chunk) == 0)) {
            t_temp = clock() - t_start;
            double time_total = ((double)t_temp)/CLOCKS_PER_SEC; // in seconds

            t_temp = clock() - t_chunk;
            double time_chunk = ((double)t_temp)/CLOCKS_PER_SEC; // in seconds

            fprintf(stderr, "Reads: %d, Chunk Time: %f, Total Time: %f\n", counter, time_chunk, time_total);

            t_chunk = clock();
        }

        if ((chunk_size > 0) && (counter % chunk_size) == 0) {
            fclose(f_out);
            file_counter = file_counter + 1;

            if (debug_mode) {
                fprintf(stderr, "Generated: %s\n", file_out);
            }

            asprintf(&file_out, "%s_%d.fastq", file_out_orig, file_counter);
            
            if (debug_mode) {
                fprintf(stderr, "Generating: %s\n", file_out);
            }

            f_out = fopen(file_out, "w");
        }
    }

    kseq_destroy(seq_2);
    fclose(f_2);
    gzclose(fp_2);

    kseq_destroy(seq_1);
    fclose(f_1);
    gzclose(fp_1);

    fclose(f_out);

    if (debug_mode) {
        fprintf(stderr, "Generated: %s\n", file_out);
        t_temp = clock() - t_start;
        double time_total = ((double)t_temp)/CLOCKS_PER_SEC; // in seconds

        fprintf(stderr, "Done. Reads: %d, Total Time: %f\n", counter, time_total);
    }

    return 0;
}