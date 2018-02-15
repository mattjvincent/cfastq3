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
    fprintf(stderr, "Usage: %s [-d] [-c chunk_size] [-o output.fastq] i1.fastq.gz r1.fastq.gz r2.fastq.gz\n", argv[0]);
}


int main(int argc, char **argv) {
    bool debug_mode = false;
    int chunk_size = 0;
    char *filename_output_original = NULL;
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
                asprintf(&filename_output_original, "%s", optarg);
                break;
            default:
                print_usage(argv);
                exit(EXIT_FAILURE);
            }
    }

    /* sample barcode */
    char* filename_i1 = NULL;

    /* cell barcode + UMI */
    char* filename_r1 = NULL;

    /* reads */
    char* filename_r2 = NULL;

    /* output file */
    char* filename_output = NULL;

    if ((argc - optind) == 3) {
        filename_i1 = argv[optind];
        filename_r1 = argv[optind+1];
        filename_r2 = argv[optind+2];
    } else {
        print_usage(argv);
        exit(EXIT_FAILURE);
    }


    if ((chunk_size != 0) && (!filename_output_original)) {
        fprintf(stderr, "-c cannot be specified when using stdout");
        exit(EXIT_FAILURE);
    }

    if (debug_mode) {
        fprintf(stderr, "I1: %s\n", filename_i1);
        fprintf(stderr, "R1: %s\n", filename_r1);
        fprintf(stderr, "R2: %s\n", filename_r2);
        if (filename_output_original) {
            fprintf(stderr, "Output: %s\n", filename_output_original);
        }
        fprintf(stderr, "Chunk size: %d\n", chunk_size);
    }

    // timers for debugging
    clock_t t_start;
    clock_t t_chunk;
    clock_t t_temp;

    t_start = clock();
    t_chunk = clock();


    /* setup I1 */
    gzFile fp_i1;
    kseq_t *seq_i1;
    FILE *f_i1 = fopen(filename_i1, "rb");
    if (f_i1 == NULL) {
        fprintf(stderr, "Can't open: %s\n", filename_i1);
        exit(EXIT_FAILURE);
    }
    fp_i1 = gzdopen(fileno(f_i1), "r");
    seq_i1 = kseq_init(fp_i1);

    /* setup R1 */
    gzFile fp_r1;
    kseq_t *seq_r1;
    FILE *f_r1 = fopen(filename_r1, "rb");
    if (f_r1 == NULL) {
        fprintf(stderr, "Can't open: %s\n", filename_r1);
        exit(EXIT_FAILURE);
    }
    fp_r1 = gzdopen(fileno(f_r1), "r");
    seq_r1 = kseq_init(fp_r1);
    
    /* setup R2 */
    gzFile fp_r2;
    kseq_t *seq_r2;
    FILE *f_r2 = fopen(filename_r2, "rb");
    if (f_r2 == NULL) {
        fprintf(stderr, "Can't open: %s\n", filename_r2);
        exit(EXIT_FAILURE);
    }
    fp_r2 = gzdopen(fileno(f_r2), "r");
    seq_r2 = kseq_init(fp_r2);

    FILE *f_out = stdout;
    if (filename_output_original) {
        if (chunk_size > 0) {
            asprintf(&filename_output, "%s_0.fastq", filename_output_original);
        } else {
            asprintf(&filename_output, "%s", filename_output_original);
        }
        f_out = fopen(filename_output, "w");
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

    char cr[17];
    char cy[17];
    char ur[11];
    char uy[11];
    char *bc = NULL;   // I1.seq
    char *qt = NULL;   // I1.qual
    memset(cr, '\0', 17);
    memset(cy, '\0', 17);
    memset(ur, '\0', 11);
    memset(uy, '\0', 11);

    //name, comment, seq, qual; )
    int counter = 0;
    int file_counter = 0;
    while (kseq_read(seq_r1) >= 0) {
        kseq_read(seq_r2);
        kseq_read(seq_i1);

        memcpy(cr, &seq_r1->seq.s[0], 16);
        memcpy(cy, &seq_r1->qual.s[0], 16);
        memcpy(ur, &seq_r1->seq.s[16], 10);
        memcpy(uy, &seq_r1->qual.s[16], 10);

        bc = seq_i1->seq.s;
        qt = seq_i1->qual.s;
        
        fprintf(f_out, "@%s|||BC||||||QT||||||CR|||%s|||CY|||%s|||UR|||%s|||UY|||%s|||BC|||%s|||QT|||%s %s\n", seq_r1->name.s, cr, cy, ur, uy, bc, qt, seq_r2->comment.s);
        fprintf(f_out, "%s\n+\n%s\n", seq_r2->seq.s, seq_r2->qual.s);

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
                fprintf(stderr, "Generated: %s\n", filename_output);
            }

            asprintf(&filename_output, "%s_%d.fastq", filename_output_original, file_counter);
            
            if (debug_mode) {
                fprintf(stderr, "Generating: %s\n", filename_output);
            }

            f_out = fopen(filename_output, "w");
        }
    }

    kseq_destroy(seq_i1);
    fclose(f_i1);
    gzclose(fp_i1);

    kseq_destroy(seq_r1);
    fclose(f_r1);
    gzclose(fp_r1);

    kseq_destroy(seq_r2);
    fclose(f_r2);
    gzclose(fp_r2);

    fclose(f_out);

    if (debug_mode) {
        fprintf(stderr, "Generated: %s\n", filename_output);

        t_temp = clock() - t_start;
        double time_total = ((double)t_temp)/CLOCKS_PER_SEC; // in seconds

        fprintf(stderr, "Done. Reads: %d, Total Time: %f\n", counter, time_total);
    }

    return 0;
}