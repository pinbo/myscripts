// to compile: gcc -Wall -g vcf2table.c -o vcf2table
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>

// #define MAX_LINE_LEN 1024
#define MAX_FIELDS 100
#define MAX_ALLELES 10

// This will only have effect on Windows with MSVC
#ifdef _MSC_VER
    #define _CRT_SECURE_NO_WARNINGS 1
    #define restrict __restrict
#endif

/*
POSIX getline replacement for non-POSIX systems (like Windows)
Differences:
    - the function returns int64_t instead of ssize_t
    - does not accept NUL characters in the input file
Warnings:
    - the function sets EINVAL, ENOMEM, EOVERFLOW in case of errors. The above are not defined by ISO C17,
    but are supported by other C compilers like MSVC
*/
int64_t my_getline(char **restrict line, size_t *restrict len, FILE *restrict fp) {
    // Check if either line, len or fp are NULL pointers
    if(line == NULL || len == NULL || fp == NULL) {
        errno = EINVAL;
        return -1;
    }
    
    // Use a chunk array of 128 bytes as parameter for fgets
    char chunk[128];

    // Allocate a block of memory for *line if it is NULL or smaller than the chunk array
    if(*line == NULL || *len < sizeof(chunk)) {
        *len = sizeof(chunk);
        if((*line = malloc(*len)) == NULL) {
            errno = ENOMEM;
            return -1;
        }
    }

    // "Empty" the string
    (*line)[0] = '\0';

    while(fgets(chunk, sizeof(chunk), fp) != NULL) {
        // Resize the line buffer if necessary
        size_t len_used = strlen(*line);
        size_t chunk_used = strlen(chunk);

        if(*len - len_used < chunk_used) {
            // Check for overflow
            if(*len > SIZE_MAX / 2) {
                errno = EOVERFLOW;
                return -1;
            } else {
                *len *= 2;
            }
            
            if((*line = realloc(*line, *len)) == NULL) {
                errno = ENOMEM;
                return -1;
            }
        }

        // Copy the chunk to the end of the line buffer
        memcpy(*line + len_used, chunk, chunk_used);
        len_used += chunk_used;
        (*line)[len_used] = '\0';

        // Check if *line contains '\n', if yes, return the *line length
        if((*line)[len_used - 1] == '\n') {
            return len_used;
        }
    }

    return -1;
}

void parse_vcf(FILE *input, FILE *output) {
    // char line[MAX_LINE_LEN];
    char *fields[MAX_FIELDS];
    char *alleles[MAX_ALLELES + 1];  // REF + ALT alleles
    char *line = NULL;
    size_t len = 0;

    // while (fgets(line, sizeof(line), input)) {
    while(my_getline(&line, &len, input) != -1) {
        if (line[0] == '#') {
            if (line[1] != '#') {
                // fprintf(output, "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
                fprintf(output, "CHROM\tPOS\tID\tREF\tALT\tQUAL");

                char *token = strtok(line, "\t\n");
                for (int i = 0; token != NULL; i++) {
                    if (i >= 9) { // Sample columns start from the 10th column
                        fprintf(output, "\t%s", token);
                    }
                    token = strtok(NULL, "\t\n");
                }
                fprintf(output, "\n");
            }
        } else {
            // This is a data line
            int field_count = 0;
            char *token = strtok(line, "\t\n");
            while (token != NULL) {
                fields[field_count++] = token;
                token = strtok(NULL, "\t\n");
            }

            // Print basic VCF fields
            for (int i = 0; i < 9; i++) {
                if (i == 6 || i == 7 || i == 8)
                    continue;
                if (i > 0) {
                    fprintf(output, "\t");
                }
                fprintf(output, "%s", fields[i]);
            }

            // Parse REF and ALT alleles
            alleles[0] = fields[3]; // REF allele
            int allele_count = 1;
            token = strtok(fields[4], ",");
            while (token != NULL) {
                alleles[allele_count++] = token;
                token = strtok(NULL, ",");
            }

            // Process genotype fields
            for (int i = 9; i < field_count; i++) {
                fprintf(output, "\t");

                // Extract the GT field (assume it's the first field in the sample column)
                char *sample = fields[i];
                char *gt = strtok(sample, ":");

                // Convert GT to ATGC codes
                int allele1, allele2;
                if (sscanf(gt, "%d/%d", &allele1, &allele2) == 2 || sscanf(gt, "%d|%d", &allele1, &allele2) == 2) {
                    // fprintf(output, "%s%s", alleles[allele1], alleles[allele2]); // biallele output
                    // single allele output
                    if (allele1 == allele2) fprintf(output, "%s", alleles[allele1]);
                    else fprintf(output, "%s", "H");
                } else {
                    fprintf(output, "N"); // Handle missing or complex genotypes
                }
            }
            fprintf(output, "\n");
        }
    }
    free(line);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage:\n %s input.vcf output.txt\nOR\nzcat input.vcf.gz | %s - output.txt", argv[0], argv[0]);
        return 1;
    }

    // FILE *input = fopen(argv[1], "r");
    FILE *input = stdin;
    if(strcmp(argv[1], "-")) input = fopen(argv[1], "r");

    if (!input) {
        perror("Error opening input file");
        return 1;
    }

    FILE *output = fopen(argv[2], "w");
    if (!output) {
        perror("Error opening output file");
        fclose(input);
        return 1;
    }

    parse_vcf(input, output);

    if (input != stdin) {
        fclose(input);
    }
    fclose(output);

    return 0;
}
