// to compile: gcc -Wall -g vcf2table.c -o vcf2table
// use default getline(); works with gcc in cygwin
// 2024-07-05: use kvec for dynamic vector
// 2024-07-07: split string into known length vector
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void splitsub2(char *str, const char *delim, char **array) // known length
{
    array[0] = str;
    size_t dl = strlen(delim); // delim length
    char *p = str;
    int i = 1;
    while ((p = strstr(p, delim)) != NULL) {
        *p = '\0';
        p += dl;
        array[i++] = p;
    }
}

// #define MAX_LINE_LEN 1024
#define MAX_FIELDS 1000
#define MAX_ALLELES 10

void parse_vcf(FILE *input, FILE *output) {
    char *line = NULL;
    size_t len = 0;
    int nfields = 0; // number of fields
    // char *fields[10];
    char **fields;
    fields = (char **)malloc(sizeof(char *));
    // fields = (char **)realloc(fields, sizeof(char *) * 10);
    char *alleles[MAX_ALLELES + 1];  // REF + ALT alleles

    // while (fgets(line, sizeof(line), input)) {
    while(getline(&line, &len, input) != -1) {
        // if (len == 0) continue; // skip blank lines
        if (line[0] == '#') {
            if (line[1] != '#') {
                // fprintf(output, "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
                fprintf(output, "CHROM\tPOS\tID\tREF\tALT\tQUAL");

                char *token = strtok(line, "\t\n");
                for (int i = 0; token != NULL; i++) {
                    nfields++;
                    if (i >= 9) { // Sample columns start from the 10th column
                        fprintf(output, "\t%s", token);
                    }
                    token = strtok(NULL, "\t\n");
                }
                fprintf(output, "\n");
                printf("nfields is %d\n", nfields);
                fields = (char **)realloc(fields, sizeof(char *) * nfields); // increase the size if necessary
            }
        } else {
            // This is a data line
            // int field_count = 0;
            // char **fields = splitsub(line, "\t", &field_count);
            splitsub2(line, "\t", fields);
            // printf("fields[0] is %s\n", fields[0]);
            // Print basic VCF fields
            fprintf(output, "%s", fields[0]);
            for (int i = 1; i < 6; i++) {
                fprintf(output, "\t%s", fields[i]);
            }

            // Parse REF and ALT alleles
            char *ref = fields[3]; // REF allele
            // int allele_count = 1; // alt allele count
            // char **alleles = splitsub(fields[4], ",", &allele_count);
            splitsub2(fields[4], ",", alleles);
            // Process genotype fields
            // printf("nfields is %d\n", nfields);
            for (int i = 9; i < nfields; i++) {
                // printf("field %d is %s\n", i, fields[i]);
                fprintf(output, "\t");
                // Extract the GT field (assume it's the first field in the sample column)
                char *sample = fields[i];
                char *gt = strtok(sample, ":");
                // Convert GT to ATGC codes
                int allele1, allele2;
                if (sscanf(gt, "%d/%d", &allele1, &allele2) == 2 || sscanf(gt, "%d|%d", &allele1, &allele2) == 2) {
                    // fprintf(output, "%s%s", alleles[allele1], alleles[allele2]); // biallele output
                    // single allele output
                    if (allele1 == allele2) {
                        if (allele1 == 0) fprintf(output, "%s", ref);
                        else fprintf(output, "%s", alleles[allele1 - 1]);
                    } else fprintf(output, "%s", "H");
                } else {
                    fprintf(output, "N"); // Handle missing or complex genotypes
                }
            }
            // free(fields);
            // free(alleles);
            fprintf(output, "\n");
        }
    }
    free(line);
    free(fields);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage:\n %s input.vcf output.txt\nOR\nzcat input.vcf.gz | %s - output.txt\n", argv[0], argv[0]);
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
