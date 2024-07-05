# vcf2table
Convert a biallelic vcf file to a SNP table

- `vcf2table_mygetline.c` use a portable `get_line` funcion to read a whole line. Use codes from https://github.com/sol-prog/fgets-getline-usage-examples
- `vcf2table.c` use the default `get_line` funcion to read a whole line. It works with gcc in cygwin, so I suppose it works in windows OS now. **Use this one**

## compile
`gcc -Wall -g vcf2table.c -o vcf2table`

## Usage
```sh
./vcf2table wheatcap.vcf out.txt
# or from stdin
cat wheatcap.vcf | ./vcf2table - out.txt
```