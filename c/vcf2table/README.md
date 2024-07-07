# vcf2table
Convert a biallelic vcf file to a SNP table

- `vcf2table_mygetline.c` use a portable `get_line` funcion to read a whole line. Use codes from https://github.com/sol-prog/fgets-getline-usage-examples
- `vcf2table.c` use the default `get_line` funcion to read a whole line. It works with gcc in cygwin, so I suppose it works in windows OS now. **Use this one**

## compile
```sh
# for local use
gcc -Wall -O3 vcf2table.c -o vcf2table

# for WebAssembly
# need to install emsdk
# https://emscripten.org/docs/getting_started/downloads.html

emcc -Wall  vcf2table.c -o vcf2table.js -s EXPORTED_RUNTIME_METHODS=["callMain"] -s ALLOW_MEMORY_GROWTH=1 -O3

```

## Usage
```sh
./vcf2table wheatcap.vcf out.txt
# or from stdin
cat wheatcap.vcf | ./vcf2table - out.txt
```

## Update
- 2024-07-05: use dynamic vector from kvec; the original version used fixed length vector, which cannot handle vcf files with too many samples.
