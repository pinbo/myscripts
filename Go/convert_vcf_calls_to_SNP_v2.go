// skip checking comments after the header, but did not increase the speed.
// seems if statement is okay.
// 2. use bufio.NewWriter: doubled the speed
// 3. minor editing: not sure whether increase the speed or not.
// 4. edited gt2snp, increased the speed a lot.

package main
import (
    "bufio"
    "fmt"
    "os"
    s "strings"
    "strconv"
)

func main () {
	// my program starts
	if len(os.Args) < 3 {
		fmt.Println("Please provide 2 arguments: input vcf file, output file name")
		os.Exit(1)
	}
	input := os.Args[1]
	output := os.Args[2]
	geno_starts := 10 // 10th column starts the genotype information

	infile, err := os.Open(input)
	check(err)

	outfile, err := os.Create(output)
	check(err)

	defer infile.Close()
	defer outfile.Close()

	scanner := bufio.NewScanner(infile)
	w := bufio.NewWriter(outfile)

	for scanner.Scan() {
		//fmt.Println(scanner.Text())
		line := scanner.Text() // "\n" is already trimmed
		if s.HasPrefix(line, "#CHROM") {
			ll := s.Split(line, "\t")
			outline := s.Join(append(ll[0:5], ll[(geno_starts - 1):]...), "\t")
			w.WriteString(outline + "\n")
			break
		}
	}

	for scanner.Scan() {
		line := scanner.Text()
		ll := s.Split(line, "\t")
		GTs := ll[(geno_starts - 1):] // all GT calls
		ref := ll[3] // reference allele
		alt := ll[4] // alternative allele(s)
		SNPs := make([]string, len(GTs))
		for pos, gt := range GTs {
			SNPs[pos] = gt2snpv2(ref, alt, gt)
		}
		outline := s.Join(append(ll[0:5], SNPs...), "\t")
		w.WriteString(outline + "\n")
	}
	w.Flush()
}

func check(e error) {
    if e != nil {
        panic(e)
    }
}

func gt2snpv2 (ref string , alt string, gt string) string {
	// suppose unphased calls, that is separated by "/", not by "|"
	switch gt {
	case "0/0":
		return ref
	case "1/1":
		return alt
	case "0/1":
		return "H"
	case "./.":
		return "N"
	default: // if more than 2 alleles
		alleles := s.Split(alt, ",")
		ll := s.Split(gt, "/")
		a := ll[0]
		b := ll[1]
		if a == b {// homozygous
			i, _ := strconv.Atoi(a)
			return alleles[i-1]
		} else {
			return "H"
		}
	}
}
