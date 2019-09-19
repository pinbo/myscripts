// skip checking comments after the header, but did not increase the speed.
// seems if statement is okay.

package main
import (
    "bufio"
    "fmt"
	"os"
	"log"
)
import s "strings"

//var p = fmt.Println

func gt2snp (allele_list []string , gt string) string {
	// suppose unphased calls, that is separated by "/", not by "|"
	// only GT fields with 4 types: 0/0, 0/1, 1/1, or ./.
	switch gt {
	case "0/0":
		return allele_list[0]
	case "1/1":
		return allele_list[1]
	case "0/1":
		return "H"
	case "./.":
		return "N"
	}
	return gt // in case other cases, for example 2/2 when more than 2 alleles
}

func main () {
	if len(os.Args) < 3 {
		fmt.Println("Please provide 2 arguments: input vcf file, output file name")
		os.Exit(1)
	}
	input := os.Args[1]
	output := os.Args[2]
	geno_starts := 10 // 10th column starts the genotype information

	infile, err := os.Open(input)
	if err != nil {
		log.Fatal(err)
	}

	outfile, err := os.Create(output)
	//w := bufio.NewWriter(outfile)

	defer infile.Close()
	defer outfile.Close()

	scanner := bufio.NewScanner(infile)
	w := bufio.NewWriter(outfile)

	for scanner.Scan() {
		//fmt.Println(scanner.Text())
		line := scanner.Text() // s.TrimSuffix(scanner.Text(), "\n")
		if s.HasPrefix(line, "#CHROM") {
			ll := s.Split(line, "\t")
			outline := s.Join(append(ll[0:5], ll[(geno_starts - 1):]...), "\t")
			//fmt.Fprintln(outfile, outline)
			w.WriteString(outline + "\n")
			break
		}
	}

	for scanner.Scan() {
		line := scanner.Text()
		ll := s.Split(line, "\t")
		GTs := ll[(geno_starts - 1):] // all GT calls
		ref := []string{ll[3]}
		alt := s.Split(ll[4], ",") // there may be more than one alternative alleles
		alleles := append(ref,  alt...) // this way, 0 will be ref, and 1, 2 ... will be alternative allele
		SNPs := make([]string, len(GTs))
		for pos, num := range GTs {
			SNPs[pos] = gt2snp(alleles, num)
		}
		outline := s.Join(append(ll[0:5], SNPs...), "\t")
		//fmt.Fprintln(outfile, outline)
		w.WriteString(outline + "\n")
	}
	w.Flush()
}