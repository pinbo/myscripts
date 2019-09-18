package main
import (
    "bufio"
    "fmt"
	"os"
	"log"
)
import s "strings"
import "strconv"

//var p = fmt.Println

func gt2snp (allele_list []string , gt string) string {
	c := "" // snp calls
	if s.Contains(gt, ".") {// missing "./." or "."
		c = "N" 
	} else {
		ll := s.Split(gt, "/")
		a := ll[0]
		b := ll[1]
		if a == b {// homozygous
			i, _ := strconv.Atoi(a)
			c = allele_list[i] 
		} else {
			c = "H"
		} // heterozygous 
	}
	return c
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
	for scanner.Scan() {
		//fmt.Println(scanner.Text())
		line := scanner.Text() // s.TrimSuffix(scanner.Text(), "\n")
		if s.HasPrefix(line, "#CHROM") {
			ll := s.Split(line, "\t")
			outline := s.Join(append(ll[0:5], ll[(geno_starts - 1):]...), "\t")
			fmt.Fprintln(outfile, outline)
		}
		if line != "" && !s.HasPrefix(line, "#") {
			ll := s.Split(line, "\t")
			GTs := ll[(geno_starts - 1):] // all GT calls
			ref := []string{ll[3]}
			alt := s.Split(ll[4], ",") // there may be more than one alternative alleles
			alleles := append(ref,  alt...) // this way, 0 will be ref, and 1, 2 ... will be alternative allele
			var SNPs []string
			for _, num := range GTs {
				SNPs = append(SNPs, gt2snp(alleles, num))
			}
			outline := s.Join(append(ll[0:5], SNPs...), "\t")
			//n4, err := w.WriteString("buffered\n")
			fmt.Fprintln(outfile, outline)
		}
	}
}