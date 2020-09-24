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
	//mincov := 1
	minhom := 2
	minhet := 5
	//fmt.Println("Test int/int 2/5")
	//fmt.Println(minhom/minhet)
	minhetper := 0.2
	minlibs := 10 // minimum libs with good coverage

	//info, _ := os.Stdin.Stat()
	var scanner *bufio.Scanner
	if input == "-" {
		scanner = bufio.NewScanner(os.Stdin)
	} else {
		infile, err := os.Open(input)
		check(err)
		defer infile.Close()
		scanner = bufio.NewScanner(infile)
	}
	outfile, err := os.Create(output)
	check(err)
	defer outfile.Close()
	w := bufio.NewWriter(outfile)
	w.WriteString("Chrom\tPos\tRef\tAlt\tMQ\tTotCov\tLib\tHo/He\tWTCov\tMACov\t#libs\n")
	// parse commments
	var libNames []string
	for scanner.Scan() {
		//fmt.Println(scanner.Text())
		line := scanner.Text() // "\n" is already trimmed
		if s.HasPrefix(line, "#CHROM") {
			ll := s.Split(line, "\t")
			libNames = ll[9:]
			break
		}
	}
	// continue for other data
	for scanner.Scan() {
		//fmt.Println(scanner.Text())
		line := scanner.Text() // "\n" is already trimmed

		ll := s.Split(line, "\t")
		//outline := s.Join(append(ll[0:6], ll[9:]...), "\t")
		//w.WriteString(outline + "\n")
		chrom := ll[0]
		pos := ll[1]
		ref0 := ll[3] // reference allele, CS allele
		alt0 := ll[4] // alternative allele, either Kronos allele or true mutations
		if s.Contains(alt0, ",") || len(ref0) == len(alt0){
			continue
		} // skip SNPs and multiple alts
		//ref := ll[3] // Kronos allele
		//alt := ll[4] // mutant
		//format := ll[8]
		infos := getInfo(ll[7])
		MQ := infos["MQ"] // mapping quality
		DP := infos["DP"] // total depth
		Geno, refDP, altDP, totalDP, nlib, nmutlib, _, _, firstMutPos := parse3(ll[9:])
		if nmutlib == 1 && nlib >= minlibs {// mutation only in one lib and at least minlibs have coverage
			// since the reference genome is CS, but mutants are from Kronos, so need to see which is WT allele
			// this is for indels, not easy to do, so here I just consider CS == Kronos alleles
			// if nmutlib > nwtlib {
			// 	ref, alt, nmut, nwt = alt, ref, nwt, nmut // switch ref and alt
			// }
			// Chrom Pos Ref MQ TotCov Lib Ho/He WTCov MACov #libs
			homhet := "hom"
			if s.Contains(Geno[firstMutPos], "0") {
				homhet = "het"
			}
			nref := refDP[firstMutPos]
			nalt := altDP[firstMutPos]
			ntotal := totalDP[firstMutPos]
			if (homhet == "hom" && ntotal >= minhom) || (homhet == "het" && ntotal >= minhet && float64(nalt)/float64(ntotal) >= minhetper){
				outline := s.Join([]string{chrom, pos, ref0, alt0, MQ, DP, libNames[firstMutPos], homhet, strconv.Itoa(nref), strconv.Itoa(nalt), strconv.Itoa(nlib)}, "\t")
				w.WriteString(outline + "\n")
			}
		}
	}
	w.Flush()
}


func check(e error) {
    if e != nil {
        panic(e)
    }
}

// FORMAT is: GT:AD:DP:GQ:PL  1/1:0,2:2:6:90,6,0
// I just need the first 3
// suppose only 1 alternative allele
// only 1 sample has the mutation
func parse3 (GTs []string) ([]string, []int, []int, []int, int, int, int, int, int) {
	//ll := s.Split(line, "\t")
	//GTs := ll[9:]
	size := len(GTs) // number of libs
	Geno := make([]string, size)
	refDP := make([]int, size)
	altDP := make([]int, size)
	totalDP := make([]int, size)
	nlib := 0 // number of libs with coverage
	nmutlib := 0 // number of libs with mutations (0), including 0/1
	nwtlib := 0 // number of libs with wt allele (1), including 0/1
	nmissing := 0
	firstMutPos := 0
	for pos, ii := range GTs {
		ff := s.Split(ii, ":")
		Geno[pos] = ff[0] // first field is always GT
		ad := s.Split(ff[1], ",")
		nref, _ :=  strconv.Atoi(ad[0])
		nalt, _ :=  strconv.Atoi(ad[1])
		refDP[pos] = nref
		altDP[pos] = nalt
		gt := ff[0]
		dp, _ := strconv.Atoi(ff[2])
		totalDP[pos] = dp
		if s.Contains(gt, "."){
			nmissing += 1
		} else {
			nlib += 1
			if s.Contains(gt, "1") {// a mutation
				nmutlib += 1
				firstMutPos = pos
			}
			if s.Contains(gt, "0") {// a mutation
				nwtlib += 1
			}
		}
	}

	return Geno, refDP, altDP, totalDP, nlib, nmutlib, nwtlib, nmissing, firstMutPos
}

// process info
// AC=2;AF=1;AN=2;DP=8;ExcessHet=3.0103;FS=0;MLEAC=9;MLEAF=1;MQ=22.67;QD=32.5;SOR=1.179
func getInfo (line string) map[string]string {
	ll := s.Split(line, ";")
	m := make(map[string]string)
	for _, info := range ll {
		ii := s.Split(info, "=")
		m[ii[0]] = ii[1]
	}
	return m
}
