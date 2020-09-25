// Writen by Junli Zhang, zhjl86@gmail.com, 2020-09-24
// extract SNPs and indels meet a certain criteria from the Tilling GATK vcf output
// build the app: go build Tilling_indel_snp_filter.go
// get help: ./Tilling_indel_snp_filter -h
// usage example: 
// zcat mapped_PC69PC70.vcf.gz | ./Tilling_indel_snp_filter -i - -minhet 5 -minhom 3 -o filterd_variants_het5hom3_PC69PC70.tsv


package main
import (
    "bufio"
    //"fmt"
    "os"
    s "strings"
	"strconv"
	"flag"
)

func main () {
	input := flag.String("i", "example.vcf", "input vcf file name")
	output := flag.String("o", "filterd_variants.tsv", "output file name")
	minhom := flag.Int("minhom", 2, "Minimum coverage for consideration of homozygous variants")
	minhet := flag.Int("minhet", 2, "Minimum coverage for consideration of heterozygous variants")
	minhetper := flag.Float64("minhetper", 0.2, "Minimum percentage for consideration of het")
	minlibs := flag.Int("minlibs", 10, "Minimum number of libraries covered to be considered a valid position")
	indelonly := flag.Bool("indelonly", false, "whether to only get indels")
	minmq := flag.Float64("minmq", 20.0, "Minimum mapping quality")
	maxmutlib := flag.Int("maxmutlib", 1, "maximum number of libaries with the mutations at a position")
	flag.Parse()

	//info, _ := os.Stdin.Stat()
	var scanner *bufio.Scanner
	if *input == "-" {
		scanner = bufio.NewScanner(os.Stdin)
	} else {
		infile, err := os.Open(*input)
		check(err)
		defer infile.Close()
		scanner = bufio.NewScanner(infile)
	}
	outfile, err := os.Create(*output)
	check(err)
	defer outfile.Close()
	w := bufio.NewWriter(outfile)
	//w.WriteString("Chrom\tPos\tRef\tAlt\tMQ\tTotCov\tLib\tHo/He\tWTCov\tMACov\t#libs\n")
	w.WriteString("Chrom\tPos\tRef\tTotCov\tWT\tMA\tLib\tHo/He\tWTCov\tMACov\tType\tLCov\t#libs\tInsertType\n")
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
		chrom := ll[0]
		pos := ll[1]
		ref0 := ll[3] // reference allele, CS allele
		alt0 := ll[4] // alternative alleles, either Kronos allele or true mutations
		altList := s.Split(alt0, ",")
		if len(altList) > 2 {// more than two alternative alleles
			continue
		}
		infos := getInfo(ll[7])
		AN, _ := strconv.Atoi(infos["AN"]) // Total number of alleles in called genotypes
		if AN < *minlibs * 2 {// too many missing data
			continue
		}
		//MQ := infos["MQ"]
		MQ, _ := strconv.ParseFloat(infos["MQ"], 64) // mapping quality
		if MQ < *minmq {
			continue
		}
		AC := s.Split(infos["AC"], ",") // Allele count in genotypes, for each ALT allele. Seems alternative alleles are ordered from more to less.
		totalAlt := 0 // total number of altenative alleles
		for _, ii := range AC {
			i, _ := strconv.Atoi(ii)
			totalAlt += i
		} 
		refKronos := 0 // Kronos wt allele
		altKronos := 1
		if len(altList) == 1 { // only 1 alternative allele
			if totalAlt == AN {// all lines are homozygous alt, so skip
				continue
			} else if totalAlt > AN - totalAlt {// alt allele is more than wt allele, then CS alt allele is Kronos wt allele
				refKronos = 1
				altKronos = 0
				ref0 = ll[4]
				alt0 = ll[3]
			} 
		} else {// two alternative alleles
			if totalAlt == AN {// no CS ref allele
				if AC[0] > AC[1] {// the big number is for the ref allele
					refKronos = 1
					altKronos = 2
					ref0 = altList[0]
					alt0 = altList[1]
				} else {
					refKronos = 2
					altKronos = 1
					ref0 = altList[1]
					alt0 = altList[0]
				}
			} else { // there is 0, so 3 alleles
				continue
			}
		}
		tt := "++" // variant type: ++ is indel; else AG, TC etc
		inserttype := "."
		if s.Contains(ll[4], "*"){
			tt = "++"
		} else {
			if len(ref0) == len(alt0) {
				if len(ref0) == 1 {
					tt = ref0 + alt0
				} else {
					ref1, alt1 := seqdif(ref0, alt0)
					if len(ref1) == 1 {
						tt = ref1 + alt1
					} else {
						tt = "++"
					}
				}
			} else {
				tt = "++"
			}
		}
		if tt == "++" {
			inserttype = tt
		}
		if *indelonly && inserttype == "." {
			continue
		} // skip SNPs if need indels only
		DP := infos["DP"] // total depth
		Geno, refDP, altDP, totalDP, nlib, nmutlib, _, _, firstMutPos := parse3(ll[9:], refKronos, altKronos, *minhetper)
		if nmutlib <= *maxmutlib && nlib >= *minlibs {// mutation only in one lib and at least minlibs have coverage
			homhet := "hom"
			if s.Contains(Geno[firstMutPos], strconv.Itoa(refKronos)) {
				homhet = "het"
			}
			nref := refDP[firstMutPos]
			nalt := altDP[firstMutPos]
			ntotal := totalDP[firstMutPos]
			if (homhet == "hom" && ntotal >= *minhom) || (homhet == "het" && ntotal >= *minhet){
				//outline := s.Join([]string{chrom, pos, ref0, alt0, MQ, DP, libNames[firstMutPos], homhet, strconv.Itoa(nref), strconv.Itoa(nalt), strconv.Itoa(nlib)}, "\t")
				outline := s.Join([]string{chrom, pos, ll[3], DP, ref0, alt0, libNames[firstMutPos], homhet, strconv.Itoa(nref), strconv.Itoa(nalt), tt, strconv.Itoa(ntotal), strconv.Itoa(nlib), inserttype}, "\t")
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
func parse3 (GTs []string, refKronos int, altKronos int, minhetper float64) ([]string, []int, []int, []int, int, int, int, int, int) {
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
		nref, _ :=  strconv.Atoi(ad[refKronos])
		nalt, _ :=  strconv.Atoi(ad[altKronos])
		refDP[pos] = nref
		altDP[pos] = nalt
		gt := ff[0]
		dp, _ := strconv.Atoi(ff[2])
		totalDP[pos] = dp
		if s.Contains(gt, "."){
			nmissing += 1
		} else {
			nlib += 1
			if s.Contains(gt, strconv.Itoa(altKronos)) && float64(nalt) > float64(nref + nalt) * minhetper {// a mutation
				nmutlib += 1
				firstMutPos = pos
			}
			if s.Contains(gt, strconv.Itoa(refKronos)) {// a mutation
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

// find differences of two sequences: from the first difference to the last difference
// seq1 and seq2 are string of the same length
func seqdif (seq1 string, seq2 string) (dif1 string, dif2 string){
	first := 0
	last := 0
	for i, _ := range seq1 {
		if seq1[i] != seq2 [i] {
			if first == 0 {
				first = i
				last = i
			} else {
				last = i
			}
		}
	}
	return seq1[first:(last + 1)], seq2[first:(last + 1)]
}