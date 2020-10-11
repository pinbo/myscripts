// Writen by Junli Zhang, zhjl86@gmail.com, 2020-09-24
// extract SNPs and indels meet a certain criteria from the Tilling GATK vcf output
// build the app: go build Tilling_indel_snp_filter.go
// get help: ./Tilling_indel_snp_filter -h
// usage example: 
// zcat mapped_PC69PC70.vcf.gz | ./Tilling_indel_snp_filter -i - -minhet 5 -minhom 3 -o filterd_variants_het5hom3_PC69PC70.tsv


package main
import (
    "bufio"
    "fmt"
    "os"
    s "strings"
	"strconv"
	"flag"
)
var version = "undefined"
func main () {	
	fmt.Println("version " + version)
	input := flag.String("i", "example.vcf", "input vcf file name")
	output := flag.String("o", "filterd_variants.tsv", "output file name")
	minhom := flag.Int("minhom", 2, "Minimum coverage for consideration of homozygous variants")
	minhet := flag.Int("minhet", 2, "Minimum coverage for consideration of heterozygous variants")
	minhetper := flag.Float64("minhetper", 0.2, "Minimum percentage for consideration of het")
	minlibs := flag.Int("minlibs", 10, "Minimum number of libraries covered to be considered a valid position")
	indelonly := flag.Bool("indelonly", false, "whether to only get indels")
	minmq := flag.Float64("minmq", 20.0, "Minimum mapping quality")
	maxmutlib := flag.Int("maxmutlib", 1, "maximum number of libaries with the mutations at a position")
	moreinfor := flag.Bool("moreinfor", false, "whether to print more information for debugging in the end of each line")
	nophased := flag.Bool("nophased", false, "whether to include phased calls, which are mostly because multiple SNPs in one read")
	addQual := flag.Bool("addQual", false, "whether to print SNP QUAL in the end")
	minWT := flag.Int("minWT", 1, "Minimum coverage for consideration of a homozygous wild type genotype")
	minQual := flag.Float64("minQual", 100.0, "Minimum vcf quality")
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
		qual, _ := strconv.ParseFloat(ll[5], 64) 
		if qual < *minQual {
			continue
		}
		chrom := ll[0]
		pos := ll[1]
		ref0 := ll[3] // reference allele, CS allele
		alt0 := ll[4] // alternative alleles, either Kronos allele or true mutations
		//altList := s.Split(alt0, ",")
		alleleList := s.Split(ref0 + "," + alt0, ",")
		if len(alleleList) > 2 {// more than two alternative alleles
			continue
		}
		infos := getInfo(ll[7])
		if infos["AC"] == infos["AN"] {//all mutant allele
			continue
		}
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
		var counts []int
		counts = append(counts, totalAlt)
		for _, ii := range AC {
			i, _ := strconv.Atoi(ii)
			totalAlt += i
			counts = append(counts, i)
		} 
		counts[0] = AN - totalAlt // now counts has counts for 0, 1, 2 
		maxcount, refKronos := max(counts) // the max allele is refKronos
		ref0 = alleleList[refKronos]
		if AN - maxcount > 4 && AN - maxcount > 2 * *maxmutlib {// too many genotypes with alterantive alleles
			continue
		}
		// refKronos := 0 // Kronos wt allele
		// altKronos := 1
		// if len(altList) == 1 { // only 1 alternative allele
		// 	if totalAlt == AN {// all lines are homozygous alt, so skip
		// 		continue
		// 	} else if totalAlt > AN - totalAlt {// alt allele is more than wt allele, then CS alt allele is Kronos wt allele
		// 		refKronos = 1
		// 		altKronos = 0
		// 		ref0 = ll[4]
		// 		alt0 = ll[3]
		// 	} 
		// } else {// two alternative alleles
		// 	if totalAlt == AN {// no CS ref allele
		// 		if AC[0] > AC[1] {// the big number is for the ref allele
		// 			refKronos = 1
		// 			altKronos = 2
		// 			ref0 = altList[0]
		// 			alt0 = altList[1]
		// 		} else {
		// 			refKronos = 2
		// 			altKronos = 1
		// 			ref0 = altList[1]
		// 			alt0 = altList[0]
		// 		}
		// 	} else { // there is 0, so 3 alleles
		// 		//continue
		// 		nref := AN - totalAlt
		// 		if nref < totalAlt { // CS allele is the WT allele
		// 			refKronos = 0

		// 		} else {

		// 		}
		// 	}
		// }
		
		DP := infos["DP"] // total depth
		//Geno, refDP, altDP, totalDP, nlib, nmutlib, _, _, firstMutPos := parse3(ll[9:], refKronos, *minhetper)
		wholeGeno, refDP, altDP, totalDP, nlib, nmutlib, _, _, mutPos, snphomhet, altKronosList := parse3(ll[9:], refKronos, *minhom, *minhet, *minhetper, *nophased, *minWT)
		if nmutlib > 0 && nmutlib <= *maxmutlib && nlib >= *minlibs {// mutation only in one lib and at least minlibs have coverage
			for _, mut := range mutPos {
				// homhet := "hom"
				// if s.Contains(Geno[mut], strconv.Itoa(refKronos)) {
				// 	homhet = "het"
				// }
				homhet := snphomhet[mut]
				altKronos := altKronosList[mut]
				alt0 = alleleList[altKronos]
				tt := "++" // variant type: ++ is indel; else AG, TC etc
				inserttype := "."
				// if s.Contains(ll[4], "*"){
				if ref0 == "*" {
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
				nref := refDP[mut]
				nalt := altDP[mut]
				ntotal := totalDP[mut]
				//if (homhet == "hom" && ntotal >= *minhom) || (homhet == "het" && ntotal >= *minhet){
				//outline := s.Join([]string{chrom, pos, ref0, alt0, MQ, DP, libNames[firstMutPos], homhet, strconv.Itoa(nref), strconv.Itoa(nalt), strconv.Itoa(nlib)}, "\t")
				outline := s.Join([]string{chrom, pos, ll[3], DP, ref0, alt0, libNames[mut], homhet, strconv.Itoa(nref), strconv.Itoa(nalt), tt, strconv.Itoa(ntotal), strconv.Itoa(nlib), inserttype}, "\t")
				if *moreinfor {
					outline += "\t" + s.Join([]string{wholeGeno[mut], ll[5], ll[7]}, "\t")
				}
				if *addQual {//add QUAL in the end
					outline += "\t" + ll[5]
				}
				w.WriteString(outline + "\n")
				//}
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
func parse3 (GTs []string, refKronos int, minhom int, minhet int, minhetper float64, nophased bool, minWT int) ([]string, []int, []int, []int, int, int, int, int, []int, []string, []int) {
	//ll := s.Split(line, "\t")
	//GTs := ll[9:]
	size := len(GTs) // number of libs
	//Geno := make([]string, size) // 0/1
	wholeGeno := make([]string, size) // the whole information for the genotype
	refDP := make([]int, size)
	altDP := make([]int, size)
	totalDP := make([]int, size)
	snphomhet := make([]string, size)
	altKronosList := make([]int, size)
	nlib := 0 // number of libs with coverage
	nmutlib := 0 // number of libs with mutations (0), including 0/1
	nwtlib := 0 // number of libs with wt allele (1), including 0/1
	nmissing := 0
	//firstMutPos := 0
	var mutPos []int
	for pos, ii := range GTs {
		wholeGeno[pos] = ii
		ff := s.Split(ii, ":")
		//Geno[pos] = ff[0] // first field is always GT
		//dp, _ := strconv.Atoi(ff[2]) // total depth
		//totalDP[pos] = dp
		ad := s.Split(ff[1], ",")
		nref, _ :=  strconv.Atoi(ad[refKronos])
		//nalt, _ :=  strconv.Atoi(ad[altKronos])
		// nalt := dp - nref
		// refDP[pos] = nref
		// altDP[pos] = nalt
		gt := ff[0]
		if nophased && s.Contains(gt, "|") {
			continue
		}
		wt := strconv.Itoa(refKronos)
		homhet := "hom"
		if s.Contains(gt, "."){
			nmissing += 1
		} else {
			// nlib += 1
			if (gt == wt + "/" + wt || gt == wt + "|" + wt) && nref >= minWT {//homo WT
				nwtlib += 1
			} else {// with mutation
				if s.Contains(gt, wt) {
					homhet = "het"
				}
				// ad := s.Split(ff[1], ",")
				// nref, _ :=  strconv.Atoi(ad[refKronos])
				altKronos, _ := strconv.Atoi(gt[2:3])
				if altKronos == refKronos {
					altKronos, _ = strconv.Atoi(gt[0:1])
				}
				altKronosList[pos] = altKronos
				nalt, _ :=  strconv.Atoi(ad[altKronos])
				refDP[pos] = nref
				altDP[pos] = nalt
				dp := nref + nalt
				totalDP[pos] = dp
				if (homhet == "hom" && nalt >= minhom) || (homhet == "het" && nalt >= minhet && float64(nalt) >= float64(dp) * minhetper){
					nmutlib += 1
					// get mutation position
					mutPos = append(mutPos, pos)
					snphomhet[pos] = homhet
				}
			}
		}
	}
	nlib = nwtlib + nmutlib

	return wholeGeno, refDP, altDP, totalDP, nlib, nmutlib, nwtlib, nmissing, mutPos, snphomhet, altKronosList
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

// find the maximum integars and its positions
func max(v []int)(int, int){
	m := 0
	p := 0
	for i, e := range v {
		if i==0 || e > m {
			m = e
			p = i
		}
	}
	return m, p
}