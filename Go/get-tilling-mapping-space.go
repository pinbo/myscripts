// Writen by Junli Zhang, zhjl86@gmail.com, 2021-01-11
// get the mapping space from MAPS mpileup output: parsed_Kronos_mpileup.txt.gz
// just to replace this command: 
// zcat parsed_Kronos_mpileup.txt.gz | awk '{n=0;for(i=7;i<=NF;i+=4) {if (n>4) break; if ($i==".") n++} if (n<5) {ss++; aa[$3]++}}END{print ss; for(i in aa) print i,aa[i]}'
// build the app: go build get-tilling-mapping-space.go
// get help: ./get-tilling-mapping-space -h
// usage example: time zcat parsed_Kronos_mpileup.txt.gz | head -n 1000000 | get-tilling-mapping-space -i -
// Turns out is slower than awk. So just awk commands above.


package main
import (
	"bufio"
	"compress/gzip"
    "fmt"
    "os"
    s "strings"
	"flag"
	// "regexp"
)
var version = "1.0"
func main () {
	fmt.Println("version " + version)
	input := flag.String("i", "parsed_Kronos_mpileup.txt.gz", "MAPS mpileup output file")
	// output := flag.String("o", "mapping-space.txt", "output file name")
	maxmiss := flag.Int("maxmiss", 4, "maximum number of libaries without coverage")
	// minlib := flag.Int("minlib", 20, "minimum number of libaries with coverages")
	flag.Parse()

	//info, _ := os.Stdin.Stat()
	var scanner *bufio.Scanner
	if *input == "-" { // from zcat ouptut
		scanner = bufio.NewScanner(os.Stdin)
	} else {//a gz file
		infile, err := os.Open(*input)
		check(err)
		defer infile.Close()
		rawContents, err := gzip.NewReader(infile)
		check(err)
		defer rawContents.Close()
		scanner = bufio.NewScanner(rawContents)
	}
	// outfile, err := os.Create(*output)
	// check(err)
	// defer outfile.Close()
	// w := bufio.NewWriter(outfile)
	// w.WriteString("Chrom\tPos\n")
	// parse file
	mspace := 0 // mapping space
	ntCount := make(map[string]int)  // a dictionary to counts of A, T, G, C
	for scanner.Scan() {
		//fmt.Println(scanner.Text())
		line := scanner.Text() // "\n" is already trimmed
		ll := s.Split(line, "\t")
		nMiss := 0
		for i := 6; i < len(ll); i+=4 {
			if ll[i] == "." {nMiss += 1}
			if nMiss > *maxmiss {break}
		}
		if nMiss <= *maxmiss {
		// nn := regexp.MustCompile("\\t\\d") // numbers
		// matches := nn.FindAllStringIndex(line, -1)
		// fmt.Println(len(matches))
		// if len(matches) > *minlib {
			mspace += 1
			ntCount[ll[2]] += 1
			// w.WriteString(ll[0]+ "\t" + "ll[1]" + "\n")
		}
	}
	// w.Flush()
	// print on screen A/T/G/C counts
	fmt.Printf("Total mapping space: %d\n", mspace)
	for k, v := range ntCount {
        fmt.Printf("%s -> %d\n", k, v)
	}
}

// check error
func check(e error) {
	if e != nil {
		panic(e)
	}
}
