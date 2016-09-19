# combine multiple files with the same header
awk 'NR == 1 || FNR > 1'  ld_sitebyall*.txt > ld_subset_by_all.txt

# Find the maximum values in 5th column for each distinct values in 1st column using awk
# for me to find the chromosome length in the genetic map:
# $1 is chromosome number, $5 is position
awk '{if (a[$1] < $5) a[$1] = $5} END {for (i in a) print i, a[i]}' infile

# still to find the chromosome length, same as before, but needs a sorted genetic map: 
# p is the chromosome name in the previous line, pl is the chromosome name + position 
# in the previous line; NR>1 is just to let it compare from the 2nd line
# END{print pl}: print the last line after finishing, because no one will compare with it.
awk 'NR > 1 && $1 != p {print pl } { p = $1; pl = $1"\t"$5 } END{print pl}' infile

# Read matches from the first file, and remove lines match them in the 2nd file
# The NEXT statement forces awk to immediately stop processing the current record and go on to the next record.
# This means that no further rules are executed for the current record, and the rest of the current rule’s action isn’t executed.
# so the cmd below means: if NR==FNR (e.g the first file), a[$0], else (the other files), if !($2 in a) then print
awk 'NR==FNR {a[$0];next} !($2 in a)' seq_in_both100.txt dist_vs_1RS_0616_300L_90identity.tsv > dist_vs_identity_no100.tsv

# another example to find unique record: here I need to find information of markers of interest
# xclip -o output my marker list, column 2 of match90k.txt is the marker name
xclip -o | awk 'NR==FNR {a[$0];next} ($2 in a)' - match90K.txt | xclip


# If you are using pdftotext you can use the -layout flag to preserve the layout of the text on the pages in your input pdf file:
# There is also -table for table layouts specifically, works great
pdftotext -layout input.pdf output.txt


# Print the next two (i=2) lines after the line matching regexp:
awk '/regexp/{i=2;next;}{if(i){i--; print;}}' file.txt

# Print the line and the next two (i=2) lines after the line matching regexp:
awk '/regexp/{i=2+1;}{if(i){i--; print;}}' file.txt


# first 5 hit for each blast
# it seems ';' between two individual arguments are not necessary
grep -v "^#" 143940316743.tsv | awk 'p!=$1{n=0}; n<5{p=$1; n++; print}' > first5hits.txt
# same as above, but without ;
grep -v "^#" 143940316743.tsv | awk 'p!=$1{n=0} n<5{p=$1; n++; print}' > first5hits.txt
# OR
grep -v "^#" 143940316743.tsv | awk 'p!=$1{i=5} i>0{i--; p=$1; print}' | less -S
#OR use i as a boolen
grep -v "^#" 143940316743.tsv | awk 'p!=$1{i=5} i{i--; p=$1; print}' | less -S


## Just the first hit for each blast
grep -v "#" snp90_blast_phsudo_molecule.txt | sort -k1,1 -u -m | less -S

# sort by multiple columns
# example below: first sort by the first field (alphabeta), then sort by 12th field reversely and numerically
sort -k1,1 -k12,12rn infile

# simple calculator using perl command line
# -l: add new line; -e: exacuate command
perl -le "print sqrt(100)"

# OR use Rscript, very simple
Rscript -e ' sqrt(100)'

# remove leading line numbers from the history command
history | tail -20 | awk '{$1="";print $0}'
history | tail -20 | awk '{$1=""}1'


# get strings in a quote: "(.*?)"
xclip -o | sed -r 's/Line "(.*?)" not found./\1/g'
xclip -o | perl -pe 's/Line "(.*?)" not found./\1/g'

