#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Examples of read and write files, dictionaries, split
#  Brach_search.py
#  "/home/junli/Research/Davis2013/Rust/WSU/Phytozome/PhytozomeV9/Rice/Brach_search2.2.py"
#  Copyright 2014 Junli Zhang <zhjl86@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
### !!!!!!! ###
# Now I use Bradi0009s00210 (gene) as identifier of dictionary rather than
# Bradi0009s00210.1 (transcript) name.
# For exact match, use transcript name.


defDict1 = dict()
defDict2 = dict()
defDict3 = dict()
defDict4 = dict()
pfmDict = dict()
pthrDict = dict()

# PFAM, PANTHER, GO et al numbers
for line in open("Bdistachyon_192_annotation_info.txt","r"):
    line = line.strip('\n')
    identifier = line.split('\t')[1] # Bradi0009s00210
    pfm = line.split('\t')[4] # PF02518,PF00512
    pthr = line.split('\t')[5] # PTHR33333,PTHR44444
    definition = line.split('\t')[12] # Arabidopsis annotation
    defDict1[identifier] = definition
    defDict2[identifier] = pfm
    defDict3[identifier] = pthr

# BLAST2GO
for line in open("Bdistachyon_192_defline.txt","r"):
    line = line.strip('\n')
    identifier1 = line.split('\t')[0] # Bradi0009s00210.1
    identifier2 = identifier1.split('.')[0] # Bradi0009s00210
    definition = line.split('\t')[1] # Description
    defDict4[identifier2] = definition

# Pfam
for line in open("Pfam-A.clans.tsv","r"):
    line = line.strip('\n') # only strip newlines
    identifier = line.split('\t')[0] # PF02518
    definition = line.split('\t')[4] # annotation
    pfmDict[identifier] = definition

# Panther
for line in open("PANTHER9.0_HMM_classifications.txt","r"):
    line = line.strip('\n')
    identifier = line.split('\t')[0] # PF02518
    definition = line
    pthrDict[identifier] = definition

# Outfile
out = open("Brachy.result2.txt", "w")
out.write('Gene\tTranscript\tFunction\n')

for line in open("ch3b.txt","r"):
    identifier = line.strip()
    #print(identifier)
    identifier = identifier.split('\t')[0]
    identifier2 = identifier.split('.')[0]
    if identifier2 in defDict1.keys():
        definition = defDict1[identifier2]
        out.write(identifier + '\t' + identifier2 + '\t' + definition)
        
        # seems I only need to use "if key in Dict"
        # no need to use the Dict.key() method
        if identifier2 in defDict4.keys():
            out.write('\t' + defDict4[identifier2])
        else:
            out.write('\t\t')
        
        pfm = defDict2[identifier2]
        if len(pfm) > 0:
            for pp in pfm.split(','):
                if pp in pfmDict.keys():
                    out.write('\t' + pp + '\t' + pfmDict[pp])
        
        pthr = defDict3[identifier2]
        if len(pthr) > 0:
            for pp in pthr.split(','):
                if pp in pthrDict.keys():
                    out.write('\t' + pthrDict[pp])
        out.write('\n')
    else:
		out.write(identifier + '\t' + identifier2 + '\n')



