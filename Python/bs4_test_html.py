#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  bs4.py
#  
#  Copyright 2016 Junli Zhang <junli@wei-ThinkPad-X1-Carbon>
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

## Learn beautiful soup

from bs4 import BeautifulSoup
import urllib2
import re

html_doc = """
<html><head><title>The Dormouse's story</title></head>
<body>
<p class="title"><b>The Dormouse's story</b></p>

<p class="story">Once upon a time there were three little sisters; and their names were
<a href="http://example.com/elsie" class="sister" id="link1">Elsie</a>,
<a href="http://example.com/lacie" class="sister" id="link2">Lacie</a> and
<a href="http://example.com/tillie" class="sister" id="link3">Tillie</a>;
and they lived at the bottom of a well.</p>

<p class="story">...</p>
"""

#soup = BeautifulSoup(html_doc, 'html.parser')
#print soup.name
#soup.find_all("a", limit=2) # for me I am sure there is only one result, so I can set limit=1 or I can just find()

#from datetime import datetime
#print datetime.now()

#url = "https://npgsweb.ars-grin.gov/gringlobal/accessiondetail.aspx?id=1427114"
#content = urllib2.urlopen(url).read()
#print datetime.now()
#soup = BeautifulSoup(content)

#print datetime.now()
#receive = soup.find("th", text = "NPGS received:")
#print receive.string

#print datetime.now()
#print receive.find_next_sibling("td").string
#print datetime.now()

def xstr(s):
    return '' if s is None else str(s).strip()

def get_trait(trait_name, traitlist, valuelist):
	if trait_name in traitlist:
		n = traitlist.index(trait_name)
		value = valuelist[n - 1] # first cell is TH
	else:
		value = "NA"
	return value


def pull_grin(acc):
	url = "https://npgsweb.ars-grin.gov/gringlobal/accessiondetail.aspx?id=" + str(acc)
	content = urllib2.urlopen(url).read()
	soup = BeautifulSoup(content, "lxml")
	pinumber = soup.find_all("h1")[0].string
	region = soup.find("th", text = "Collected from:").find_next_sibling("td").string
	#print region
	receive = soup.find("th", text = "NPGS received:").find_next_sibling("td").string
	#print receive
	pedigree = soup.find("th", text = "Pedigree:").find_next_sibling("td").string
	#print pedigree
	status = soup.find("th", text = "Improvement status:").find_next_sibling("td").string
	#print status Growth Habit
	
	## trait table
	tt = soup.find(id = "ctl00_cphBody_tblCropTrait")
	#print tt.find_all("tr")[1].find_all("th")
	allrows = tt.find_all("tr")
	traits = [item.text for item in allrows[1].find_all("th")]
	values = [item.text for item in allrows[2].find_all("td")]
	habit = get_trait("Growth Habit", traits, values)
	height = get_trait("Plant Height", traits, values)
	## put into a list
	info = [xstr(pinumber), str(acc), xstr(region), xstr(receive), xstr(pedigree), xstr(status), xstr(habit), xstr(height)]
	#print info
	info2 = '\t'.join(info)
	#print info2
	return info2


#print pull_grin(1170349)
#print pull_grin(1427114)

# change the input and output file names here
inputfile = "rust875ID_test.txt"
outputfile = "testout.txt"


accid = []
for line in open(inputfile, "r"):
	line = line.strip()
	if line:
		accid.append(int(line))

# Outfile
out = open(outputfile, "w")
out.write('\t'.join(["Accession", "AccessionID", "Region", "Receive_date", "Pedigree", "Status", "Growth_habit", "Plant_height"]) + "\n")

for id in accid:
	info = pull_grin(id)
	out.write(info + "\n")

print "Finished successfully!\n"






