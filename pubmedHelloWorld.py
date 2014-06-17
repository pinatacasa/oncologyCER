import urllib, urllib2
import xml.etree.ElementTree as ET
from Bio import Entrez

def tSF():
    Entrez.email = "franklinzyang@gmail.com"

    # first we query entrez using esearch and store the keys in history
    queryString = "radiation oncology comparison"
    handle = Entrez.esearch(db="pubmed", term=queryString, retmax=1000, usehistory="y")    
    record = Entrez.read(handle)
    keys = (record['WebEnv'], record['QueryKey'])
    
    # next we query entrez using efetch and find all keys of relevant articles
    handle = Entrez.efetch(db="pubmed", query_key=keys[1], WebEnv=keys[0],
		    retstart=0, retmax=100000, rettype="uilist", retmode="text")
    # record = Entrez.read(handle)
    
    # testfetch is a concatenated string of all ids from the esearch query
    val = ''
    countVals = 0
    for line in handle:
	val += line.rstrip('\n') + ','
	countVals += 1 #USED FOR DEBUGGING
    # truncate last value
    val = val[:-1]

# DEBUGGING VALUES RETRIED FROM FIRST EFETCH QUERY
    print val

    handle = Entrez.efetch("pubmed", id=val, retmode="xml")
    
    f = open("fetched_XML_RAW", "w")
    for line in handle:
	f.write(line)
    f.close()

    # return val

def search(output):
    values = dict(term="radation", db="pubmed", usehistory="y")
    baseurl = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    url = baseurl + "esearch.fcgi?"

    data = urllib.urlencode(values)
    url = url + data
    
    req = urllib2.Request(url)
    rsp = urllib2.urlopen(req)
    
    outputXML = rsp.read()

    # writing xml into file so it can be parsed
    f = open(output, 'w')
    f.write(outputXML)
    f.close()

    print(outputXML)


# helloworld conducted a basic search query using Entrez. now we parse xml
def getIdList(outputXMLFile):
	tree = ET.parse(outputXMLFile)
	root = tree.getroot()
	# pull out ID values from XML file and add to list
	idList = [ ]
	for idVal in root.findall('IdList'):
		# maybe we don't need to recurse
		for currId in idVal.findall('Id'):
			# the currId are formatted as strings, so we cast
			idList.append(int(currId.text))	
	
	return idList

# returns (webEnv, queryKey) tuple
def getWebAndKey(outputXMLFile):
	tree = ET.parse(outputXMLFile)
	root = tree.getroot()

	webEnv = root.find('WebEnv').text
	queryKey = root.find('QueryKey').text
	return (webEnv, queryKey)

def fetch(webEnv, queryKey, output):

    #retrieve data in batches of 500
    url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    url += "efetch.fcgi?"
    values = dict(db="pubmed", query_key=queryKey, WebEnv=webEnv, 
		    retstart=0, retmax=1000, rettype="abstract", retmode="text")

    data = urllib.urlencode(values)
    url += data
    req = urllib2.Request(url)
    rsp = urllib2.urlopen(req)

    # writing xml into file so it can be parsed
    f = open(output, 'w')
    f.write(rsp.read())
    f.close()

    return rsp.read()

def getSummary(webEnv, queryKey):
	url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
	url += "efetch.fcgi?"
	values = dict(db="pubmed", query_key=queryKey, WebEnv=webEnv, rettype="abstract", retmode="text")

	data = urllib.urlencode(values)
	url += data


	print(url)
	req = urllib2.Request(url)
	rsp = urllib2.urlopen(req)

	return(rsp.read())
