import urllib, urllib2
import xml.etree.ElementTree as ET

def helloworld(output):
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

def getSummary(webEnv, queryKey):
	url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
	url += "efetch.fcgi?db=pubmed&"
	values = dict(db="pubmed", query_key=queryKey, WebEnv=webEnv, rettype="abstract", retmode="text")

	data = urllib.urlencode(values)
	url += data


	print(url)
	req = urllib2.Request(url)
	rsp = urllib2.urlopen(req)

	return(rsp.read())
