import urllib, urllib2
import xml.etree.ElementTree as ET
from Bio import Entrez

def searchFetchMultiple():
    # 1) create list of queries, iterate through list and esearch/efetch
    # 2) esearch/efetch generates list of ids, make sure ids are unique
    # 3) pass unique ids in the form of a string into function that generates
    #	xml data for each id

    # STEP ONE
    queryList = [ ]
    baseQuery = "radiation[Title/Abstract/ AND oncology[Title/Abstract]"
    queryList.append(baseQuery + " AND retrospective[Title/Abstract]")
    queryList.append(baseQuery + " AND case control")
    queryList.append(baseQuery + " AND matched pair")
    queryList.append(baseQuery + " AND SEER")
    
    currOutputFilePrefix = "fetchXML_DATA_RAW_"
    currOutputFileCount = 0
    articleIDList = [ ]
    for query in queryList:
	# create filename
	currOutputFile = currOutputFilePrefix + str(currOutputFileCount)

	articleIDList = searchFetch(query, currOutputFile, articleIDList)
	# parseXMLIntoYears(currOutputFile + "OUTPUT")

	currOutputFileCount += 1

    efetchToXML('fetchXML_DATA_COMPILED', articleIDList)
    parseXMLIntoYears('fetchXML_DATA_COMPILED')

# input: queryString (query), outputFile (name of file)
#	articleIDList (list of unique articleIDs)
# returns list of IDS
def searchFetch(queryString, outputFile, articleIDList):
    Entrez.email = "franklinzyang@gmail.com"

    # first we query entrez using esearch and store the keys in history
    # queryString = "radiation oncology[Title/Abstract] AND retrospective[Title/Abstract]"
    handle = Entrez.esearch(db="pubmed", term=queryString, usehistory="y")    
    record = Entrez.read(handle)

    print record

    keys = (record['WebEnv'], record['QueryKey'])
    
    # next we query entrez using efetch and find all keys of relevant articles
    handle = Entrez.efetch(db="pubmed", query_key=keys[1], WebEnv=keys[0],
		    retstart=0, retmax=100000, rettype="uilist", retmode="text")
    # record = Entrez.read(handle)
    
    # testfetch is a concatenated string of all ids from the esearch query
    countVals = 0
    for articleID in handle:
	# skip article if it has already been seen
	if articleID in articleIDList:
	    continue

	articleIDList.append(articleID)

    return articleIDList

def efetchToXML(outputFile, articleIDList):
    
    ids = ''
    # generate list of ids from articleIDList
    for articleID in articleIDList:
	ids += articleID.rstrip('\n') + ','
    # truncate last value
    ids = ids[:-1]

    # TODO: Come up with a better name for handle
    handle = Entrez.efetch("pubmed", id=ids, retmode="xml")
    
    f = open(outputFile, "w")
    for line in handle:
        f.write(line)
    f.close()

# returns list of IDS
def searchFetchDeprecated(queryString, outputFile):
    Entrez.email = "franklinzyang@gmail.com"

    # first we query entrez using esearch and store the keys in history
    # queryString = "radiation oncology[Title/Abstract] AND retrospective[Title/Abstract]"
    handle = Entrez.esearch(db="pubmed", term=queryString, usehistory="y")    
    record = Entrez.read(handle)

    print record

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
    # print val

    handle = Entrez.efetch("pubmed", id=val, retmode="xml")
    
    f = open(outputFile, "w")
    for line in handle:
        f.write(line)
    f.close()

# <MedlineCitation><DateCreated><Year>
def parseXMLIntoYears(xmlFile):
    tree = ET.parse(xmlFile)
    root = tree.getroot()
    # pull out ID values from XML file and add to list
    idList = [ ]
    yearMap = dict()
    articleSet = root.findall('PubmedArticle')

    countArticle = 0
    for articleParent in articleSet:
	article = articleParent.find('MedlineCitation').find('Article')
	articleTitle = article.find('ArticleTitle').text
	countArticle += 1

	yearContainer = article.find('Journal').find('JournalIssue').find('PubDate').find('Year')

	# if yearContainer not found, date was not formatted properly.
	# there's nothing we can do right now... just skip
	if yearContainer is None:
	    continue

	year = int(yearContainer.text)
	
	#check to see if this year already has articles in dict
	if year not in yearMap:
	    yearMap[year] = 1
	else:
	    yearMap[year] += 1
    
    print countArticle

    f = open(xmlFile + "_RESULTS" + ".txt", "w")
    for year in yearMap:
	f.write(str(year) + "\t" + str(yearMap[year]) + "\n")
    f.close()
    return yearMap

    

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
