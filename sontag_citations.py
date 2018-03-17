from collections import defaultdict
import json
import xml.etree.ElementTree as ET

from Bio import Entrez


# usage
# sontag_citations.count_authors(
#     sontag_citations.fetch_author_lists(
#         sontag_citations.search('David Sontag[Author]')
#     )
# )


def count_authors(author_lists):
    author_count = defaultdict(int)
    for author_list in author_lists:
        for author in author_list:
            try:
                author_count['{lastname} ({initials})'.format(
                    lastname=author['LastName'],
                    initials=author['Initials'][0]
                )] += 1
            except KeyError:
                # should never do this
                pass
            
    return sorted(author_count.items(), key=lambda k_v: k_v[1], reverse=True)


def fetch_author_lists(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'franklin.z.yang@gmail.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    # These datastructures are so weird. results has a PubmedArticle attribute which is a list.
    # Iterate over PubmedArticle list
    return [r['MedlineCitation']['Article']['AuthorList'] for r in results['PubmedArticle']
            if 'AuthorList' in r['MedlineCitation']['Article']]

    # returns a list of list of authors, by publication


# return list of ids for articles that match query
def search(query):
    Entrez.email = 'franklin.z.yang@gmail.com'
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='1000',
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    return results['IdList']
