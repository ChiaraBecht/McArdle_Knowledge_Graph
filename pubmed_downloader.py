from Bio import Entrez

# define E-Mail
Entrez.email = "chiara.becht@web.de"

def count_total_records(keyword):
    """
    Step 1: Retrieve the total number of records for the keyword
    For a provided keyword titles and abstracts within the PubMed database
    is searched and the total number of found articles is returned.

    :params:
        keyword: search term

    :return:
        total_records: number of abstracts where search term was found in either title or abstract
    """
    handle = Entrez.esearch(db = "pubmed",
                            term = f"{keyword}[Title/Abstract]",
                            api_key = "c4507f85c841d7430a209603112dba418607")
    record = Entrez.read(handle)
    total_records = int(record["Count"])
    print("total number of hits for the keyword:", total_records)
    return total_records


def retrieve_pmids(total_records, batch_size, keyword):
    """
    Step 2: obtain all PMIDS in batches and store them in a list.

    :params:
        total_records: number of publications that contain the keyword
        batch_size: number of publications to be accessed from PubMed
                    per acquisition step
        keyword: search term

    :return:
        pmids: list with collected pmids matching the keyword
    """
    pmids = []
    for start in range(0, total_records, batch_size):
        handle = Entrez.esearch(db="pubmed",
                                term=f"{keyword}[Title/Abstract]",
                                retstart=start,
                                api_key = 'c4507f85c841d7430a209603112dba418607',
                                retmax=batch_size)
        record = Entrez.read(handle)
        for item in record["IdList"]:
            pmids.append(item)
    print("length of pmid list:", len(pmids))
    return pmids

# Step 3: download the abstracts
'''abstracts = []
titles = []
for pmid in pmids:
    handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
    record = Entrez.read(handle)
    title = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleTitle"]
    titles.append(title)
    abstract = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
    abstract_text = " ".join(abstract)
    abstracts.append(abstract_text)

print(abstracts)

#def save_abstracts_to_file(pmids, titles, abstracts):
for pmid, title, abstract in zip(pmids, titles, abstracts):
    filename = pmid + '.txt'
    text = title + ' ' + abstract
    with open(filename, "w", encoding="utf-8") as file:
        file.write(text)

print("file downloaded")

def keyword_search_pmids(keyword):
    handle = Entrez.esearch(db="pubmed", term=keyword, api_key = 'c4507f85c841d7430a209603112dba418607', retmax=100)  # Adjust retmax as needed
    record = Entrez.read(handle)
    pmids = record["IdList"]
    print(pmids)

def download_abstract_text():
    pass'''

if __name__ == '__main__':
    keyword = "mcardle"
    batch_size = 100
    total_records = count_total_records(keyword = keyword)
    pmids = retrieve_pmids(total_records, batch_size, keyword)