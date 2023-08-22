from Bio import Entrez
import multiprocessing as mp

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

def download_title_abstract_loop(pmids):
    """
    Step 3: download the titles and abstracts and save as text files.

    :params:
        pmids: list with pmids to download
    """
    abstracts = []
    titles = []
    for pmid in pmids:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        record = Entrez.read(handle)
        title = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleTitle"]
        titles.append(title)
        abstract = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
        abstract_text = " ".join(abstract)
        abstracts.append(abstract_text)

    for pmid, title, abstract in zip(pmids, titles, abstracts):
        filename = pmid + '.txt'
        text = title + ' ' + abstract
        with open(filename, "w", encoding="utf-8") as file:
            file.write(text)

def download_title_abstract(pmid):
    """
    Step 3: download the titles and abstracts and save as text files.

    :params:
        pmids: list with pmids to download
    """
    # obtain text
    try:
        print(f'obtaining text for abstract {pmid}')
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        record = Entrez.read(handle)
        title = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleTitle"]
        abstract = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
        abstract_text = " ".join(abstract)

        # save to file
        print("saving file")
        filename = pmid + '.txt'
        print(filename)
        text = title + ' ' + abstract_text
        print(text)
        with open(filename, "w", encoding="utf-8") as file:
            file.write(text)
        print("file downloaded")
        print(20*'--')
        return filename
    except:
        print("no text obtained for {pmid}")
        filename = pmid + '.txt'
        return filename

    


if __name__ == '__main__':
    keyword = "mcardle"
    batch_size = 100
    total_records = count_total_records(keyword = keyword)
    pmids = retrieve_pmids(total_records, batch_size, keyword)

    # set up mulitprocessing task to download 10 referenced articles concurrently
    cpus = mp.cpu_count()
    if cpus > 4:
        cpus = 4
    with mp.Pool(cpus) as pool:
        results = pool.map(download_title_abstract, pmids)