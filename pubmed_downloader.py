from Bio import Entrez
import multiprocessing as mp
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# define E-Mail
Entrez.email = "chiara.becht@web.de"

# set data path
cwd = Path(__file__).parent.absolute()
data_dir_p = cwd / 'data'

def data_directory(data_dir_p):
    """
    Create the output directory in the same directory this script is located
    """
    if not(data_dir_p.exists()):
            data_dir_p.mkdir(parents=True, exist_ok=False)

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
        text = title + ' ' + abstract_text
        print(text)
        with open(f'{data_dir_p}/{pmid}.txt', "w", encoding="utf-8") as file:
            file.write(text)
        print("file downloaded")
        print(20*'--')
        return (pmid, "success")#, year)

    except:
        print("no text obtained for {pmid}")
        return (pmid, "fail")#, "-")

def extract_year(pmid):
    """
    parse year from each article

    :params:
        pmid: pubmed id of article to be processed
    
    :return: 
        tuple of pmid and year
    """ 
    handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
    record = Entrez.read(handle)
    if record["PubmedArticle"]:
        print(pmid)
        try:
            year = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][0]["Year"]
            return (pmid, year)
        except:
            try:
                year = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"]
                return (pmid, year)
            except:
                date = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["MedlineDate"]
                year = date[:4]
                return (pmid, year)
    elif record["PubmedBookArticle"]:
        print(pmid)
        year = record["PubmedBookArticle"][0]["BookDocument"]["ContributionDate"]["Year"]
        return (pmid, year)
    else:
        return (pmid, '-')


def download_stats(total_records):
    """
    bar chart that shows count of all articles matching the keyword, the downloaded text files
    and the meta data download.

    :param:
        total_records: number of articles
    """
    texts = pd.read_csv('download_track.csv', names=['pmid', 'indicator'])
    texts
    # Create a Seaborn countplot
    sns.set(style="whitegrid")
    ax = sns.countplot(x=texts["indicator"])
    # Save the countplot to a file
    plt.savefig("data_download_stat.png")

    metadata = pd.read_csv('years.csv', names = ['pmid', 'year'])
    # Create a Seaborn countplot
    sns.set(style="whitegrid")
    ax = sns.countplot(x=metadata["year"])
    # Save the countplot to a file
    plt.savefig("metadata_download_stat.png")

    count_text = (texts["indicator"] == 'success').sum()
    count_metadata = len(metadata.index)- (metadata["year"] == '-').sum()
    labels = ['abstracts', 'metadata', 'records']
    data = [count_text, count_metadata, total_records]
    sns.set(style="whitegrid")
    ax = sns.barplot(x=labels, y=data)
    plt.savefig("download_stat_comp.png")

if __name__ == '__main__':
    keyword = "mcardle"
    batch_size = 100

    # create output directory
    data_directory(data_dir_p)

    # count total records
    total_records = count_total_records(keyword = keyword)

    # retrieve pubmed ids
    pmids = retrieve_pmids(total_records, batch_size, keyword)

    # set up mulitprocessing task to download titles and abstracts concurrently
    cpus = mp.cpu_count()
    if cpus > 4:
        cpus = 4
    with mp.Pool(cpus) as pool:
        results = pool.map(download_title_abstract, pmids)
    print(results)

    with open("download_track.csv", "w") as file:
        for item in results:
            file.write(f"{item[0]},{item[1]}\n")
    
    for pmid in pmids:
        res_tup = extract_year(pmid)
        with open("years.csv", "a") as file:
                file.write(f"{res_tup[0]},{res_tup[1]}\n")
    
    download_stats(total_records)