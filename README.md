# McArdle_Knowledge_Graph

This repository contains code to construct a Knowledge Graph from Literature filtered for keyword of choice. I have chosen the disease McArdle, to keep the number of papers found at PubMed manageable for my computational resources.

Parts of the module:
- PubMed downloader: retrieve all title and abstracts that match the keyword as txt files, retrieve interesting metadata (year, authors)
- text mining: annotate the NER
