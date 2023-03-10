import slack
from Bio import Entrez
import datetime
from credentials import SLACK_API_TOKEN, ENTREZ_EMAIL, SLACK_CHANNEL
from biomedrxiv_search_function import biomedrxivsearch


client = slack.WebClient(token=SLACK_API_TOKEN)

# use yesterday to make sure we don't miss new papers
yesterday = datetime.datetime.now() - datetime.timedelta(1)

yesterday_str = f'{yesterday.year}/{yesterday.month}/{yesterday.day}'
Entrez.email = ENTREZ_EMAIL
handle = Entrez.esearch(
    db='pubmed',
    term='((cross-linking OR (crosslinking OR (CLMS OR (XL-MS OR CX-MS)))) AND protein) OR AP-MS '
         'OR interactomics OR "protein-protein interactions"',
    mindate=yesterday_str,
    maxdate=yesterday_str
)
pubMed_results = Entrez.read(handle)

bioRxiv_results_df = biomedrxivsearch(
    start_date=yesterday - datetime.timedelta(1),
    end_date=yesterday,
    journal='biorxiv',
    subjects=['Bioinformatics', 'Biochemistry', 'Molecular Biology', ],
    kwd=['crosslinking', 'cross-linking', 'CLMS', 'XL-MS', 'CX-MS',
         'AP-MS, interactomics, protein-protein interactions'],
    kwd_type='any',
    max_records=50, max_time=1000,
    cols=['title', 'authors', 'url']
)

n_pubmed_results = len(pubMed_results["IdList"])

if n_pubmed_results > 0:
    txt = f"I found {n_pubmed_results} paper(s) published yesterday ({yesterday_str}):\n\n"
    # PubMed
    for res_id in pubMed_results["IdList"]:
        handle = Entrez.esummary(db="pubmed", id=res_id)
        record = Entrez.read(handle)
        try:
            link = f"https://doi.org/{record[0]['DOI']}"
        except KeyError:
            link = f"https://www.ncbi.nlm.nih.gov/pubmed/{res_id}"
        txt += f"*{record[0]['Title']}*\n" \
               f"{record[0]['FullJournalName']} - {link}\n"
        txt += ", ".join([a for a in record[0]['AuthorList']])
        txt += "\n\n"
    handle.close()
    client.chat_postMessage(channel=SLACK_CHANNEL, text=txt)
else:
    client.chat_postMessage(channel=SLACK_CHANNEL,
                            text=f"No papers published yesterday ({yesterday_str}):")


n_bioRxiv_results = len(bioRxiv_results_df)

if n_bioRxiv_results > 0:
    txt = f"I found {n_bioRxiv_results} preprint(s) published yesterday ({yesterday_str}):\n\n"
    # BioRxiv
    if len(bioRxiv_results_df) > 0:
        for i, row in bioRxiv_results_df.iterrows():
            # post result
            txt += f"*{row['title']}*\n" \
                  f"bioRxiv - {row['url']}\n"
            txt += ", ".join(row['authors'])
            txt += "\n\n"
    client.chat_postMessage(channel=SLACK_CHANNEL, text=txt)
else:
    client.chat_postMessage(channel=SLACK_CHANNEL,
                            text=f"No preprints published yesterday ({yesterday_str}):")
