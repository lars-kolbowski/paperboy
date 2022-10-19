import slack
from Bio import Entrez
import datetime
from biomedrxiv_search_function import biomedrxivsearch


# client = slack.WebClient(token=os.environ['SLACK_API_TOKEN'])
slack_token = "***REMOVED***"
client = slack.WebClient(token=slack_token)

# use yesterday to make sure we don't miss new papers
yesterday = datetime.datetime.now() - datetime.timedelta(1)

yesterday_str = f'{yesterday.year}/{yesterday.month}/{yesterday.day}'
Entrez.email = '***REMOVED***'
handle = Entrez.esearch(
    db='pubmed',
    term='((cross-linking OR (crosslinking OR (CLMS OR (XL-MS OR CX-MS)))) AND protein) OR AP-MS '
         'OR interactomics OR "protein-protein interactions"',
    mindate=yesterday_str,
    maxdate=yesterday_str
)
pubMed_results = Entrez.read(handle)

# Biorxiv API only returns published papers
# https://api.biorxiv.org/details/biorxiv/2018-08-21/2018-08-28/0
# import requests
# import json
#
# journal_list = ['Bioinformatics', 'Biochemistry', 'Molecular Biology', ],
#
# url = f'https://api.biorxiv.org/pubs/biorxiv/' \
#       f'{yesterday.year}-{yesterday.month}-16/{yesterday.year}-{yesterday.month}-{yesterday.day}'
# response = requests.get(url)
# if response.status_code == 200:
#     res = json.loads(response.content)
# else:
#     raise Exception(response.status_code)
#
# try:
#     for res in res['collection']:
#         if res['published_journal'] not in journal_list:
#             continue
#
#         # post result
#         link = f"https://doi.org/{res['preprint_doi']}"
#         txt = f"*{res['preprint_title']}*\n" \
#               f"bioRxiv - {link}\n" \
#               f"{res['preprint_authors']}"
#         client.chat_postMessage(channel='#papers', text=txt)
#
# except KeyError:
#     pass

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

    client.chat_postMessage(
        channel='#papers',
        text=f"I found {n_pubmed_results} crosslinking paper(s) published yesterday ({yesterday_str}):"
    )
    # BioRxiv
    if len(bioRxiv_results_df) > 0:
        for i, row in bioRxiv_results_df.iterrows():
            # post result
            txt = f"*{row['title']}*\n" \
                  f"bioRxiv - {row['url']}\n" \
                  f"{row['authors']}"
            client.chat_postMessage(channel='#papers', text=txt)

    # PubMed
    for res_id in pubMed_results["IdList"]:
        handle = Entrez.esummary(db="pubmed", id=res_id)
        record = Entrez.read(handle)

        link = f"https://doi.org/{record[0]['DOI']}"
        txt = f"*{record[0]['Title']}*\n" \
              f"{record[0]['FullJournalName']} - {link}\n"
        txt += ", ".join([a for a in record[0]['AuthorList']])
        client.chat_postMessage(channel='#papers', text=txt)
    handle.close()
else:
    client.chat_postMessage(channel='#papers',
                            text=f"No papers published yesterday ({yesterday_str}):")


n_bioRxiv_results = len(bioRxiv_results_df)

if n_bioRxiv_results > 0:
    client.chat_postMessage(
        channel='#papers',
        text=f"I found {n_bioRxiv_results} preprint(s) published yesterday ({yesterday_str}):"
    )
    # BioRxiv
    if len(bioRxiv_results_df) > 0:
        for i, row in bioRxiv_results_df.iterrows():
            # post result
            txt = f"*{row['title']}*\n" \
                  f"bioRxiv - {row['url']}\n" \
                  f"{row['authors']}"
            client.chat_postMessage(channel='#papers', text=txt)

else:
    client.chat_postMessage(channel='#papers',
                            text=f"No preprints published yesterday ({yesterday_str}):")
