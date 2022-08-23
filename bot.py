import slack
from Bio import Entrez
import datetime

# client = slack.WebClient(token=os.environ['SLACK_API_TOKEN'])
slack_token = "***REMOVED***"
client = slack.WebClient(token=slack_token)

# use yesterday to make sure we don't miss new papers
yesterday = datetime.datetime.now() - datetime.timedelta(1)

yesterday_str = f'{yesterday.year}/{yesterday.month}/{yesterday.day}'
Entrez.email = '***REMOVED***'
handle = Entrez.esearch(
    db='pubmed',
    term='(cross-linking OR (crosslinking OR (CLMS OR (XL-MS OR CX-MS)))) AND protein',
    mindate=yesterday_str,
    maxdate=yesterday_str
)
results = Entrez.read(handle)
n_results = len(results["IdList"])
if n_results > 0:

    client.chat_postMessage(
        channel='#papers',
        text=f"I found {n_results} crosslinking paper(s) published yesterday ({yesterday_str}):"
    )

    for id in results["IdList"]:
        handle = Entrez.esummary(db="pubmed", id=id)
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
