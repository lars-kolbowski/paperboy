# paperboy

### Setup
Clone the repository

`git clone https://github.com/lars-kolbowski/paperboy.git`

Install pipenv

`cd paperboy && pipenv install`

Create a new App in Slack: https://api.slack.com/authentication/basics, give the App access to your desired channel and copy the `Bot User OAuth Token`

Fill your Slack token, channel and an email address for the Entrez search in the credentials.py.template and save the file as credentials.py

Modify `run_paperboy.sh`, replace `<PATH-TO-REPOSITORY>` with the paperboy path.

Create a cronjob
 `crontab -e`
 
 Add the following line to enable the cronjob (The bot will run every day at 12pm)
 `0 12 * * * <PATH_TO_REPOSITORY>/run_paperboy.sh`
