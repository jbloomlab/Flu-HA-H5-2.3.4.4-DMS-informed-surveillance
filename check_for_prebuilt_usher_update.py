"""Check if UShER pre-built tree has been updated.

This script is designed to be run as a cron job to indicate if the
UShER pre-built tree has been updated, which means the repo should
be re-run.

Note it reads the Slack webhook from a file not tracked in this repo.

After creating that webhook file and giving it a valid Slack webhook,
add a line like the following via ``crontab -e``:

    0 5 * * * python3 <full_path>/check_for_prebuilt_usher_update.py

"""


import json
import os
import urllib.request


currdir = os.path.dirname(__file__)


with open(os.path.join(currdir, "results/usher_prebuilt/version_info.txt")) as f:
    existing_version = f.read()

response = urllib.request.urlopen("https://hgdownload.gi.ucsc.edu/hubs/GCF/000/864/105/GCF_000864105.1/UShER_NC_007362.1/fluA.GCF_000864105.1.NC_007362.1.latest.version.txt")
new_version = response.read().decode("utf-8")

if existing_version == new_version:
    msg = "UShER H5N1 version unchanged"
else:
    msg = f"UShER H5N1 version has **changed**"
    
msg += "\n\nnew_version=" + new_version + "\nexisting_version=" + existing_version + "\n"

print(msg)

with open(os.path.join(currdir, "_Slack_webhook_url.txt")) as f:
    webhook_url = f.read()

payload = {"text": msg}
data = json.dumps(payload).encode("utf-8")
req = urllib.request.Request(
    webhook_url,
    data=data,
    headers={'Content-Type': 'application/json'}
)

with urllib.request.urlopen(req) as response:
    if response.status == 200:
        print("Message sent successfully.")
    else:
        raise ValueError(f"Failed to send message. Status code: {response.status}")
