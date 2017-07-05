import re
import urllib.request
import os.path

# Full reference genome sequences available: https://www.ncbi.nlm.nih.gov/genome/genomes/11681

def extract_srr_ids(html):
    pattern = "<td>(SRR[0-9]+)</td>"

    with open(html,'r+') as f:
        html_string = f.read()

    return set(re.findall(pattern, html_string))

def ncbi_download_sra(
        srr_list,
        destination_dir = "data/sra/",
        ftp_root = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/"):
    for srr in srr_list:
        srr_base = srr[:6]

        filename = srr + ".sra"
        filepath = destination_dir+filename
        url = ftp_root + srr_base + "/" + srr + "/" + filename

        if os.path.exists(filepath) == False:
            print('Downloading '+filename)
            response = urllib.request.urlretrieve(url, filepath)
        else:
            print('Skipping '+filename)

if __name__ == "__main__":
    srr_list = extract_srr_ids("data/srr_list.html")
    ncbi_download_sra(srr_list)