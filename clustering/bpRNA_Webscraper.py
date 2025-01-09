import requests
from bs4 import BeautifulSoup
import re

"""
Family:
https://bprna.cgrb.oregonstate.edu/search.php?query=bpRNA_RFAM_1
Right after :
<b>Reference Name:</b> 

Download URL example :
https://bprna.cgrb.oregonstate.edu/dbnFiles/bpRNA_RFAM_1.dbn

#Name: bpRNA_RFAM_1
#Length:  117 
#PageNumber: 1
CCCGGUGACUAUAGAGAGAGGGCCACACCCGUUCCCAUCCCGAACACGGAAGUUAAGCCUCUCAUCGCUGAUGGUACUAUGUGGUUCGCUGCAUGGGAGAGUAGGACGUUGCCGGGU
.((((((((....((.(((((...((..((((((.......))..))))..))....)))))..))(((.((..(.((....((....)).....)).).)).))).))))))))..
"""


def single_page_scraper(page_family,page_download):
    scrape_fam=BeautifulSoup(page_family.text, 'html.parser')
    family_name=scrape_fam.find(string="Reference Name:").find_next(string=True)
    name=family_name[8:]
    family=family_name[:8]
    scrape_down=BeautifulSoup(page_download.text, 'html.parser')
    sequence_dotbracket=scrape_down.find(string=re.compile("#Name: ")).splitlines(keepends=False)
    for i,elt in enumerate(sequence_dotbracket):
        if i == 3:
            sequence=elt
        elif i == 4:
            dotbracket=elt
    
    return family,sequence,dotbracket
    


def site_scraper(struct_nb):
    "Scrapes the entire website of *struct_nb* structures"
    tbw=""
    for i in range(1,struct_nb+1):
        
        URL_fam='https://bprna.cgrb.oregonstate.edu/search.php?query=bpRNA_RFAM_'+str(i)
        URL_down='https://bprna.cgrb.oregonstate.edu/dbnFiles/bpRNA_RFAM_'+str(i)+'.dbn'
        page_fam = requests.get(URL_fam)
        page_down= requests.get(URL_down)
        if page_fam.status_code==200 and page_down.status_code==200:
            family, sequence, dotbracket = single_page_scraper(page_fam,page_down)
            tbw+=">"+family+"_"+str(i)+"\n"+sequence+"\n"+dotbracket+"\n"
        print("\r"+str(round((i/struct_nb)*100,3))+"%",end="\r")
    return tbw

def fasta_generator(nb):
    "Creates a fasta file with all the scrape data to then be used in clustering"
    tbw=site_scraper(nb)
    f_created=open("data_clustering_BPRNA_scraped_"+str(nb)+".dat",'a')
    f_created.write(tbw)
    f_created.close()
