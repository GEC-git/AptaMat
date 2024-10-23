import requests
from bs4 import BeautifulSoup
import re

"""
Example of URL
'https://mirbase.org/hairpin/MI0000183'
"""

def single_page_scraper(page):
    scrape=BeautifulSoup(page.text, 'html.parser')
    
    if scrape.find(id="hairpinSequence") is not None:
        seq_db=scrape.find(id="hairpinSequence").find_all('span',limit=2)
    else:
        return None,None,None
    
    for i,elt in enumerate(seq_db):
        if i == 0:
            sequence=elt.string
        if i == 1:
            dotbracket=elt.string
    if scrape.find(href=re.compile('/browse/results/')) is not None:
        family=scrape.find(href=re.compile('/browse/results/')).find('i').string
    else:
        return None,None,None
    
    #print(page.url,family,sequence,dotbracket)
    return family, sequence, dotbracket

def site_scraper(struct_nb):
    tbw=""
    for i in range(1,struct_nb):
        stri=str(i)
        dif=7-len(stri)
        rep=dif*"0"+stri
        URL='https://mirbase.org/hairpin/MI'+rep
        page = requests.get(URL)
        if page.status_code==200:
            family, sequence, dotbracket = single_page_scraper(page)
            if (family,sequence, dotbracket) != (None,None,None):
                tbw+=">"+family+"_"+str(i)+"\n"+sequence+"\n"+dotbracket+"\n"
        print("\r"+str(round((i/struct_nb)*100,3))+"%",end="\r")
    return tbw

def fasta_generator(nb):
    tbw=site_scraper(nb)
    f_created=open("data_clustering_scraped_"+str(nb)+".dat",'a')
    f_created.write(tbw)
    f_created.close()
