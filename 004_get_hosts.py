#!/usr/bin/python3

from time import gmtime, strftime, sleep
import argparse
import requests
import json
import bs4
import re
from urllib.parse import quote_plus


def safe_get_page(page):

    ok = False
    unsuccessful = 0
    while not ok and unsuccessful != 3:
        req = requests.get(page)
        if req.status_code != 200:
            print('Error at page: {} {}'.format(page, unsuccessful))
            unsuccessful += 1
        else:
            ok = True

    if unsuccessful == 3:
        return ''
    else:
        tmp1 = re.sub(r'(>)(\s*)(.)', r'\1\3', req.text)
        tmp2 = re.sub(r'(.)(\s*)(<)', r'\1\3', tmp1)
        return tmp2


def safe_get_lineage(table, level):
    try:
        result = table.find_all('a', attrs={"alt": level})[0].get_text().strip()
    except IndexError:
        print('SafeGetError in {}'.format(level))
        result = None
    return result


def parse_page(s):

    genus, species = None, None
    data_list = s.find_all("td", attrs={"valign": "top"})

    if len(data_list) < 1:
        kingdom = None
        return kingdom, genus, species

    data = data_list[0]
    with open('tmp.html', 'w') as f:
        f.write(data.prettify())

    kingdom = safe_get_lineage(data, 'superkingdom')
    genus = safe_get_lineage(data, 'genus')
    species = safe_get_lineage(data, 'species')

    if kingdom != 'Bacteria':
        return kingdom, genus, species

    scientific_name = data.find_all('em', string='Scientific name:')[0].next.next.next
    rank = data.find_all('em', string='Rank:')[0].next.next

    scientific_name = scientific_name.strip('\"')

    if rank == 'genus':
        genus = scientific_name
    if rank == 'species':
        species = scientific_name

    return kingdom, genus, species


def get_taxonomy(host_name, dictionary):

    if host_name in dictionary:
        genus = dictionary[host_name][0]
        species = dictionary[host_name][1]
    else:
        genus = None
        species = None

        default_page = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lvl=0&name={}'
        page = default_page.format(quote_plus(host_name))
        print(page)

        content = safe_get_page(page)
        soup = bs4.BeautifulSoup(content, 'lxml')
        small = soup.find_all('small')

        if len(small) > 0 and 'Did you mean' in str(small[0]):
            print('TYPE: did you mean')
            for i in range(1, len(small)):
                page = default_page.format(quote_plus(small[i].get_text()))
                print(page)
                content = safe_get_page(page)
                result = parse_page(bs4.BeautifulSoup(content, 'lxml'))
                superkingdom, genus, species = result[0], result[1], result[2]
                if superkingdom == 'Bacteria':
                    break
        elif len(small) > 0 and 'for references in articles please use' in str(small[0]):
            print('TYPE: right')
            result = parse_page(soup)
            superkingdom, genus, species = result[0], result[1], result[2]
        else:
            print('TYPE: unexpected')

        dictionary[host_name] = [genus, species]

    return genus, species


parser = argparse.ArgumentParser(description='Unify taxonomy of downloaded phages according to ncbi taxonomy browser tool')
parser.add_argument('-i', '--conversion_file', required=True)
parser.add_argument('-o', '--output_file', required=True)
parser.add_argument('-d', '--dictionary')
args = parser.parse_args()

with open(args.conversion_file) as f:
    lines = f.readlines()

if args.dictionary:
    str2host = json.load(open(args.dictionary))
else:
    str2host = dict()

output = open(args.output_file, 'w')

for i in range(len(lines)):
    phage_id, accession_id, phage_name, host1, host2 = lines[i].rstrip().split('\t')

    print(phage_id, accession_id, phage_name, host1, host2, i)

    for host in [host1, host2]:
        if host == 'NO_DATA':
            genus, species = None, None
            output.write('{}\t{}\t{}\t{}\n'.format(phage_id, genus, species, host))
            output.flush()
        else:
            genus, species = get_taxonomy(host, str2host)
            output.write('{}\t{}\t{}\t{}\n'.format(phage_id, genus, species, host))
            output.flush()

# dict_file = '{}.generated.json'.format(strftime("%y-%m-%d", gmtime()))
# json.dump(str2host, open(dict_file, 'w'), indent=4, sort_keys=True)

output.close()
