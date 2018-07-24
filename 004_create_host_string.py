#!/usr/bin/python3

import sys


def get_hosts_string(hosts):
    if len(hosts) == 0:
        return 'no_host'
    else:
        return ';'.join(list(set(hosts)))


if len(sys.argv) != 2:
    print('Usage:', sys.argv[0], '<003_deduplicated.genomes.conversion>')
    exit()

with open(sys.argv[1]) as f:
    records = f.readlines()

records.sort()

phage = records[0].split('\t')[0]
hosts = []

for record in records:
    record_phage = record.split('\t')[0]
    record_host1 = record.split('\t')[3].strip().strip('\'').lower().replace(' ', '_')
    record_host2 = record.split('\t')[4].strip().strip('\'').lower().replace(' ', '_')

    if record_phage != phage:

        print('{}\t{}'.format(phage, get_hosts_string(hosts)))

        phage = record_phage
        hosts = []

    if record_host1 != 'no_data':
        hosts.append(record_host1)
    if record_host2 != 'no_data':
        hosts.append(record_host2)

print('{}\t{}'.format(phage, get_hosts_string(hosts)))
