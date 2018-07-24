#!/usr/bin/python3

import sys

if len(sys.argv) != 4:
    print('Usage', sys.argv[0], '<hosts> <results> <stats>')
    exit()

# create dictionary phage-host
with open(sys.argv[1]) as f:
    lines = f.readlines()

phage_to_host = dict()
for line in lines:
    phage = line.split()[0]
    hosts = line.split()[1]

    phage_to_host[phage] = hosts

# edit results and get stats
with open(sys.argv[2]) as f:
    lines = f.readlines()

correct = 0
false_positive = 0
false_negative = 0

results_out = []

for line in lines:
    phage = line.split()[0]
    spec = line.split()[1]
    prediction = int(line.split()[2][1:-1])
    hosts = phage_to_host[phage]

    if spec in hosts:
        reality = 1
    else:
        reality = 0

    if prediction == reality:
        correct += 1
    else:
        if prediction == 1:
            false_positive += 1
        if prediction == 0:
            false_negative += 1

    results_out.append('{}\t{}\t{}\t{}\n'.format(phage, spec, hosts, prediction))

with open(sys.argv[2] + '.ph', 'w') as f:
    f.writelines(results_out)

all_records = correct + false_positive + false_negative

with open(sys.argv[3], 'w') as f:
    f.write('all\tcorrect\tf_pos\tf_neg\tcorrect_p\tf_pos_p\tf_neg_p\n')
    f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(all_records, correct, false_positive, false_negative,
                                                  correct/all_records,
                                                  false_positive/all_records,
                                                  false_negative/all_records))
