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
    host = line.split()[1].strip()

    if phage not in phage_to_host:
        phage_to_host[phage] = [host]
    else:
        phage_to_host[phage].append(host)
        phage_to_host[phage].sort()

# edit results and get stats
with open(sys.argv[2]) as f:
    lines = f.readlines()

tp = 0
tn = 0
fp = 0
fn = 0

results_out = []
spec = lines[0].split()[1]

for line in lines:
    phage = line.split()[0]
    prediction = int(line.split()[2][1:-1])
    hosts = phage_to_host[phage]

    if spec in hosts:
        reality = 1
    else:
        reality = 0

    if prediction == 1 and reality == 1:
        tp += 1
    elif prediction == 1 and reality == 0:
        fp += 1
    elif prediction == 0 and reality == 1:
        fn += 1
    elif prediction == 0 and reality == 0:
        tn += 1
    else:
        print('Warning: values out of possible options')

    results_out.append('{}\t{}\t{}\t{}\n'.format(phage, spec, hosts, prediction))

with open(sys.argv[2] + '.ph', 'w') as f:
    f.writelines(results_out)

all_records = tp + tn + fp + fn

accuracy = (tp + tn) / (tp + tn + fp + fn)
sensitivity = tp / (tp + fn)
specificity = tn / (tn + fp)
informedness = sensitivity + specificity - 1    # https://en.wikipedia.org/wiki/Confusion_matrix

with open(sys.argv[3], 'w') as f:
    f.write('model\tall\tt_pos\tt_neg\tf_pos\tf_neg\taccuracy\tsensitivity\tspecificity\tinformedness\n')
    f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(spec, all_records,
                                                  tp, tn, fp, fn,
                                                  accuracy, sensitivity, specificity, informedness))
