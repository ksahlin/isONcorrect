import random 
import os
import errno
from collections import defaultdict

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def custom_histogram(data, outfolder, name='histogram.png', x='x-axis', y='y-axis', title=None, params = {"bins" : 100}):

    plt.hist(data,**params)
    plt.xlabel(x)
    plt.ylabel(y)
    if title:
        plt.title(title)
    if "label" in params:
        plt.legend()
        print("here")
    plt.savefig(os.path.join(outfolder, name))
    plt.clf()

def iteritems(d):
    'Factor-out Py2-to-3 differences in dictionary item iterator methods'
    try:
         return d.iteritems()
    except AttributeError:
         return d.items()

def mkdir_p(path):
    print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def read_config(config_file):
    config = {}
    for line in open(config_file, "r").readlines():
        fields, values = line.strip().split(":")
        try:
            values = [int(values)]
        except:
            values = values.split(" ")
        fields = fields.split(",")
        print(values)
        if len(fields) == len(values):
            for i, field in enumerate(fields):
                try:
                    config[field] = int(values[i])
                except ValueError:
                    config[field] = float(values[i])
        else:
            config[fields[0]] = [int(value) for value in values]

    print(config)
    return config

def transpose(dct):
    d = defaultdict(dict)
    for key1, inner in dct.items():
        for key2, value in inner.items():
            d[key2][key1] = value
    return d


def read_fasta(fasta_file):
    fasta_seqs = {}
    k = 0
    temp = ''
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            fasta_seqs[accession] = ''
            k += 1
        elif line[0] == '>':
            yield accession, temp
            temp = ''
            accession = line[1:].strip().split()[0]
        else:
            temp += line.strip()
    yield accession, temp


def substitution(nucl):
    return random.choice(list(set(["A", "G", "C", "T"]).difference(nucl)))
    # if nucl == "A":
    #     return random.choice(["G","C","T"])
    # elif nucl == "C":
    #     return random.choice(["A","G","T"])
    # elif nucl == "G":
    #     return random.choice(["A","C","T"])
    # elif nucl == "T":
    #     return random.choice(["A","C","G"])
    # else:
    #     return random.choice(["A", "G", "C", "T"])

def insertion(insertion_rate):
    insertion = [random.choice(["A", "G", "C", "T"])]

    # less likely for longer insertions, geometric distribution
    while True:
        if random.uniform(0,1) < insertion_rate:
            insertion.append(random.choice(["A", "G", "C", "T"]))
        else:
            break
    return insertion

def deletion_length(deletion_rate):
    # less likely for longer deletions, geometric distribution
    deletion_size = 1
    while True:
        if random.uniform(0,1) < deletion_rate:
            deletion_size += 1
        else:
            break
    return deletion_size

def weighted_choice(choices):
   total = sum(w for c, w in choices)
   r = random.uniform(0, total)
   upto = 0
   for c, w in choices:
      if upto + w >= r:
         return c
      upto += w
   assert False, "Shouldn't get here"

def mutate_sequence(sequence, mutation_rate, insertion_rate, deletion_rate):
    new_sequence = []
    error_log = []
    total_mut_length = 0
    total_ins_length = 0
    total_del_length = 0
    choices = [("m", mutation_rate), ("i", insertion_rate), ("d", deletion_rate)]
    for i, nucl in enumerate(sequence):
        if random.uniform(0,1) < (1 - (mutation_rate + insertion_rate + deletion_rate)):
            new_sequence.append(nucl)
        else:
            choice = weighted_choice(choices)
            if choice == "m":
                mut_base = substitution(nucl)
                new_sequence.append(mut_base)
                total_mut_length += 1

            elif choice == "i":
                new_sequence.append(nucl)
                insert = insertion(insertion_rate)
                new_sequence.append(insert)
                total_ins_length += len(insert)

            elif choice == "d":
                total_del_length += 1       
                continue
    total_error_length = total_del_length + total_ins_length + total_mut_length

        # if random.uniform(0,1) < mutation_rate:
        #     mut_base = substitution(nucl)
        #     new_sequence.append(mut_base)
        #     error_log.append(("{0}_{1}_{2}|".format(i, nucl, "".join(mut_base))))
        #     total_error_length += 1

        # elif random.uniform(0,1) < deletion_rate:
        #     del_length = deletion_length(deletion_rate)
        #     del new_sequence[-del_length:]
        #     new_sequence.append(nucl)
        #     error_log.append(("{0}_{1}|".format(i, del_length)))
        #     total_error_length += del_length
        #     total_indel_length += del_length
        #     total_del_length += del_length


        # elif random.uniform(0,1) < insertion_rate:
        #     insert = insertion(insertion_rate)
        #     new_sequence.append(insert)
        #     new_sequence.append(nucl)
        #     error_log.append(("{0}_{1}|".format(i, "".join(insert))))
        #     total_error_length += len(insert)
        #     total_indel_length += len(insert)

        # else:
        #     new_sequence.append(nucl)

    new_sequence_flattened = [item for sublist in new_sequence for item in sublist]
    return "".join(nucl for nucl in new_sequence_flattened), error_log, total_error_length, total_ins_length, total_del_length
