import re
from operator import itemgetter
from functools import reduce
import gzip

def reverse_complement(string):
    comp = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'M': 'K',
        'K': 'M',
        'Y': 'R',
        'R': 'Y',
        'S': 'S',
        'W': 'W',
        'N': 'N'
    }
    return ''.join(list(map(lambda c: comp[c], string))[::-1])

def parse_rest(rest, file_format='.gtf'):
    # Parses the 'Attributes' field of the gff
    if file_format == '.gff3':
        return dict(map(lambda x: x.rstrip('\n').split('='), rest.split(';')))
    elif file_format == '.gtf':
        # Sorry, this is a disgusting file format
        # Who the heck puts semicolons in a semicolon-delimited file?!
        PATTERN = re.compile(r'''((?:[^;"']|"[^"]*"|'[^']*')+)''')
        rest = PATTERN.split(rest.rstrip('\n'))[1::2]
        try:
            return dict(map(lambda x: x.rstrip('\n')\
                                       .lstrip(' ')\
                                       .replace('"', '')\
                                       .split(' ', 1), rest))
                            # rest.rstrip(';').split(';')[:-1]))
        except Exception as e:
            print([r.rstrip('\n').lstrip(' ').replace('"', '').split(' ', 1) for r in rest])
            print(rest)


def parse_fasta(path):
    scaffolds = {}
    sequences = []
    current = ''
    with gzip.open(path, 'r') as fh:
        for l in fh:
            line = l.decode('utf-8')
            line = line.rstrip('\n')
            if line.startswith('>'):
                if len(current) > 0:
                    scaffolds[current] = ''.join(sequences)
                sequences = []
                current = line.split()[0][1:]
            else:
                sequences.append(line)

        scaffolds[current] = ''.join(sequences)
    return scaffolds

"""
Per https://stackoverflow.com/a/20062829
"""
def flatten(ivs):
    return reduce(lambda ls, iv: ls + [iv[0], iv[1]], ivs, [])

def unflatten(endp):
    return [(endp[i], endp[i + 1]) for i in range(0, len(endp) - 1, 2)]

def merge(a_ivs, b_ivs, op):
    a_endpoints = flatten(a_ivs)
    b_endpoints = flatten(b_ivs)
    if len(a_ivs) == 0 or len(b_ivs) == 0:
        return unflatten(a_ivs)

    sentinel = max(a_endpoints[-1], b_endpoints[-1]) + 1
    a_endpoints += [sentinel]
    b_endpoints += [sentinel]

    a_index = 0
    b_index = 0

    res = []

    scan = min(a_endpoints[0], b_endpoints[0])
    while scan < sentinel:
        in_a = not ((scan < a_endpoints[a_index]) ^ (a_index % 2))
        in_b = not ((scan < b_endpoints[b_index]) ^ (b_index % 2))
        in_res = op(in_a, in_b)

        if in_res ^ (len(res) % 2): res += [scan]
        if scan == a_endpoints[a_index]: a_index += 1
        if scan == b_endpoints[b_index]: b_index += 1
        scan = min(a_endpoints[a_index], b_endpoints[b_index])

    return unflatten(res)

def interval_diff(a, b):
    return merge(a, b, lambda in_a, in_b: in_a and not in_b)

def interval_union(a, b):
    return merge(a, b, lambda in_a, in_b: in_a or in_b)

def interval_intersect(a, b):
    return merge(a, b, lambda in_a, in_b: in_a and in_b)

def merge_intervals(intervals):
    sorted_intervals = sorted(intervals, key=itemgetter(0))
    if not sorted_intervals:  # no intervals to merge
        return
    low, high = sorted_intervals[0]
    for iv in sorted_intervals[1:]:
        if iv[0] <= high:
            high = max(high, iv[1])
        else:
            yield low, high
            low, high = iv
    yield low, high

def len_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def overlap(iv1, iv2):
    start = end = 0
    if iv2[0] <= iv1[0] <= iv2[1]:
        start = iv1[0]
    elif iv1[0] <= iv2[0] <= iv1[1]:
        start = iv2[0]
    if iv2[0] <= iv1[1] <= iv2[1]:
        end = iv1[1]
    elif iv1[0] <= iv2[1] <= iv1[1]:
        end = iv2[1]
    return (start, end)


def collapse_N(string, k=31):
    # Collapses arbitrary length regions of N-bases into regions of length k
    curr = False
    mod = []
    counter = 0
    for c in string:
        if c == 'N' and curr and counter >= k:
            continue
        elif c == 'N':
            curr = True
            counter += 1
        else:
            curr = False
            counter = 0
        mod.append(c)
    return ''.join(mod)
