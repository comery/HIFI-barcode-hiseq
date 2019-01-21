def sortLengthSizeTrim(seqs):
    count = {}
    for s in seqs:
        if s in count.keys():
            count[s] += 1
        else:
            count[s] = 1

    seqs = sorted(seqs, key=lambda k : (len(k), count[k]), reverse=False)
    print(seqs)

seqs = ['ATCGAGAGGTAGTATG',
       'ATCGAGAGGTAGTAT',
       'ATCGAGAGGTAGTAT',
       'ATCGAGAGGTAGTg',
       'ATCGAGAGGTAGTa',
       'ATCGAGAGGTAGTa']

sortLengthSizeTrim(seqs)
