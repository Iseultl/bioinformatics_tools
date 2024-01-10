from collections import defaultdict

def read_cog(cog_annotation):
    f = open(cog_annotation, 'r')
    data = f.readlines()
    f.close()

    tuple_lst = []

    for line in data[1:]:
        columns = line.split('\t')
        gene = columns[0]
        cog = columns[7]
        tups = (cog, gene)
        tuple_lst.append(tups)
    cog_dct = defaultdict(list)
    count = 0
    for key, val in tuple_lst:
        cog_dct[key].append(val)
        count = count + 1
    return(cog_dct) 