def gff_summary(gff):
    outfile = open(gff)
    gff = outfile.read()
    outfile.close
    summary = []
    gff_lines = gff.split('\n')

    for line in gff_lines:
        dct = {}
        columns = line.split('\t')
        if ('	gene	' in line):
            
            details = columns[8].split(';')

            for item in details:
                if 'Name=' in item:
                    ID = item.replace('Name=', "")

            start = int(columns[3])
            stop = int(columns[4])
            strand_1  = columns[6]
            dct = {'gene': ID, 'start': start, 'stop': stop, 'strand': strand_1}
            
            summary.append(dct)

    return(summary)