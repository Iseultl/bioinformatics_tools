def fasta_dict(fasta):
    f = open(fasta, 'r')
    data = f.readlines()
    f.close()

    Fasta = ""

    for item in data[1:]:
        Fasta = Fasta + item.replace('\n', '') 

    fasta_lst = [*Fasta]
    positions = [*range(1, len(fasta_lst)+1)]

    #Make a dictionary who's keys are nt positions and values are nts

    fasta_dict = {positions[i]: fasta_lst[i] for i in range(len(positions))}

    return(fasta_dict)

def read_fasta(fasta_file):

    #Read the fasta file and store the nucleotide sequence as a string 

    fasta = open(fasta_file, 'r')
    nt_lst = fasta.readlines()[1:]
    fasta_str = ''
    for elem in nt_lst:
        elem = elem.replace('\n', '')
        fasta_str = fasta_str + elem
    
    return(fasta_str)

    