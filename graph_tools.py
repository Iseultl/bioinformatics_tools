import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def heatmap(matrix, basename, output, genes):
    fig_1, axs = plt.rcParams["figure.figsize"]=16,8
    sns.set_style('whitegrid')
    
    map = sns.heatmap(matrix, annot=False, fmt="g", vmin = -25, vmax = 0, yticklabels = True)
    map.set_yticks(range(1, len(matrix), 4))
    map.set_yticklabels(genes)
    map.set_xticks(range(0, 91, 5))
    map.set_xticklabels(range(-15, 76, 5))
    plt.title('Matrix for '+ basename)
    plt.xlabel('Location Around Start Site')
    #plt.ylabel('Window Size')

    plt.savefig(output + '/' + basename + '_matrix.pdf')