import numpy as np
# Data viz pkgs
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

NC_mapping = {'A':0, 'C':1, 'G':2, 'T':3 }

def join_str(seq, up_to=None):
    '''
    start index at 1 since we want to exclude the description;
    this will join the x amount of rows into one long string
    '''
    seqstr = ''.join(seq[1:]).replace('\n','')
    if up_to:
        assert 1 < int(up_to) < len(seqstr)
        seqstr = seqstr[:up_to]
    return seqstr

def assign_numeric_array(seqstr):
    '''
    this initialises a zero filled numpy 1d array that has the
    same length as the number of nucleic acids in the sequence
    string. it will then map the string to integer numbers from
    the dictionary so that we can do 'faster' runs on a numpy
    array
    '''
    numeric_seq = np.zeros((len(seqstr)),dtype=int)
    for i,nucA in enumerate(seqstr):
        numeric_seq[i] = NC_mapping[nucA]

    return numeric_seq


def generate_comparison_array(numeric_seq, numeric_seq2):
    '''
    initialises a zero-filled 2d array with a shape equal to the
    length of the sequence string. it will iterate thru the array
    and set 0 if the nucleic acids equal each other; otherwise leave as 1
    '''
    compr_array = np.zeros((numeric_seq.shape[0],
                            numeric_seq.shape[0]),
                            dtype=int)
    for i,s1 in enumerate(numeric_seq):
        for j,s2 in enumerate(numeric_seq2):
            compr = 1
            if s1 == s2:
                compr = 0
            compr_array[i][j] = compr

    return compr_array

def plot_dotplot(compr_array, seqstr1, seqstr2, figsize=None):
    if figsize:
        assert len(figsize) == 2
        if isinstance(figsize, tuple) or isinstance(figsize, list):
            xsize = figsize[0]
            ysize = figsize[1]
        else:
            raise TypeError('Figure size argument must be a list or tupe, i.e. (x,y)')
    else:
        xsize = 12
        ysize = 12

    fig, ax = plt.subplots(figsize=(xsize,ysize))
    im = ax.imshow(compr_array)
    # individual nucleic acids
    xs = list(seqstr2)
    ys = list(seqstr1)
    # the second sequence is placed on x-axis by convention
    ax.set_xticks(np.arange(len(xs)))
    ax.set_yticks(np.arange(len(ys)))
    # label the axes with the respective list entries
    ax.set_xticklabels(xs)
    ax.set_yticklabels(ys)
    # Rotate the tick labels and set their alignment.
    ax.set_title("sequence dot plot - length of nucleic acids={}x{}".format(len(xs), len(ys)))
    fig.tight_layout()
    return fig, ax
