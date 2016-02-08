# -*- coding: utf-8 -*-
"""

Max Schommer

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G') 
    'C'
    >>> get_complement('T')
    'A'
    """
    #I added these tests so that the check is comprehensive

    if nucleotide == 'A': 
    	return 'T'
    elif nucleotide == 'C':
    	return 'G'
    elif nucleotide == 'T':
    	return 'A'
    elif nucleotide == 'G':
    	return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    #This is sufficient, because there should not be an empty input into this function: the input should be checked from the function that calls get_reverse_complement.
    length = len(dna)
    
    rev = ''
    new = ''
    for n in range(length-1, -1, -1):  # This goes through the list backwards and then puts each element in a new list called rev in reverse order.
        rev = rev + dna[n]
    for n in range(0,length):           # This takes the complement of the reverse list
        new = new + get_complement(rev[n]) 
    return new                          # Returns reverse complement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    #This is comprehensive because the string is assumed to begin with a start codon, so a check is not needed.
    
    length = len(dna) - len(dna)%3
    x = 0
    for n in range(0,length,3):
        potentialstopcodon = dna[n]+dna[n+1]+dna[n+2]
        if potentialstopcodon == 'TAG' or potentialstopcodon == 'TAA' or potentialstopcodon == 'TGA':
            x = n
            
            break
        elif (length - n) < 4:
            return dna
    return dna[:x]


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("LOLLOLATGCATGAATGTAGATAGLOLATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe('LOL')
    []
    """

    #This makes sure that you do not include junk DNA in your gene. The extra test checks if you look for a start codon, and not just assume that after a stop codon is a start codon.


   
    n = 0
    check = searchforstartcodon(dna) #This makes sure that the ORF starts with a start codon
    if check == None:               # Here we rule out the case where there is no more dna left with an ORF
        length = 0
    else: 
        dna = check                 # If it passes the test, then dna is defined as the remaining dna after the junk dna is parsed through. 

        length = len(dna) - len(dna)%3
    strings = []
    while n < length:               # This loop goes through the gene in groups of three, and calculates the orfs by adding on strings to a list.
        stringleft = searchforstartcodon(dna[n:])

        if stringleft == None:
            break
        else: # Found start codon
            gene = rest_of_ORF(stringleft)
            strings.append(gene)
            n = len(dna)-len(stringleft)+len(gene)
    return strings

def searchforstartcodon(dna):
    """ Searches for the start codon and then returns all dna after and including the start codon
        dna: a DNA sequence
        returns: a dna sequence withouth the junk dna in the beginning
    >>> searchforstartcodon("LOLLOLATGCATGAATGTAGATAGATGTGCCC")
    'ATGCATGAATGTAGATAGATGTGCCC'
    """
    #This makes sure that the junk dna is ignored between potential genes. As you can see, the 'LOLLOL' is removed from the output.




    stringleft = dna
    for i in range(0, len(dna)-len(dna)%3, 3):
        potentialstartcodon = dna[i]+dna[i+1]+dna[i+2]
        #print potentialstartcodon
        if stringleft == 'ATG':
                return None
        if potentialstartcodon == 'ATG':

            stringleft = stringleft[i:]
            return stringleft

        else:
            continue


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("LOLLOLLOLATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    #I modified the doctest to check for the string not beginning with a start codon.
    # TODO: implement this
    orf1 = find_all_ORFs_oneframe(dna)
    orf2 = find_all_ORFs_oneframe(dna[1:])
    orf3 = find_all_ORFs_oneframe(dna[2:])
    allorfs = orf1+orf2+orf3
    return allorfs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    #This is complete because all of the special cases have been checked by the smaller functions.
    bothstrands = find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
    return bothstrands


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    length = []
    bothstrands = find_all_ORFs_both_strands(dna)                       # This sorts the resulting string by shortest to longest, and then returns the last of the string. 
    bothstrands = sorted(bothstrands, key = lambda seq : len(seq))
    end = len(bothstrands)-1

    return bothstrands[end]

def ORFsinOrder(dna, longest):
    """Finds all orfs and returns a list of them in order of shortest to longest.
    """
    length = []                                                         #Does the same thing as the previous function, but returns all strings that are longer than a threshold.
    bothstrands = find_all_ORFs_both_strands(dna)
    bothstrands = sorted(bothstrands, key = lambda seq : len(seq)) 
    end = len(bothstrands)-2
    for i in range(end, 0 , -1):

        print len(bothstrands[i])
        if len(bothstrands[i]) < longest:

            return bothstrands[i:]
    
        

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    counter = 0
    for i in range(0, num_trials):  
        shuffle = shuffle_string(dna)
        length = len(longest_ORF(shuffle))
        if length >= counter:                             # Creates a counter that dynamically changes to equal the length of the longest string found.
            counter = length
        else:
            continue
    return counter


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """

    proteins = ''
    for i in range(0, len(dna)-len(dna)%3, 3):
        index = 0
        lookmeup = dna[i]+dna[i+1]+dna[i+2]
        n = 0
        proteins = proteins + aa_table[lookmeup]
    return proteins


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    len(dna)
    num_trials = 1500                                                                   
    allorfs = find_all_ORFs_both_strands(dna)                                                 
    longest = len(longest_ORF(dna))
    amirandom = longest_ORF_noncoding(dna, num_trials)
    finallist =  ORFsinOrder(dna , amirandom)
    translated = []
    for i in range(0, len(finallist)-1):                                                #This loop translates all of the resulting potential genes into their final format. 

        translated.append(coding_strand_to_AA(finallist[i]))
    return translated





if __name__ == "__main__": 
    import doctest
    doctest.testmod()

mygenes = gene_finder(load_seq("./data/X73525.fa"))                         #mygenes is the final variable that contains all the potential genes above a threshold. 