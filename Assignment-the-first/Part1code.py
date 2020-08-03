#!/usr/bin/env python


import argparse
from matplotlib import pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description='Make kmer spectrum graph')
    parser.add_argument("-f", "--inputfile", type=str, help='The file used as input.')
    parser.add_argument("-i", "--imagename", type=str, help='The name by which to save the file.')
    parser.add_argument("-r", "--readlength", type=int, help='read sequence length')
    parser.add_argument("-t", "--title", type=str, help='Graph title')
    return parser.parse_args()

parse_args = get_args()
title = parse_args.title
infile = parse_args.inputfile
imgname = parse_args.imagename
readlen = parse_args.readlength


def init_list(array, value=0.0):
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with read-length values of 0.0.'''
    StartList = []
    for i in range(readlen):
        StartList.append(value)
    return StartList

def convert_phred(letter):
    """Converts a single character into a phred score"""
    return (ord(letter)-33)


def populate_list(file):
    """populate_list will call the init_list function to create a list called 'mean_scores'\
    open your FASTQ file and loop through every record, convert the Phred quality score for each letter in a line
    and add it to an ongoing sum in mean_scores for other scores in the same position of the line.
    """
    EmptyList = []
    mean_scores = init_list(EmptyList)
    with open(file, "r") as fh:
        LN = 0
        for line in fh:
            LN+=1
            if LN % 4 == 0:
                char = 0
                for character in line.rstrip():
                    mean_scores[char] += convert_phred(character)
                    char+=1
    return mean_scores, LN


def calculate_mean(file):
    """
    calculate_mean will call the init_list function to create a list called 'mean_scores'\
    open your FASTQ file and loop through every record, convert the Phred quality score for each letter in a line
    and add it to an ongoing sum in mean_scores for other scores in the same position of the line.
    It will then divide that sum by the total number of records, which is the total number of lines, divided by four.
    """
    EmptyList = []
    mean_scores = init_list(EmptyList)
    EmptyList1 = []
    Mean_scores = init_list(EmptyList1)
    with open(file, "r") as fh:
        LN = 0
        for line in fh:
            LN+=1
            if LN % 4 == 0:
                char = 0
                for character in line.rstrip():
                    mean_scores[char] += convert_phred(character)
                    char+=1
        Records = (LN/4) #Counting the number of records in the file.
        Counterrr = 0
        for i in range(len(mean_scores)):
            Mean_scores[i] = (mean_scores[i]/Records)
            Counterrr+=1
    plswork = range(readlen)
    return plswork, Mean_scores  




base_pos, Meannn = calculate_mean(infile)




plt.bar(base_pos, Meannn, align='center',alpha=0.5,color='cyan')
plt.ylabel('Mean Quality Value')
plt.xlabel('Base position')
plt.title(title)


plt.savefig(imgname)



