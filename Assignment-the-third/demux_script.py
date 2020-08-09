#!/usr/bin/env python

import argparse
import gzip
import itertools
import re

def get_args():
    parser = argparse.ArgumentParser(description='Given the two index files and two read files, demultiplexes them, while also filtering out records with index hopping and low quality indices.')
    parser.add_argument("-1", "--read1", type=str, help='The fastq file which contains read 1')
    parser.add_argument("-2", "--index1", type=str, help='The fastq file which contains index 1')
    parser.add_argument("-3", "--index2", type=str, help='The fastq file which contains index 2')
    parser.add_argument("-4", "--read2", type=str, help='The fastq file which contains read 2')
    parser.add_argument("-i", "--indexfile", type=str, help='The file which contains the indexes, where the indexes must be in the fifth column.')
    return parser.parse_args()

args = get_args()
read1 = args.read1
index1 = args.index1
index2 = args.index2
read2 = args.read2
indexfile = args.indexfile

def get_indexlist(indexfile):
    """Input an a text file with indexes, and return the indexes as elements of a list"""
    indexlist = []
    barcodes = re.compile("[ATCG]{2,}")                                                        #Put we will be wanting to pull out, so we can grab it out soon.
    with open(indexfile,'r') as txt:                                                            #Open the indexes files.
        for line in txt:                                                                        #Go through it, line by line.
            line = line.rstrip()                                                                #Take out the whitespace on the right end.
            break_it_up = line.split('\t')                                                      #Separate it by tabs
            if barcodes.search(line):
                indexlist.append(break_it_up[4])                                                #As long as the indexes are in the fifth column, we pull them out and add them to our list.
    return indexlist

def rev_comp(indexlist):
    """Given an index list, return a dictionary where the key is the index, and the value is the reverse complement of said index"""
    RevDict = {'A':'T', 'G':'C', 'C':'G', 'T':'A', 'N':'N'}
    rev_comp_dict = {}
    indexesfwdback = list(indexlist)
    for i in range(len(indexlist)):
        rev_ele = indexlist[i][::-1]
        rev_com = ''
        for j in rev_ele:
            rev_com+=str(RevDict[j])
        indexesfwdback.append(rev_com)                                  #Make a list contains both the indexes and the reverse complement of the indexes.
        rev_comp_dict[indexlist[i]] = rev_com                           #Generate a dictionary where the key is the index, and the value is the reverse complement of the index.
    return rev_comp_dict, indexesfwdback

def convert_phred(letter):
    """Converts a single character into a phred score"""
    return (ord(letter)-33)

def avg_qscore(phred_sequence):
    """Given a sequence of phred scores, return the average value of the quality scores"""
    Sum = 0
    Length = 0
    for letter in phred_sequence:
        Sum+= convert_phred(letter)
        Length+= 1
    Avg_qual = (Sum/Length)
    return Avg_qual

def possible_perms(indexlist):
    """Given a list of indexes, return a dictionary with permutations of the list elements as keys, and the value as 0 so it may be incremented later."""
    permut_dict = {}
    for permutation in itertools.permutations(indexlist,2):                         #Gives permutation of elements in list, however element1, element1 will have to be added.
        permut_dict[permutation] = 0
    for index in indexlist:
        permut_dict[(index, index)] = 0
    permut_dict['Unk_lowQ'] = 0
    permut_dict['Index_hopped'] = 0
    permut_dict['Correctly_matched'] = 0
    return permut_dict


def open_write_files(indexlist):
    """Given an index list, this function names and opens all of the files to which you shall write."""
    file_list = {}
    for i in range(len(indexlist)):                                              #Here, we are getting the names of our files.
        filename = indexlist[i] + "_Read1" + ".fastq"                #So we have nicely named files.
        filehandle = open(filename,"a")                                    #'a' option is "append", so it creates files if they are not yet created, and opens them for appending otherwise.
        file_list[filename] = filehandle                                   #Genericize this so we can more easily write to it later.
        filename = indexlist[i] + "_Read2" + ".fastq"                #So we have nicely named files for read 2.
        filehandle = open(filename,"a")                                    #'a' option is "append", so it creates files if they are not yet created, and opens them for appending otherwise.
        file_list[filename] = filehandle                                   #Once again, genericize this so we can more easily write to it later.
    #Now, this is great, but we still need to make files for our index hopping, 
        filename = "indexhopped_Read1.fastq"
        filehandle = open(filename,"a")
        file_list[filename] = filehandle
        filename = "indexhopped_Read2.fastq"
        filehandle = open(filename,"a")
        file_list[filename] = filehandle
        filename = "unk_lowqual_Read1.fastq"
        filehandle = open(filename,"a")
        file_list[filename] = filehandle
        filename = "unk_lowqual_Read2.fastq"
        filehandle = open(filename,"a")
        file_list[filename] = filehandle
    return file_list                                                       #Returns the file_list dictionary, which contains the names of all the files as the keys, and their filehandles as the values.

def close_written_files(file_list):
    """The antithesis of the previous function; this function closes all the files which had been opened by the open_write_files function."""
    for files in file_list.values():
        files.close()
    return "Files closed."



def main_demultiplexing(read1,read2,index1,index2,indexfile):
    """Calls the output of other functions, demultiplexes the records, and all. Also tracks of index-hoppings, correctly paired, and low quality/unknown indexes."""
    indexlist = get_indexlist(indexfile)
    rev_comp_dict, indexesfwdback = rev_comp(indexlist)
    permut_dict = possible_perms(indexlist)
    file_list = open_write_files(indexlist)                     #Open all of our files to which we'll write.
    with gzip.open(read1, "rt") as r1f, gzip.open(index1, "rt") as i1f, gzip.open(index2, "rt") as i2f, gzip.open(read2, "rt") as r2f:
        #open main files that contain sequences (still zipped)
        LN = 0
        recordHolder = [[],[],[],[]]
        for r1fLine, i1fLine, i2fLine, r2fLine in zip(r1f, i1f, i2f, r2f):

            LN += 1
            r1fLine = r1fLine.strip()
            i1fLine = i1fLine.strip()
            i2fLine = i2fLine.strip()
            r2fLine = r2fLine.strip()

            recordHolder[0].append(r1fLine)
            recordHolder[1].append(i1fLine)
            recordHolder[2].append(i2fLine)
            recordHolder[3].append(r2fLine)

            if LN % 4 == 0:
                indices2header = str(recordHolder[1][1]) + "-" + str(rev_comp(recordHolder[2][1]))
                recordHolder[0][0] = str(recordHolder[0][0]) + ":" + indices2header
                recordHolder[3][0] = str(recordHolder[3][0]) + ":" + indices2header

                RecordHolder = recordHolder                                             #Somehow this fixes my errors, but if I switch all of them to a single variable the errors pop back up? I'm not sure why...
                i1 = RecordHolder[1][1]
                i2 = RecordHolder[2][1]

                if i1 not in indexlist or i2 not in indexesfwdback or 'N' in i1 or 'N' in i2:
                    permut_dict['Unk_lowQ']+= 1                                            #Update our counter, so we can keep track of the low quality or unknown records
                    for i in range(len(RecordHolder[3])):                                #Iterate through the lines of the record, because as we are no longer manipulating the information, we can now...
                        RecordHolder[3][i] = RecordHolder[3][i] +'\n'                    #...add newlines to the lines so we can write the records to a file
                        RecordHolder[0][i] = RecordHolder[0][i] +'\n'
                    file_list['unk_lowqual_Read1.fastq'].writelines(RecordHolder[0])     #Actually write them to files.
                    file_list['unk_lowqual_Read2.fastq'].writelines(RecordHolder[3])
                elif avg_qscore(RecordHolder[1][3]) < 30 or avg_qscore(RecordHolder[2][3]) < 30:
                    permut_dict['Unk_lowQ']+= 1                                            #Update our counter, so we can keep track of the low quality or unknown records
                    for i in range(len(RecordHolder[3])):                                #Iterate through the lines of the record, because as we are no longer manipulating the information, we can now...
                        RecordHolder[3][i] = RecordHolder[3][i] +'\n'                    #...add newlines to the lines so we can write the records to a file
                        RecordHolder[0][i] = RecordHolder[0][i] +'\n'
                    file_list['unk_lowqual_Read1.fastq'].writelines(RecordHolder[0])     #Actually write them to files.
                    file_list['unk_lowqual_Read2.fastq'].writelines(RecordHolder[3])
                elif rev_comp_dict[i1] == i2:                                             #So, if the second index is the reverse complement of the first index,
                    permut_dict[(i1,i1)]+= 1                                             #Update our counter, so we can keep track of the pairs
                    permut_dict['Correctly_matched']+= 1                                   #And of our overall correct matches.
                    for i in range(len(RecordHolder[3])):                                #Iterate through the lines of the record, because as we are no longer manipulating the information, we can now...
                        RecordHolder[3][i] = RecordHolder[3][i] +'\n'                    #...add newlines to the lines so we can write the records to a file
                        RecordHolder[0][i] = RecordHolder[0][i] +'\n'
                    named_file1 = i1 + '_Read1.fastq'                                       #Then we find the name of the file for the first one,
                    file_list[named_file1].writelines(RecordHolder[0])                   #And write the read1 record to the appropriate file.
                    named_file2 = i1 + '_Read2.fastq'                                       #Similarly, we find the file name for the second one,
                    file_list[named_file2].writelines(RecordHolder[3])                   #And write to that file.
                else:                                                                    #As we have now taken out all of the records where the indexes weren't in our index list, and also the ones where they were in our index and matched properly, then...
                    i2rc = rev_comp(i2)                                                  #Now we get the reverse complement of the 2nd index, so we have what would initially be its forward sequence
                    permut_dict[(i1,i2rc)]+= 1                                           #And increment its specific counter
                    permut_dict['Index_hopped']+= 1                                        #And increment its generic counter
                    for i in range(len(RecordHolder[3])):                                #Iterate through the lines of the record, because as we are no longer manipulating the information, we can now...
                        RecordHolder[3][i] = RecordHolder[3][i] +'\n'                    #...add newlines to the lines so we can write the records to a file
                        RecordHolder[0][i] = RecordHolder[0][i] +'\n'
                    file_list['indexhopped_Read1.fastq'].writelines(RecordHolder[0])     #Then we write the records to the files
                    file_list['indexhopped_Read2.fastq'].writelines(RecordHolder[3])
            RecordHolder = [[],[],[],[]]                                                 #Reset our record holder to zero so that it will be ready for the next set of records given it
    close_written_files(file_list)                                                       #Now we close all our files have opened
    print('Entity'+'\t'+'Times encountered')
    for key,value in permut_dict.items():
        print(key,'\t',value)
    return "Successfully demultiplexed"


main_demultiplexing(read1,read2,index1,index2,indexfile)                                 #Now we actually do it! 
