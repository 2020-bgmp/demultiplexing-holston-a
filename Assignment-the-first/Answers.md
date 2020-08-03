# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 |
| 1294_S1_L008_R4_001.fastq.gz | Read 2 |

2. Per-base NT distribution
    1. ![R1](https://raw.githubusercontent.com/2020-bgmp/demultiplexing-holston-a/master/Assignment-the-first/R1_mean_qual_base_pos.png)
    ![R2](https://raw.githubusercontent.com/2020-bgmp/demultiplexing-holston-a/master/Assignment-the-first/R2_mean_qual_base_pos.png)
    ![R3](https://raw.githubusercontent.com/2020-bgmp/demultiplexing-holston-a/master/Assignment-the-first/R3_mean_qual_base_pos.png)
    ![R4](https://raw.githubusercontent.com/2020-bgmp/demultiplexing-holston-a/master/Assignment-the-first/R4_mean_qual_base_pos.png)
    
    2. ```A quality score of 30 would be a good cutoff, because for a quality score of 30 it means it is 99.9% probability of it being correct (one in one thousand chance of an incorrect base call). Per the graphs of the mean quality score per base position, for each read and index each positions averages at least 30, thus if the read is below average quality, it would be cut off.```
    3. ```Command used: ls -1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R[23]_001.fastq.gz | while read FASTQ; do echo $FASTQ; zcat $FASTQ | sed -n "2~4p" | grep "N" | wc -l; done```

```
Output obtained:
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
3976613
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
3328051
Thus, a grand total of 3976613 + 3328051 = 7304664 indexes with undetermined base calls.
```
    
## Part 2
1. Define the problem
    We have four input files (all reads; 2 are indexes and 2 biological) and 24 known indexes. We need to demultiplex the reads so that the biological reads with the corresponding index go in the same respective positions in our forward and reverse FASTQ output files. We also need to filter out the reads where index-hopping or low quality indexes occur and put them into separate FASTQ files.
2. Describe output
        The output will be 52 FASTQ files. Half of these will be forward, and half reverse. So, for each half, there will be 26 files, and 24, related to the 24 indexes, of these files will be where reads that matched in an index-pair will go. Another file will contain reads where index-hopping occurred, and the last file will be where reads with low quality or unknown indexes reside. For both forward and reverse files with correct index pair matching, appended to the header will be the sequence of the indexes. Along with the files being output, there will also be information output which will tell how many read-pairs were properly matched (for each of the 24 indexes), how many times index-hopping occurred, and how many indexes were unknown or had low quality; this will likely be output in a text file.
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [4 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
Start with a shebang, specifying where to find the version of Python to use.

Next, we will import argparse, itertools, and numpy

def get_args():
    """Use argparse to input the files we want."""
    Here we will specify what we are doing with our argparse argument. For this, we will use argparse to specify our input files; four reads, and a file with indexes.
    return parser.parse_args()

def get_indexlist(txtfile):
    """Input an a text file with indexes, and return the indexes as elements of a list"""
    here we will input a textfile, cut it so we only grab the column with the indexes, strip each newline and put each index as an element in a list.
    return indexlist

def open_write_files(indexlist):
    """Given an index list, this function names and opens all of the files to which you shall write."""
    This will open all the files to which we will be writing
    return index1f.fastq index1r.fastq index2f.fastq, etc

def rev_comp(indexlist):
    """Given an index list, return a dictionary where the key is the index, and the value is the reverse complement of said index"""
    Generate a dictionary where they key is the index, and the value is the reverse complement of the index.
    Also return a list where the list contains both the indexes and the reverse complement of the indexes. 
    return rev_comp_dict indexesfwdback_list

def possible_perms(indexesfwdback_list):
    """Given a list of indexes and their reverse complements, return a dictionary with permutations of the list elements as keys, and the value as 0 so it may be incremented later.
    Use itertools to create a dictionary where the keys are all the possible permutations of 2 elements from the list added together, and where the value is zero.
    return permut_dict

def main_function(read1,read2,read3,read4):           this will be renamed to something more descriptive than "main function" at a later date.
    """This is where the real work begins."""
    Open all of the files WITHOUT unzipping them.
        Iterate through them by each read with itertools
            Compare the reverse complement of the index from read2 to the index in read3.
            if they match, append index-index to the headers, and write the forward record to the indexf file, and the reverse to the reverse file.
                increment the value of the index-index key from the permut_dict by 1.
            if they don't match, well, still append index-index, but write them to  indexhopf and indexhopr instead of the above files.
                increment the value of the index-index key from the permut_dict by 1.
            Are there N's in the indexes, so that they don't match and aren't in the dictionary?
                If the above is true, then still append index-index, but write them to  unk_lowqf and unk_lowqr instead of the above files.
    Close all of the files.
    return the files, as well as the number keys and values of the permut_dict


5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement


