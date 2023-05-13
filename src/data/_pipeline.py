import os
from pathlib import Path
import pathlib

def bigwig_to_bedgraph_rep(repath):
    """
    Given a directory of data for a replicate, convert all bigwig files to bedgraph files.
    """
    for track_dir in os.listdir(repath):
        track_fs = os.listdir(repath + "/" + track_dir)
        if len(track_fs) > 1:
            for file in track_fs:
                if file.endswith(".bigwig"):
                    os.system("bigWigToBedGraph " + repath + "/" + track_dir + "/" + file + " " + repath + "/" + track_dir + "/" + file.replace(".bigwig", ".bedgraph"))
                    print("converted " + file + " to bedgraph")

def create_genomedata_rep(repath, seqfile):
    """
    Given a directory of
        - sequence: fasta
        - binary tracks: bedGraph
    data for a replicate,
    create a corresponding genomedata archive.
    """