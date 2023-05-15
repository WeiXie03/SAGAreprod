import os
from pathlib import Path
import pathlib

def bigwig_to_bedgraph_dir(dir_path):
    """
    Convert all bigwig files in `dir_path` to bedgraph files.
    """
    for track_dir in os.listdir(dir_path):
        track_fs = os.listdir(dir_path + "/" + track_dir)
        if len(track_fs) > 1:
            for file in track_fs:
                if file.endswith(".bigwig"):
                    os.system("bigWigToBedGraph " + dir_path + "/" + track_dir + "/" + file + " " + dir_path + "/" + track_dir + "/" + file.replace(".bigwig", ".bedgraph"))
                    print("converted " + file + " to bedgraph")

def create_genomedata_ct(repath: str, seqfile: str) -> None:
    """
    Given a directory of
        - alignment: bam
        - binary tracks: bedGraph
    data for a cell type,
    and a sequence file for the genome,
    create a corresponding genomedata archive.
    """
    for track_dir in os.listdir(repath):
        track_fs = os.listdir(repath + "/" + track_dir)
        track_list = ""
        # need both bedgraph and bam
        if len(track_fs) > 1:
            track_list += "-t {}={}".format(track_dir, repath + "/" + track_dir + "/" + track_fs[0])

if __name__ == "__main__":
    # Create genomedata archives for all replicates
    for rep in os.listdir("data/"):
        print("creating genomedata archive for " + rep)
        bigwig_to_bedgraph_dir("data/" + rep)
        create_genomedata_ct("data/" + rep, "data/sequence/hg38.fa")
        print("created genomedata archive for " + rep)