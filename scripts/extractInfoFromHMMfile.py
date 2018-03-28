#!/usr/bin/env python3
# vim: set fileencoding=utf-8 :
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package
# in LICENSE.txt.
"""Crawl through all hmm files stored in all subdirectories and extract the information stored in the annotation section"""

import os
from os.path import join
import fnmatch
import json


# allowed keywords to parse
#
# CAVEAT: Currently does parse only single-line content
keywords = ("HMMER",
            "NAME",
            "ACC",
            "DESC",
            "LENG",
            "MAXL",
            "ALPH",
            "RF",
            "MM",
            "CONS",
            "CS",
            "MAP",
            "DATE",
#           "COM", # can be multi-line
            "NSEQ",
            "EFFN",
            "CKSUM",
            "GA",
            "TC",
            "NC",
#           "STATS", # multi-line
#           "HMM",   # multi-line
#           "COMPO"  # multi-line
            )


def read_and_process_file(root: str, filename: str):
    """ Read HMM profile file and parse metadata into hmmannotation dictionary
    
        hmmannotation["NAME"]="abcd"
    
        All metadata, including directory of file, filename is stored in this dictionary,
        which then is used to generate the various output formats
        
        
        Arguments:
            root: Directory containing HMM profile file
            filename: Filename of HMM profile file
        
        Returns:
             hmmannotation dictionary"""

    #CAVEATS: Multi-line comments are ignored; HMM file parsing

    hmmannotation = {}
    with open(join(root, filename), "r") as f:
        hmmannotation["HMMDirectory"] = root
        hmmannotation["HMMfile"] = filename
        
        # read first line and parse HMMer profile file version; as there may be weird files, we have to check, whether they have the HMMer header
        readln = f.readline()
        if not readln in ("\n", "", "//"):
            (HMMerVersion, dummy) = readln.split(maxsplit=1)
            if not HMMerVersion.startswith("HMMER"):
                raise ValueError ("Error: " + filename + "does not start with HMMer version header line")
            
            
            hmmannotation["HMMER"] = HMMerVersion
        for readln in f:
            readln = readln.rstrip()
            if readln == "//":
                break
            (key, value) = readln.split(maxsplit=1)
            if key in keywords:
                hmmannotation[key] = value

    return hmmannotation


def write_table(directory: str, hmmannotation: dict):
    """Writes hmm-meta-tab.txt tab-separated file of annotation inferred from HMM file / other functions
    
        Arguments:
            Directory: directory to store the table in
            hmmannotation: hmmannotation dictionary
            
        Returns:
            None"""

    with open(join(directory, "hmm-meta-tab.txt"), "w") as f:
        f.writelines(["\t".join(hmmannotation) + "\n", "\t".join(hmmannotation.values()) + "\n"])



def write_json(directory: str, hmmannotation: dict):
    """Writes hmm-meta.json" JSON formatted file of annotation inferred from HMM file / other functions
    
        Arguments:
            Directory: directory to store the table in
            hmmannotation: hmmannotation dictionary
            
        Returns:
            None"""
    with open(join(directory, "hmm-meta.json"), "w") as f:
        json.dump(hmmannotation, f)


def infer_source_from_accession(hmmannotation: dict):
    """Infer data source from accession records, i.e. if accession is
    
    PFxxxxx ==> PFAM
    TIGRxxxx ==> TIGRFAMS
    
    and add this to hmmannotations dictionary
    
    Arguments:
        hmmannotations dictionary
        
    Returns:
        None"""

    HMMacc = hmmannotation.get("ACC", "")
    # Check for PFAM IDs
    if HMMacc.startswith("PF"):
        hmmannotation["SOURCE"] = "PFAM"

    elif HMMacc.startswith("TIGR"):
        hmmannotation["SOURCE"] = "TIGRFAMS"

#TODO: Use table from antiSMASH publication to refer to sources, such as BAGEL, CLUSEAN,...
    else:
        hmmannotation["SOURCE"] = "unknown"
#====================================================================
#main
#########################



# Crawl through all directories and identify hmm files
for root, dirs, files in os.walk('.'):
        # Skip the data directory as this directory will contain the collection of all hmm profiles that
        # are collected from the individual directories
        if "data" in dirs:
            dirs.remove("data")
        for filename in files:
            if fnmatch.fnmatch(join(root, filename), "*.hmm"):
                hmmannotation = read_and_process_file(root, filename)

                # Infer source
                infer_source_from_accession(hmmannotation)
                print("Directory:", root)

                # Write tabular metadata summary from HMM file
                write_table(root, hmmannotation)

                # Write JSON metadata file
                write_json(root, hmmannotation)
