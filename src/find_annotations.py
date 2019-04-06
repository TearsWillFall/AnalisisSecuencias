#!python3

import argparse
import os
import webbrowser

from sources import deltablast
from sources import uniprot
from sources import kegg
from sources import go

from plots import multi_pie
from utils.tools import select_data

###################################################################################################
# Argument parser.
parser = argparse.ArgumentParser(description="Find similar proteins using Delta Blast and suggests\
                                              GO annotations and KEGG routes associated to them.")
_ = parser.add_argument("protein",
                        action="store",
                        nargs="?",
                        help="Protein to analyze")

_ = parser.add_argument("min_identity",
                        action="store",
                        type=int,
                        help="Minimum identity score to consider (Default: 75)")

_ = parser.add_argument("min_support",
                        action="store",
                        type=int,
                        help="Minimum support degree to consider (amount of sequences sharing an\
                              annotation (Default: 10)")

_ = parser.add_argument("-f",
                        "--file",
                        action="store",
                        help="Specify dblast output file with list of similar proteins")

args = parser.parse_args()

if args.protein is None and args.file is None:
    parser.error("Specify either a protein accession or a dblast output file.")
    exit(1)
###################################################################################################


if args.file:
    print(f"* Running with file [{args.file}]", end=" ")
    name = os.path.splitext(os.path.basename(args.file))[0]

elif args.protein:
    print(f"* Running with accession [{args.protein}]", end=" ")
    name = args.protein

print(f"(Min identity score: [{args.min_identity}], Min support score: [{args.min_support}]).")


###################################################################################################
# Process start

# Read data from file or performs a delta blast.
if args.file:
    dblast_data = deltablast.read_dblast_file(args.file)

elif args.protein:
    dblast_data = deltablast.do_dblast_query(args.protein)

# Selects rows with given identity score.
selected_data = select_data(dblast_data, "% identity", args.min_identity)

# Find proteins in uniProt, extract GO list and KEGG accession.
uniprot_data = uniprot.search_uniprot(selected_data["subject acc.ver"], ["go", "database(kegg)"])

# Accumulate GO data.
go_list = go.accumulate_go(uniprot_data['Gene ontology (GO)'])

# Finds KEGG accessions, extract description and ID.
kegg_list = kegg.search_kegg(uniprot_data["Cross-reference (kegg)"])

# Draw charts.
multi_pie.plot(go_list, kegg_list, args.min_support, name)

webbrowser.open(f"{name}.html")
