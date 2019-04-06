import datetime
from yaspin import yaspin

import pandas

from Bio.Blast.Applications import NcbideltablastCommandline


#
# Performs a Delta Blast query using command line deltablast.
#
def do_dblast_query(query, evalue=0.001):

    start_time = datetime.datetime.now()

    with yaspin(text="Performing Delta Blast search...", color="cyan") as sp:
        cline = NcbideltablastCommandline(query=query,
                                          db="nr",
                                          evalue=evalue,
                                          remote=True,
                                          out=f"{query}.tsp",
                                          outfmt=7)
        cline()

        time_diff = (datetime.datetime.now() - start_time).total_seconds()

        sp.text = f"Performing Delta Blast search => Task done in {time_diff}."
        sp.ok("✔")


#
# Reads a Delta Blast results file.
#

def read_dblast_file(file_name):

    with yaspin(text="Reading Delta Blast search...", color="cyan") as sp:
        sp.ok("✔")

        start_time = datetime.datetime.now()
        data = ""

        try:
            # Read column names from comments in file.
            with open(file_name) as handle:
                for line in handle:
                    if line.startswith("# Fields: "):
                        names = list(map(str.strip, line.strip("# Fields: ").rstrip("\n").split(",")))
                        break

            # Read file as Data.Frame
            data = pandas.read_csv(file_name, sep="\t", comment="#", names=names)

        except Exception as e:
            sp.text = f"Reading Delta Blast search => Error: {e}"
            sp.color = "red"
            sp.fail("✗")
            exit(1)

        time_diff = (datetime.datetime.now() - start_time).total_seconds()

        sp.text = f"Reading Delta Blast search => Task done in {time_diff} seconds."
        sp.ok("✔")

        return data
