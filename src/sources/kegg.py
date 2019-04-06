import datetime
import pandas
import re

from yaspin import yaspin
from bioservices import KEGG


#
# Gets a list of KEGG pathways.
#
def search_kegg(accessions):
    start_time = datetime.datetime.now()

    with yaspin(text="Retrieving KEGG annotations...", color="cyan") as sp:
        raw_data = ""
        for accession in accessions.dropna():
                path = KEGG()
                res = accession.split(":")

                for k, val in path.get_pathway_by_gene(res[1], res[0]).items():
                    _id = re.search("\d+", k).group(0)
                    raw_data = f"{raw_data}map{_id}\t\"{val}\"\n"

        kegg = pandas.read_csv(pandas.compat.StringIO(raw_data), sep="\t", header=None)
        kegg.columns = ["accession", "description"]

        # Add column of counts.
        kegg["count"] = kegg.groupby("accession")["accession"].transform("count")
        kegg = (
            kegg.drop_duplicates(subset="accession")
            .sort_values(by="count", ascending=False)
            .reset_index(drop=True)
        )

        time_diff = (datetime.datetime.now() - start_time).total_seconds()

        sp.text = f"Retrieving KEGG annotations => Task done in {time_diff} seconds."
        sp.ok("✔")

        return kegg
