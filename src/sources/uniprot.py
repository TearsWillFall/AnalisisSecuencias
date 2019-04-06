import datetime
import pandas

from yaspin import yaspin
from bioservices import UniProt


#
# Yield successive n-sized chunks from l.
#
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


#
# Searches for entries in uniprot, returns merged results.
#
def search_uniprot(genes, columns):

    start_time = datetime.datetime.now()

    with yaspin(text="Performing UniProt search...", color="cyan") as sp:
        uniprot = UniProt(verbose=False)
        raw_data = ''
        headers = True

        for chunk in chunks(genes, 190):
            gene_search = "+OR+".join(list(chunk))
            new_data = uniprot.search(gene_search, frmt="tab", columns=f"entry name, {','.join(columns)}")

            # Removes first line if this is second or next batches.
            if not headers:
                new_data = new_data.split("\n")[1]

            raw_data = raw_data + "\n" + new_data
            headers = False

        data = pandas.read_csv(pandas.compat.StringIO(raw_data), sep="\t")

        time_diff = (datetime.datetime.now() - start_time).total_seconds()

        sp.text = f"Performing UniProt Search => Task done in {time_diff} seconds."
        sp.ok("✔")

        return data
