from bioservices import UniProt
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
help(NCBIWWW.qblast)
input_sequence = "PIO39916.1"

delta_blast = NCBIWWW.qblast(program = "blastp", database = "nr", sequence=input_sequence, service = "del", entrez_query = "all[filter] NOT predicted[title]")



blast_record = NCBIXML.read(delta_blast)
proteinIds=[]
for alignment in blast_record.alignments:
    proteinIds.append(alignment.hit_id.rsplit("|",2)[-2])
    for hsp in alignment.hsps:
            print("****Alignment****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print("identities:", (hsp.identities/ hsp.align_length*100),"%")
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")

uni = UniProt(verbose=False)
search= "+OR+".join(proteinIds)
data = uni.search(search, frmt="tab", columns="entry name,length,id, go, database(kegg)")
print(data)
print(proteinIds)