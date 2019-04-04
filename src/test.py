
import os.path
import sys
from src import reader

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
reader.searchDBlast('O14733')
reader.searchUniProt(reader.readDBlast())
data= reader.readUniProt()
go= reader.getGO(data)
kegg= reader.getKEGG(data)
reader.drawMultiPie([go, kegg])






