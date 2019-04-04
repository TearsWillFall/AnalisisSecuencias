
import os.path
import sys
from src import reader

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

##reader.readDBlast()
##Para hacer busquedas en DBlast pero tarda un huevo asi que no la hago.

data= reader.readUniProt(reader.searchUniProt(reader.readDBlast()))
go= reader.getGO(data)
kegg= reader.getKEGG(data)
reader.drawMultiPie([go, kegg])






