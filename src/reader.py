
from bioservices import UniProt
from bioservices import KEGG
import time
from Bio.Blast.Applications import NcbideltablastCommandline
from yaspin import yaspin
from pygal.style import Style
import pygal
import webbrowser


def searchDBlast(query="Q03963",evalue=0.001):
    '''
    This function performs a remote Delta Blast search. Requires previous installation of Blast+ (v 2.9.0+)
    :param query: A string with the name of the file with the query.
    :param evalue: A float with evalue threshold.
    :return: A tab delimited file (filename="search.tab") with Delta Blast search results.
    '''
    start_time=time.clock()
    with yaspin(text="Performing Delta Blast Search...", color="cyan") as sp:
        cline = NcbideltablastCommandline(query=query, db="nr", evalue=evalue, remote=True, out="search.tab", outfmt=7)
        cline()
        sp.text="Performing Delta Blast Search=> Task done in "+str(run_time(start_time,4))+"s"
        sp.ok("✔")
def readDBlast ():
    '''
    This function reads the Tab delimited file generated during a Delta Blast search and generates a dictionary.
    :return: A Dictionary with Delta Blast search data.
    '''
    start_time = time.clock()
    with yaspin(text="Reading Delta Blast Search...", color="cyan") as sp:
        d={}
        with open ("search.tab") as file:
            for line in file:
                line=line.rstrip()

                if line.find("Fields")!=-1:
                    keys = line.split(": ")
                    keys = keys[1].split(", ")

                elif line[0] !="#":
                    line=line.split("\t")
                    i=0
                    for key in keys:
                        if key not in d:
                            d[key]=[line[i]]
                        else:
                            d[key].append(line[i])
                        i+=1

        file.close()
        ##os.remove("search.tab")
        sp.text = "Reading Delta Blast Search=> Task done in " + str(run_time(start_time,4)) + "s"
        sp.ok("✔")
        return d

def searchUniProt(d):
    '''
    This function uses the Uniprot IDs of the Delta Blast search saved in a dictionary and returns a tab delimited file.
    :param d: A dictionary with Delta Blast search data.
    :return: A Tab Delimited File with UniProt search results.
    '''
    search_uniprot=open("search_UniProt.tab","w")
    start_time = time.clock()
    with yaspin(text="Performing UniProt Search...", color="cyan") as sp:
        uni = UniProt(verbose=False)
        search = "+OR+".join(d['subject acc.ver'][0:200])
        data = uni.search(search, frmt="tab", columns="entry name,length,id, go, database(kegg)")

        reps=((len(d['subject acc.ver'])/200.0)-1)
        reps_done=1
        while reps >0:
            next_search = "+OR+".join(d['subject acc.ver'][(200*reps_done+1):200*(reps_done+1)])
            result = uni.search(next_search, frmt="tab", columns="entry name,length,id, go, database(kegg)")
            result = result.split("\n")
            header=True
            data = data.rstrip()
            for line in result:
                if not header:
                    data +="\n"+line
                header=False
            reps_done += 1
            reps -= 1
        search_uniprot.write(data)
        sp.text = "Performing UniProt Search=> Task done in " + str(run_time(start_time,4)) + "s"
        sp.ok("✔")
def readUniProt():
    '''
    This function reads the Tab delimited file generated in the UniProt search and generates a dictionary of the data.
    :return: A dictionary with UniProt search data.
    '''
    start_time = time.clock()
    with yaspin(text="Reading UniProt Search...", color="cyan") as sp:
        header = True
        d={}
        with open("search_UniProt.tab") as data:
            for line in data:
                line = line.strip()
                line=line.split("\t")
                if header:
                    keys =line
                    header=False
                else:
                    i=0
                    for key in keys:

                        try:
                            if line[i].find("; ") != -1:
                                out = line[i].split("; ")
                            else:
                                out=line[i][:-1]
                                out=out.split(";")
                        except  IndexError:
                            out="-"
                        if key not in d:
                            d[key] = [out]
                        else:
                            d[key].append(out)
                        i += 1
            sp.text = "Reading UniProt Search=> Task done in " + str(run_time(start_time,4)) + "s"
            sp.ok("✔")
            return d

def getGO(d):
    '''
    This function retrieves GO Annotations count from a dictionary of UniProt search data.
    :param d: A dictionary with UniProt search data.
    :return: A dictionary with GO Annotation count.
    '''
    start_time = time.clock()
    with yaspin(text="Retrieving GO Annotations...", color="cyan") as sp:
        count_Go={}
        for list in d["Gene ontology (GO)"]:
            for elem in list:
                if elem !="-":
                    try:
                        count_Go[elem] +=1
                    except KeyError:
                        count_Go[elem] = 1
        sp.text = "Retrieving GO Annotations=> Task done in " + str(run_time(start_time,4)) + "s"
        sp.ok("✔")
        print(count_Go)
        return count_Go


def getKEGG(d):
    '''This function retrieves KEGG Pathways count using KEGG IDs obtained from a UniProt search dictionary.
   :param d: A dictionary with UniProt search data.
   :return: A dictionary with KEGG pathways count.
   '''
    start_time = time.clock()
    with yaspin(text="Retrieving KEGG Annotations...", color="cyan") as sp:
        count_KEGG={}
        for list in d["Cross-reference (kegg)"]:
            for elem in list:
                if elem !="-":
                    try:
                        count_KEGG[elem] +=1
                    except KeyError:
                        count_KEGG[elem] = 1
        path = KEGG()
        dic_Paths={}
        count_Paths={}
        print (count_KEGG)
        for key in count_KEGG:
            res=key.split(":")

            try:
                for k,val in path.get_pathway_by_gene(res[1], res[0]).items():
                    name = ":".join([k, val])
                    if key not in dic_Paths:
                        dic_Paths[key] =[name]
                    else:
                        dic_Paths[key].append(name)
                    path_id = split(k)
                    rename="map"+path_id
                    rename=":".join([val, rename])
                    try:
                        count_Paths[rename] += count_KEGG[key]
                    except KeyError:
                        count_Paths[rename] = count_KEGG[key]
            except AttributeError:
                break
        sp.text = "Retrieving KEGG Annotations=> Task done in " + str(run_time(start_time,4)) + "s"
        sp.ok("✔")
        return count_Paths

def drawMultiPie (dbs,limit = 0):
    '''
    This function draws a Multi Series Pie chart for KEGG and GO terms
    :param dbs: List of Dictionaries with Annotation counts for each Annotation Database.
    :param limit: Integer value of count to filter Annotations .
    :return: HTML Multi Series Pie
    '''
    custom_style = Style(
        opacity='0.8',
        opacity_hover='0.5',
        title_font_size=36,
        inner_radius = 0.75
    )
    pie_map = pygal.Pie(height=400,tooltip_border_radius=10, style=custom_style)
    pie_map.title = "Functional Annotations"
    pie_map_other = pygal.Treemap(height=400, tooltip_border_radius=10, style=custom_style)
    pie_map_other.title = "Other"
    html_file = open("merged.html", 'w')
    html_file.write("<html><head>…</head><body>" + "\n")
    draw_pie=False
    draw_pie_other=False
    for db in dbs:
        values = []
        values_other= []
        go=False
        limited=False
        if [*db][0].find("[")!=-1:
            go=True
        other=0
        for key in sorted(db, key=db.get, reverse=True):
            if db[key] <=limit:
                break
            if db[key] <=15:
                other+=db[key]
                limited=True
            if go:
                name, id = key.split(" [")
                id = id.replace(" [", "").replace("]", "")
                name = name.capitalize()
                link = 'https://www.ebi.ac.uk/QuickGO/term/' + id
                db_name = "GO"
                color='rgba(255, 45, 20, .6)'
            else:
                name, id = key.split(":")
                link = 'https://www.genome.jp/dbget-bin/www_bget?' + id
                db_name = "KEGG"
                color='rgba(68, 108, 179, .6)'

            if limited:
                values_other.append({'value': db[key], 'label': name, 'xlink': {'href': link}, 'color': color})
            else:
                values.append({'value': db[key], 'label': name, 'xlink': {'href': link},'color':color})
        if other:
            values.append({'value': other, 'label': 'Other', 'color': color})
        if values != []:
            draw_pie=True
            pie_map.add(db_name, values)
        if values_other != []:
            draw_pie_other=True
            pie_map_other.add(db_name, values_other)
        else:
            print ("No annotations were found.")
    if draw_pie:
        pie_map.render_to_file('graph1.svg')
        html_file.write("      <object type=\"image/svg+xml\" data=\"graph1.svg\"></object>" + "\n")
        if draw_pie_other:
            pie_map_other.render_to_file('graph2.svg')
            html_file.write("      <object type=\"image/svg+xml\" data=\"graph2.svg\"></object>" + "\n")
        html_file.write("</body></html>")
        webbrowser.open("merged.html")

def split (s):
    '''
    This function removes the organism identifier in KEGG IDs.
    :param s: A string with the KEGG ID to split.
    :return: A string with the numerical KEGG ID.
    '''
    head=s.rstrip("0123456789")
    tail=s[len(head):]
    return tail
def run_time (start_time,r):
    '''
    This function measures the run time of another function.
    :param start_time: A float with the start time of a function
    :param r: A integer value by which
    :return:
    '''
    end_time=time.clock()-start_time
    end_time=round(end_time,r)
    return end_time
