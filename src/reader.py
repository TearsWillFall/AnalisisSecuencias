
from bioservices import UniProt
from bioservices import KEGG
import os
import sys
import time
from Bio.Blast.Applications import NcbideltablastCommandline
from yaspin import yaspin
from pygal.style import Style
import pygal
from lxml.html import open_in_browser
import webbrowser


def searchDBlast(query="Q03963",evalue=0.001):
    start_time=time.clock()
    with yaspin(text="Performing Delta Blast Search...", color="cyan") as sp:
        cline = NcbideltablastCommandline(query=query, db="nr", evalue=evalue, remote=True, out="search.tab", outfmt=7)
        cline()
        sp.text="Performing Delta Blast Search=>Task done in "+str(run_time(start_time,4))+"s"
        sp.ok("✔")
def readDBlast ():
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
        end_time = time.clock() - start_time
        sp.text = "Reading Delta Blast Search=>Task done in " + str(run_time(start_time,4)) + "s"
        sp.ok("✔")
        return d

def searchUniProt(d):
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
        end_time = time.clock() - start_time
        sp.text = "Performing UniProt Search=>Task done in " + str(run_time(start_time,4)) + "s"
        sp.ok("✔")
        return data

def readUniProt(data):
    start_time = time.clock()
    with yaspin(text="Reading UniProt Search...", color="cyan") as sp:
        header = True
        d={}
        data=data.rstrip()
        data=data.split("\n")

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
                        if line[i].find("; ")!=-1:
                            out=line[i].split("; ")
                        else:
                            out=[line[i].replace(";","")]
                    except  IndexError:
                        out="-"
                    if key not in d:
                        d[key] = [out]
                    else:
                        d[key].append(out)
                    i += 1
        end_time = time.clock() - start_time
        sp.text = "Reading UniProt Search=>Task done in " + str(run_time(start_time,4)) + "s"
        sp.ok("✔")
        return d

def getGO(d):
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
        end_time = time.clock() - start_time
        sp.text = "Retrieving GO Annotations=>Task done in " + str(run_time(start_time,4)) + "s"
        sp.ok("✔")
        return count_Go


def getKEGG(d):
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
        for key in count_KEGG:
            res=key.split(":")
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
        end_time = time.clock() - start_time
        sp.text = "Retrieving KEGG Annotations=>Task done in " + str(run_time(start_time,4)) + "s"
        sp.ok("✔")
        return count_Paths

def drawMultiPie (dbs,limit = 30):
    custom_style = Style(
        opacity='0.8',
        opacity_hover='0.5',
        title_font_size=36,
        inner_radius = 0.75
    )
    treemap = pygal.Pie(height=400,tooltip_border_radius=10, style=custom_style)
    treemap.title = "Functional Annotations"
    treemap_other = pygal.Pie(height=400, tooltip_border_radius=10, style=custom_style)
    treemap_other.title = "Other"
    html_file = open("merged.html", 'w')
    html_file.write("<html><head>…</head><body>" + "\n")

    for db in dbs:
        values = []
        values_other= []
        go=False
        limited=False
        if [*db][0].find("[")!=-1:
            go=True
        other=0
        for key in sorted(db, key=db.get, reverse=True):
            if db[key] <= limit:
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
            values.append({'value': other, 'label': 'Other', 'xlink': {'href': link}, 'color': color})
        if values_other != []:
            treemap_other.add(db_name, values_other)
        if values != []:
            treemap.add(db_name, values)
        else:
            print ("No annotations were found.")
    treemap.render_to_file('graph1.svg')
    html_file.write("      <object type=\"image/svg+xml\" data=\"graph1.svg\"></object>" + "\n")
    treemap_other.render_to_file('graph2.svg')
    html_file.write("      <object type=\"image/svg+xml\" data=\"graph2.svg\"></object>" + "\n")
    html_file.write("</body></html>")
    webbrowser.open("merged.html")

def split (s):
    head=s.rstrip("0123456789")
    tail=s[len(head):]
    return tail
def run_time (start_time,r):
    end_time=time.clock()-start_time
    end_time=round(end_time,r)
    return end_time
