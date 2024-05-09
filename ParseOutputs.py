
from pathlib import Path
import re
import statistics

def parseProdigal(genbankinfo):
    prodigalfile = Path('./prodigaloutput.gbk')
    prodigalinfo = []
    countinfo = []
    currentPrediction = ()
    previous = False
    first = True
    comp = False

    with open(prodigalfile, 'r') as file:
        for line in file:
            if first and line.startswith("     CDS"):
                location = re.search(r".*CDS             <([0-9]*)..([0-9]*)", line)
                currentPrediction = (location.group(1), location.group(2))
                first = False
                previous = True
                pass
            elif previous:
                conf = re.search(r".*conf=([0-9]*.[0-9]*)", line)
                currentPrediction = (currentPrediction[0], currentPrediction[1], conf.group(1))
                if comp:
                    prodigalinfo.append(currentPrediction)
                    comp = False
                else:
                    prodigalinfo.append(currentPrediction)
                currentPrediction = ()
                previous = False
            else:
                if line.startswith("     CDS"):
                    if line.__contains__("complement"):
                        location = re.search(r".*CDS             complement.([0-9]*)..([0-9]*)", line)
                        currentPrediction = (location.group(1), location.group(2))
                        comp = True
                    else:
                        location = re.search(r".*CDS             ([0-9]*)..([0-9]*)", line)
                        currentPrediction = (location.group(1), location.group(2))
                        #This deals with the last line which is different from the others
                        if not location:
                            location = re.search(r".*CDS             <([0-9]*)..>([0-9]*)", line)
                    previous = True

    #Get the number of exact matches, 5' matches, 3' matches, and no matches
    exactmatches = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in prodigalinfo if one[1] == two[0] and one[2] == two[1]]
    fiveonly = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in prodigalinfo if one[1] == two[0] and one[2] != two[1]]
    threeonly = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in prodigalinfo if one[1] != two[0] and one[2] == two[1]]
    nomatch = [one for one in genbankinfo if all(one[1] != two[0] and one[2] != two[1] for two in prodigalinfo)]

    countinfo.append((len(genbankinfo)))
    countinfo.append((len(prodigalinfo)))    
    countinfo.append(len(exactmatches))
    countinfo.append(len(fiveonly))
    countinfo.append(len(threeonly))
    countinfo.append(len(nomatch))
    #for CDS medians
    lengths = [int(two[1]) - int(two[0]) for two in prodigalinfo if two[0] != '' and two[1] != '']
    countinfo.append(statistics.median(lengths))

    finallist = sorted(exactmatches + fiveonly + threeonly + nomatch, key=lambda x: x[0])
    return finallist, countinfo

def parseGlimmer(genbankinfo):
    glimmerfile = Path('./glimmeroutput.predict')
    #Get the number of exact matches, 5' matches, 3' matches, and no matches
    glimmerinfo = []
    with open(glimmerfile, 'r') as file:
        for line in file:
            if line.startswith(">"):
                continue
            
            #captures the location of the gene, skips reading frame info at start
            location = re.search(r"\s*\d+\s+(\d+)\s+(\d+)", line)
            if location:
                glimmerinfo.append((location.group(1), location.group(2)))

    exactmatches = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in glimmerinfo if one[1] == two[0] and one[2] == two[1]]
    fiveonly = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in glimmerinfo if one[1] == two[0] and one[2] != two[1]]
    threeonly = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in glimmerinfo if one[1] != two[0] and one[2] == two[1]]
    nomatch = [one for one in genbankinfo if all(one[1] != two[0] and one[2] != two[1] for two in glimmerinfo)]

    countinfo = []
    countinfo.append((len(genbankinfo)))
    countinfo.append((len(glimmerinfo)))    
    countinfo.append(len(exactmatches))
    countinfo.append(len(fiveonly))
    countinfo.append(len(threeonly))
    countinfo.append(len(nomatch))
    #for CDS medians
    lengths = [int(two[1]) - int(two[0]) for two in glimmerinfo if two[0] != '' and two[1] != '']
    countinfo.append(statistics.median(lengths))

    finallist = sorted(exactmatches + fiveonly + threeonly + nomatch, key=lambda x: x[0])
    return finallist, countinfo

def parseMGA(genbankinfo):
    mgafile = Path('./MGAout')
    mgainfo = []

    with open(mgafile, 'r') as file:
        for line in file:
                if line.startswith("#"):
                    #This skips the inital header lines of the gff file
                    continue

                #features is a list of columns in the current line of the gff file.
                features = line.strip().split("\t") #strips unnecessary white space and splits the line by tabs
                chromStart = features[1]
                chromEnd = features[2]
                finalinfo = (chromStart, chromEnd)
                mgainfo.append(finalinfo)

    exactmatches = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in mgainfo if one[1] == two[0] and one[2] == two[1]]
    fiveonly = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in mgainfo if one[1] == two[0] and one[2] != two[1]]
    threeonly = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in mgainfo if one[1] != two[0] and one[2] == two[1]]
    nomatch = [one for one in genbankinfo if all(one[1] != two[0] and one[2] != two[1] for two in mgainfo)]

    countinfo = []
    countinfo.append((len(genbankinfo)))
    countinfo.append((len(mgainfo)))    
    countinfo.append(len(exactmatches))
    countinfo.append(len(fiveonly))
    countinfo.append(len(threeonly))
    countinfo.append(len(nomatch))
    #for CDS medians
    lengths = [int(two[1]) - int(two[0]) for two in mgainfo if two[0] != '' and two[1] != '']
    countinfo.append(statistics.median(lengths))

    finallist = sorted(exactmatches + fiveonly + threeonly + nomatch, key=lambda x: x[0])
    return finallist, countinfo

def parseFGS(genbankinfo):
    fgsfile = Path('./FRAGout.gff')
    fgsinfo = []
    #Get the number of exact matches, 5' matches, 3' matches, and no matches
    with open(fgsfile, 'r') as file:
        for line in file:
                if line.startswith("#"):
                    #This skips the inital header lines of the gff file
                    continue

                #features is a list of columns in the current line of the gff file.
                features = line.strip().split("\t") #strips unnecessary white space and splits the line by tabs


                chromStart = features[3]
                chromEnd = features[4]
                finalinfo = (chromStart, chromEnd)
                fgsinfo.append(finalinfo)

    exactmatches = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in fgsinfo if one[1] == two[0] and one[2] == two[1]]
    fiveonly = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in fgsinfo if one[1] == two[0] and one[2] != two[1]]
    threeonly = [(one[0], one[1], one[2], two[0], two[1]) for one in genbankinfo for two in fgsinfo if one[1] != two[0] and one[2] == two[1]]
    nomatch = [one for one in genbankinfo if all(one[1] != two[0] and one[2] != two[1] for two in fgsinfo)]

    countinfo = []
    countinfo.append((len(genbankinfo)))
    countinfo.append((len(fgsinfo)))    
    countinfo.append(len(exactmatches))
    countinfo.append(len(fiveonly))
    countinfo.append(len(threeonly))
    countinfo.append(len(nomatch))
    #for CDS medians
    lengths = [int(two[1]) - int(two[0]) for two in fgsinfo if two[0] != '' and two[1] != '']
    countinfo.append(statistics.median(lengths))

    finallist = sorted(exactmatches + fiveonly + threeonly + nomatch, key=lambda x: x[0])    
    return finallist, countinfo