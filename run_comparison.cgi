#!/usr/local/bin/python3
import subprocess
import jinja2
import cgi
from pathlib import Path
import ParseOutputs
import mysql.connector
#import output
import json

#CGI form
form = cgi.FieldStorage()
accession = form.getvalue('accession')
accession = "GCF_000008865"

#This line tells the template loader where to search for template files
templateLoader = jinja2.FileSystemLoader(searchpath="./templates")

#This creates your environment and loads a specific template
env = jinja2.Environment(loader=templateLoader)
template = env.get_template('output.html')

conn = mysql.connector.connect(user='dlang15', password='Bioinformatics!',
                               host='localhost', database='dlang15_final')
curs = conn.cursor()
subprocess.run("touch info.txt",shell=True,stdout=subprocess.PIPE)
subprocess.run("chmod 777 info.txt",shell=True,stdout=subprocess.PIPE)
store = open("info.txt", "w")
store.write(accession)
store.close()
qry = """
SELECT GenesNum, ProdigalInfo, GlimmerInfo, FGSInfo, MGAInfo FROM Results
WHERE AccessionNum = %s
"""
curs.execute(qry, (accession,))
results = curs.fetchall()

#This only needs to run if the accession number is not in the database
if not results:
    accessionNumStore = accession
    #Removes files made from previous runs
    #I add perms to final because they disappeared once and I debugged for multiple hours
    subprocess.run("chmod 777 /var/www/html/dlang15/final",shell=True,stdout=subprocess.PIPE)
    subprocess.run("rm -f FRAG*",shell=True,stdout=subprocess.PIPE)
    subprocess.run("rm -f glimmeroutput.*",shell=True,stdout=subprocess.PIPE)
    subprocess.run("rm -f MGAout",shell=True,stdout=subprocess.PIPE)
    subprocess.run("rm -f prodigaloutput.gbk",shell=True,stdout=subprocess.PIPE)
    subprocess.run("rm -f sequence.fna",shell=True,stdout=subprocess.PIPE)
    subprocess.run("rm -f genomic.gff",shell=True,stdout=subprocess.PIPE)
    subprocess.run("rm -rf ncbi_dataset",shell=True,stdout=subprocess.PIPE)
    subprocess.run("rm -rf ncbi_dataset.zip",shell=True,stdout=subprocess.PIPE)

    #Retrieve the genome data from NCBI
    retrievecommand = f"./datasets download genome accession {accession} --include gff3,genome"
    subprocess.run(retrievecommand, shell=True,stdout=subprocess.PIPE)
    subprocess.run("chmod 777 ncbi_dataset.zip", shell=True,stdout=subprocess.PIPE)
    subprocess.run("unzip ncbi_dataset.zip", shell=True,stdout=subprocess.PIPE)
    subprocess.run("mv ./ncbi_dataset/data/GC*/* .", shell=True,stdout=subprocess.PIPE)
    subprocess.run("mv GC* sequence.fna", shell=True,stdout=subprocess.PIPE)
    subprocess.run("rm -rf ncbi_dataset", shell=True,stdout=subprocess.PIPE)
    subprocess.run("rm -f ncbi_dataset.zip", shell=True,stdout=subprocess.PIPE)
    subprocess.run("rm -f README.md",shell=True,stdout=subprocess.PIPE)
    subprocess.run("chmod 777 genomic.gff", shell=True,stdout=subprocess.PIPE)
    subprocess.run("chmod 777 sequence.fna",shell=True,stdout=subprocess.PIPE)


    gffFile = Path('./genomic.gff')
    sequenceFile = Path('./sequence.fna')

    genbankinfo = []
    genenum = 1

    #This is ignoring incomplete features
    #Gets the cds information from the gff file for the comparisons
    with open(gffFile, 'r') as file:
        for line in file:
                if line.startswith("#"):
                    #This skips the inital header lines of the gff file
                    continue

                #features is a list of columns in the current line of the gff file.
                features = line.strip().split("\t") #strips unnecessary white space and splits the line by tabs

                #We are only interested in gene features.
                if features[2] == "gene":
                    accessionNum = f"{accession}-{genenum}"
                    chromStart = features[3]
                    chromEnd = features[4]
                    finalinfo = (accessionNum, chromStart, chromEnd)
                    genbankinfo.append(finalinfo)
                    genenum += 1
                    #This finds the name of the gene from the attributes column of the gff file
                else:
                    continue

    #Set the commands for the programs
    commandMGA = f"/var/www/html/dlang15/final/mga_linux_ia64 {sequenceFile}"                                   
    commandGlimmer = f"/var/www/html/dlang15/final/glimmer3.02/bin/g3-from-scratch.csh {sequenceFile} glimmeroutput"
    commandFGS = f"/var/www/html/dlang15/final/FragGeneScan-master/run_FragGeneScan.pl -genome={sequenceFile} -out=FRAGout -complete=1 -train=complete"
    commandProdigal = f"/var/www/html/dlang15/final/prodigal -i {sequenceFile} -o prodigaloutput.gbk"

    #Run the programs + give necessary perms
    subprocess.run(commandProdigal,shell=True,stdout=subprocess.PIPE)
    subprocess.run(commandGlimmer, shell=True,stdout=subprocess.PIPE)
    subprocess.run(commandFGS, shell=True,stdout=subprocess.PIPE)
    subprocess.run("chmod 777 FRAG*",shell=True,stdout=subprocess.PIPE)
    subprocess.run("chmod 777 glimmeroutput.*",shell=True,stdout=subprocess.PIPE)
    subprocess.run("chmod 777 prodigaloutput.gbk",shell=True,stdout=subprocess.PIPE)
    subprocess.run("touch MGAout",shell=True,stdout=subprocess.PIPE)
    subprocess.run("chmod 777 MGAout",shell=True,stdout=subprocess.PIPE)

    with open("MGAout", "w") as outfile:
            subprocess.run(commandMGA, shell=True, stdout=outfile)

    prodigalFinal, prodigalCount = ParseOutputs.parseProdigal(genbankinfo)
    glimmerFinal, glimmerCount = ParseOutputs.parseGlimmer(genbankinfo)
    fgsFinal, fgsCount = ParseOutputs.parseFGS(genbankinfo)
    mgaFinal, mgaCount = ParseOutputs.parseMGA(genbankinfo)

    #Save to database
    qry = """
    INSERT INTO ProdigalPredictions (AccessionNum, PredictionID, Start, Stop)
    VALUES (%s, %s, %s, %s)
    """
    for i, (prodigalaccession, prodigal_start, prodigal_stop, *_) in enumerate(prodigalFinal, start=1):
        #remove trailing number from accession
        prodigalaccessionFix = prodigalaccession.split("-")[0]
        predictionID = f"Gene_{i}"
        prodigal_prediction = (prodigalaccessionFix, predictionID, prodigal_start, prodigal_stop)
        curs.execute(qry, prodigal_prediction)
    conn.commit()
    
    qry = """
    INSERT INTO FGSPredictions (AccessionNum, PredictionID, Start, Stop)
    VALUES (%s, %s, %s, %s)
    """
    for i, (fgsaccession, FGS_start, FGS_stop, *_) in enumerate(fgsFinal, start=1):
        fgsaccessionFix = fgsaccession.split("-")[0]
        predictionID = f"Gene_{i}"
        fgs_prediction = (fgsaccessionFix, predictionID, FGS_start, FGS_stop)
        curs.execute(qry, fgs_prediction)
    conn.commit()

    qry = """
    INSERT INTO GlimmerPredictions (AccessionNum, PredictionID, Start, Stop)
    VALUES (%s, %s, %s, %s)
    """
    for i, (glimmeraccession, glimmer_start, glimmer_stop, *_) in enumerate(glimmerFinal, start=1):
        glimmeraccessionFix = glimmeraccession.split("-")[0]
        predictionID = f"Gene_{i}"
        glimmer_prediction = (glimmeraccessionFix, predictionID, glimmer_start, glimmer_stop)
        curs.execute(qry, glimmer_prediction)
    conn.commit()

    qry = """
    INSERT INTO MGAPredictions (AccessionNum, PredictionID, Start, Stop)
    VALUES (%s, %s, %s, %s)
    """
    for i, (mgaaccession, MGA_start, MGA_stop, *_) in enumerate(mgaFinal, start=1):
        mgaaccessionFix = mgaaccession.split("-")[0]
        predictionID = f"Gene_{i}"
        MGA_prediction = (mgaaccessionFix, predictionID, MGA_start, MGA_stop)
        curs.execute(qry, MGA_prediction)
    conn.commit()

    qry = """
    INSERT INTO Results (AccessionNum, GenesNum, ProdigalInfo, GlimmerInfo, FGSInfo, MGAInfo)
    VALUES (%s, %s, %s, %s, %s, %s)
    """
    
    GenesNum = len(genbankinfo)
    ProdigalInfo = json.dumps(prodigalCount)
    GlimmerInfo = json.dumps(glimmerCount)
    FGSInfo = json.dumps(fgsCount)
    MGAInfo = json.dumps(mgaCount)
    data_results = (accessionNumStore, GenesNum, ProdigalInfo, GlimmerInfo, FGSInfo, MGAInfo)
    curs.execute(qry, data_results)
    conn.commit()

    #subprocess.run("rm -f FRAG*",shell=True,stdout=subprocess.PIPE)
    #subprocess.run("rm -f glimmeroutput.*",shell=True,stdout=subprocess.PIPE)
    #subprocess.run("rm -f MGAout",shell=True,stdout=subprocess.PIPE)
    #subprocess.run("rm -f prodigaloutput.gbk",shell=True,stdout=subprocess.PIPE)
    #subprocess.run("rm -f sequence.fna",shell=True,stdout=subprocess.PIPE)
    #subprocess.run("rm -f genomic.gff",shell=True,stdout=subprocess.PIPE)

    #subprocess.run("rm -f sequence.fna",shell=True,stdout=subprocess.PIPE)
    #subprocess.run("rm -f genomic.gff",shell=True,stdout=subprocess.PIPE)
else:
    #We have the results from the check query
    GenesNum = results[0][0]
    FGSInfo = json.loads(results[0][1])
    GlimmerInfo = json.loads(results[0][2])
    ProdigalInfo = json.loads(results[0][3])
    MGAInfo = json.loads(results[0][4])

    #These are now holding a list of # of gff genes, # of predictions, # of exact matches, # of 5' matches, # of 3' matches, # of no matches


print("Content-type: text/html\n\n")
print(template.render(prodigalCount = [GenesNum,FGSInfo,GlimmerInfo,ProdigalInfo,MGAInfo,6,7]))

conn.close()