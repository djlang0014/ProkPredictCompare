#!/usr/local/bin/python3
import subprocess
import jinja2
import cgi
from pathlib import Path
import ParseOutputs
import mysql

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

#Removes files made from previous runs
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

#Save to database if new

#subprocess.run("rm -f FRAG*",shell=True,stdout=subprocess.PIPE)
#subprocess.run("rm -f glimmeroutput.*",shell=True,stdout=subprocess.PIPE)
#subprocess.run("rm -f MGAout",shell=True,stdout=subprocess.PIPE)
#subprocess.run("rm -f prodigaloutput.gbk",shell=True,stdout=subprocess.PIPE)
#subprocess.run("rm -f sequence.fna",shell=True,stdout=subprocess.PIPE)
#subprocess.run("rm -f genomic.gff",shell=True,stdout=subprocess.PIPE)

#subprocess.run("rm -f sequence.fna",shell=True,stdout=subprocess.PIPE)
#subprocess.run("rm -f genomic.gff",shell=True,stdout=subprocess.PIPE)

print("Content-type: text/html\n\n")
print(template.render(prodigalCount = [1,2,3,4,5,6,7]))