#!/usr/local/bin/python3

import mysql.connector
import jinja2
import cgi


#CGI form
form = cgi.FieldStorage()
accession = form.getvalue('accession')

#This line tells the template loader where to search for template files
templateLoader = jinja2.FileSystemLoader(searchpath="." )

#This creates your environment and loads a specific template
env = jinja2.Environment(loader=templateLoader)
template = env.get_template('search.html')

conn = mysql.connector.connect(user='dlang15', password='Bioinformatics!',
                               host='localhost', database='dlang15_chado')
curs = conn.cursor()

qry = """
SELECT f.uniquename, product.value, location.fmin, location.fmax, location.strand 
FROM feature f 
JOIN cvterm geneproduct ON f.type_id=geneproduct.cvterm_id 
JOIN featureprop product ON f.feature_id=product.feature_id 
JOIN featureloc location ON f.feature_id=location.feature_id 
JOIN cvterm geneproductname ON product.type_id=geneproductname.cvterm_id
WHERE uniquename = %s AND geneproductname.name = 'gene_product_name'
"""

curs.execute(qry, (accession,))
result = curs.fetchall()
#convert bytes to string, the database was returning bytearrays for unique name and product value
result = [(row[0].decode(), row[1].decode(), *row[2:]) for row in result]


print("Content-type: text/html\n\n")
print(template.render(result=result))

conn.close()

