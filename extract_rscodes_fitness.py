import MySQLdb
import numpy as np
import pickle
import pandas as pd
import time
from sys import argv

filename, first = argv

print first


# Open sql database
db = MySQLdb.connect(host="localhost",
                     user="root",
                     passwd="r3pxto0fjk",
                     db="athgene")

# Prep database for queries
cur = db.cursor()


# Extract id from test_kits that are version 2
cur.execute("SELECT id FROM athgene.gene_markers;")
ID = [x[0] for x in cur.fetchall()]


cur.execute("SELECT rs_code FROM athgene.gene_markers;")
rs_code = [x[0] for x in cur.fetchall()]

cur.execute("SELECT result_identifier FROM athgene.gene_markers;")
result_identifier = [x[0] for x in cur.fetchall()]

print len(rs_code)
print len(result_identifier)

handle = open(first + "/" + "rs_code_fitness.txt", "w+")
for i in range(len(rs_code)):
    if result_identifier[i] == None:
        handle.write("%s\n" % rs_code[i])
    else:
        handle.write("%s\n" % result_identifier[i])
handle.close()
