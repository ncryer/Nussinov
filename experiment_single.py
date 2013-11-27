# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 20:17:41 2013

@author: nicolai
"""

import nussinov_single

testseqs = Nussinov()

seqs = testseqs.read_fasta("testseqs.fasta")

preds = []
pred_dict = {}

for seqname in seqs.keys():
    
    seq = seqs[seqname]
    
    predpairs = testseqs.predict(seq)
    
    pred = testseqs.returnVienna(seq,predpairs)
    
    preds.append( (seqname,seq,pred) )
    pred_dict[seqname] = pred

# output experiment file in the correct order
output = open("single_pred.dbn","w")

with open("testseqs.fasta","r") as myfile:
    myfile = myfile.readlines()
myfile = [line.strip() for line in myfile]
for line in myfile:
    if line.startswith(">"):
        output.write( line + "\n" )
        output.write( pred_dict[line] + "\n" )
output.close()
    
    
    
    