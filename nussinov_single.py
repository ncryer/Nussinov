# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 08:59:27 2013

@author: nicolai
"""
import numpy as np

class Nussinov:
    """
    This class implements Nussinov's algorithm
    to predict RNA secondary structure from a
    RNA sequence alone. 
    
    It contains the following methods:
    
    returnVienna(seq,pairs):
    Return the Vienna encoding of a prediction made by
    the Nussinov algorithm
    Input: sequence used for prediction, list of predicted basepairs
    Output: A Vienna encoded string e.g. "..(.)..((..))"
    
    read_fasta("file.fasta"):
    Read a fasta file and store the sequence(s) in a dictionary
    of the form seqDict["seqName"] = "aauu..." 

    
    score_function(x,y):
    A base-pairing score function which positively weights
    Watson-Crick base pairs. Can be modified to accomodate 
    wobble base pairs, etc. 
    
    predict(seq):
    Fills out the score matrix in O(N^3) time and O(N^2) space
    Backtracks through the matrix in O(N) time 
    and outputs a prediction in Vienna format
    
    MatthewsCC(pred,true):
    Calculate and output the Matthews correlation coefficient
    for a prediction given its true/known structure
    """
    
    def __init__(self):
        return None
    
    def returnVienna(self,seq,pairs):
        """
        Return the structure in Vienna format
        
        Input: list of base pairs
               sequence of interest
        Output: A string of Vienna encoding
        """
        
        output = ["."] * len(seq)
        
        for a, b in pairs:
            output[a] = "("
            output[b] = ")"
        return "".join(output)
    
    
    def read_fasta(self,filename):
        """
        Comments
        TODO: What this function returns
        """
        
        with open(filename,"r") as myfile:
            tempfile = myfile.readlines()
        
        # Store sequences in a dictionary for later retrieval
        seqDict = {}
        
        # Find all lines containing names
        name = ""
        for line in tempfile:
            line = line.strip()
            if line.startswith(">"):
                name = line
            seqDict[name] = line
        return seqDict
        
    def score_fun(self,x,y):
        """
        Internal function for use in prediction
        x and y are nucleotides
        This function matches the nucleotides and returns a score as an int
        """
        x = str(x.upper())
        y = str(y.upper())
        
        scores = {
                    "A" : {"U" : 1},
                    "U" : {"A" : 1,
                           "C" : 1},
                    "C" : {"G" : 1},
                    "G" : {"U" : 1},
                                
        }
        
        try:
            return scores[x][y]
        except KeyError:
            return 0
        
                
    def predict(self,seq):
        """
        Input is a nucleotide sequence
        
        K is a constant specifying how many in-between nucleotides we 
        want between basepairs. Controls "size of loops"
        
        Output is optimal path through the DP matrix
        """
        
        K = 3
        N = len(seq)
        DPmat = np.zeros((N,N))

        for n in range(K,N):
            for j in range(n,N):
                i = j-n
                # Down
                down = DPmat[i+1,j]
                # Left
                left = DPmat[i,j-1]
                # Diag - match 
                diag = DPmat[i+1,j-1] + self.score_fun(seq[i],seq[j])
                # Search k 
                mymax = 0
                for k in range(i,j-1):
                    subs = DPmat[i,k] + DPmat[k+1,j]
                    if mymax < subs:
                        mymax = subs
                DPmat[i,j] = max(down,left,diag,mymax)

        
        # Build output - best path through matrix
        # starting at position (0,N)
        path = [ (0,N-1) ]
        
        pairs = []
        
        while len(path) > 0:
            i,j = path.pop()
            
            if j - i < K:
                continue
            


            # Survey neighborhood
            left = DPmat[i,j-1]
            down = DPmat[i+1,j]
            diag = DPmat[i+1,j-1]
            
            
            match = self.score_fun(seq[i],seq[j])
            

            if DPmat[i,j] == diag + match and match > 0:
                pairs.append( (i,j) )
                path.append( (i+1,j-1) )

                
            elif DPmat[i,j] == left:
                path.append( (i,j-1) )

            elif DPmat[i,j] == down:
                path.append( (i+1,j) )

            else:
                # Find k that maximizes the local search
                mymax = 0
                maxk = 0
                for k in range(i,j-1):
                    subs = DPmat[i,k] + DPmat[k+1,j]
                    if mymax < subs:
                        mymax = subs
                        maxk = k
                path.append( (i,maxk) )
                path.append( (maxk+1,j) )
                
        return pairs

    def MatthewsCC(self,true,pred):
        """
        Examine the prediction accuracy via the MCC
        """
        
        TP = 0; TN = 0; FP = 0; FN = 0
        
        for i in range(len(true)):
            if true[i] == pred[i] and true[i] in "()":
                TP += 1
            elif true[i] == pred[i] and true[i] == ".":
                TN += 1
            else:
               if pred[i] in "()":
                   FP += 1
               if pred[i] == ".":
                   FN += 1
        
        sensitivity = TP / float((TP+FN))
        specificity = TN / float(TN+FP)
        
        print "Sensitivity: " + str(sensitivity)
        print "Specificity: " + str(specificity)
        
        MCC = (TP*TN - FP*FN) / float( np.sqrt( (TP+FN)*(TN+FP)*(TP+FP)*(TN+FN) ) )
        
        return MCC



