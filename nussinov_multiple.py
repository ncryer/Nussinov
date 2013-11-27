# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 20:57:44 2013

@author: nicolai
"""
import numpy as np

class Nussinov_multiple:
    """
    This class implements the Nussinov algorithm modified
    to gain additional prediction accuracy by interpreting 
    a multiple alignment of sequences, such that we may
    make use of compensatory mutations at the base-pairing loci
    by calculating the Mutual Information between columns in the alignment
    
    Methods:
    
    read_fasta(filename):
    Reads a multiple alignment 
    
    MI(x,y):
    calculate the mutual information between positions x and y
    
    
    """
    def __init__(self):
        # Store the multiple alignment as a list of lists
        self.alignment = None
        
        
    
    def read_alignment(self,filename):
        """
        Read in a multiple alignment in fasta format
        store it sensibly for access in self.alignment
        """
        with open(filename) as myfile:
            myfile = myfile.readlines()
        
        myfile = [line.strip() for line in myfile]
        seqs = []
        for line in myfile:
            if line.startswith(">"):
                continue
            else:
                seqs.append(line)
        
        self.alignment = [ list(seq) for seq in seqs ]
        
    
    def MI(self,i,j):
        """
        i,j refer to the ith and jth columns of a multiple alignment
        
        This method extracts cols i and j from self.alignment and
        calculates the mutual information
        """
        L = len(self.alignment)
        
        col_i = [ self.alignment[x][i] for x in range(L) ]
        col_j = [ self.alignment[x][j] for x in range(L) ]

        
        # Calculate single frequencies in each column
        prob_i = {}
        prob_j = {}
        
        for nuc in set(col_i):
            prob_i[nuc] = 0
        for nuc in set(col_j):
            prob_j[nuc] = 0
        
        for x in range(len(col_i)):
            prob_i[ col_i[x] ] += 1
            prob_j[ col_j[x] ] += 1

        # Normalization step
        norm = sum(prob_i.values())
        for key in prob_i.keys():
            prob_i[key] = prob_i[key]/float(norm)

        norm = sum(prob_j.values())
        for key in prob_j.keys():
            prob_j[key] = prob_j[key]/float(norm)

        # Calculate joint frequencies
        pairs = [ col_i[x] + col_j[x] for x in range(L) ]
        
        pairdict = {}
        
        for pair in set(pairs):
            pairdict[pair] = 0
        for pair in pairs:
            pairdict[pair] += 1
        
        
        # Normalization step
        for key in pairdict.keys():
            pairdict[key] = pairdict[key]/float(len(pairs))
        
        # Calculate Mutual Information
        terms = 0
        for x in range(len(col_i)):
            pair = col_i[x] + col_j[x]
            pairfreq = pairdict[pair]

            terms += pairfreq * np.log2(pairfreq/(prob_i[col_i[x]] * prob_j[col_j[x]]))
        
        return terms
        
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
        
                
    def predict(self):
        """
        Input is a nucleotide sequence
        
        K is a constant specifying how many in-between nucleotides we 
        want between basepairs. Controls "size of loops"
        
        Output is optimal path through the DP matrix
        """
        seq = "".join(self.alignment[1])
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
                diag = DPmat[i+1,j-1] + self.score_fun(seq[i],seq[j]) + self.MI(i,j)
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
            
            
            match = self.score_fun(seq[i],seq[j]) + self.MI(i,j)
            

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
        
    def returnVienna(self):
        """
        This method invokes a full prediction via
        the methods in this class
        
        
        Output: A string of Vienna encoding
        """
        
        output = ["."] * len(self.alignment[1])
        
        pairs = self.predict()
        
        for a, b in pairs:
            output[a] = "("
            output[b] = ")"
        return "".join(output)
    
    
hej = Nussinov_multiple()
hej.read_alignment("testseqs_aligned.fasta")
eh = hej.returnVienna()

