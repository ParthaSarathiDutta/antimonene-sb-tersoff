
import sys
import os
sys.path.append("/global/cfs/cdirs/m1917/blast_ff/troy_mcts_test")
from math import exp, sqrt, log, fabs
from random import random, choice, seed, shuffle
from time import time
from datetime import datetime
from scipy.optimize import minimize
import numpy as np
import csv

from MCTree import Tree, randomlist
from ParameterObject import ParameterData
from main1 import objective
from Node import Node
seed(datetime.now())

restartsim = True

import blast, ast

blast.logger.setLevel(30)
#model = sum(blast.model.load_from_json('model.json'))
model_dict = blast.json.read("model.json")
model = blast.model.load_from_dict(model_dict)
ranges = [ast.literal_eval(m['range'][p])
          for m in model.added for p in m['in']]
lbounds = [float(l[0]) for l in ranges]
ubounds = [float(l[1]) for l in ranges]
nParameters = len(ubounds)

if not restartsim:
    startset = [ (ux-lx)*random() + lx for lx,ux in zip(lbounds, ubounds)]
else:
    startset = [m.p[p] for m in model.read('mcts_restart.tersoff') for p in m['in']]
    #startset = restart()
print(lbounds)
print(ubounds)

# --------------------- Check for Parameter Boundaries ---------------------- #
def hitting_boundaries(startset):
    ranges = [ast.literal_eval(m['range'][p]) for m in model.added for p in m['in']]
    lbounds = np.array([l[0] for l in ranges],dtype=float)
    ubounds = np.array([l[1] for l in ranges],dtype=float)
    startset = np.array(startset,dtype=float)
    
    # adjust lbounds and ubounds to include startset
    
    ix=(lbounds>startset); lbounds[ix] = 0.9*startset[ix]
    ix=(ubounds<startset); ubounds[ix] = 1.1*startset[ix]

    # write to model
    i=0
    for m in model:
        for p in m['in']:
            m['range'][p] = f"[{lbounds[i]}, {ubounds[i]}]"
            i += 1
    
    model.save('model.json')
    return

#hitting_boundaries(startset)
# ------------------------------------------------------------------------- #

indata = ParameterData(parameters=startset, ubounds=ubounds, lbounds=lbounds)
objfunc = objective
indata.setevaluator(objfunc)
# print("PRINTING BELOW:")
# print(indata)
# print(objfunc)
def UBEnergy(nodelist, exploreconst):
    #Compute the average uniqueness factor
    uniqscore = [0.0 for x in nodelist]
    print('inside UBEnergy: START OF GETTING UNIQUENESS FACTOR')
   
    for i, node in enumerate(nodelist):
        print(f"inside UBEnergy: ITERATING NODE {i} FOR UNIQUENESSFACTOR")
        score = node.getuniquenessdata(nodelist)
        print(f"inside UBEnergy - score =  {score}")
        uniqscore[i] = score
    print(f"inside UBEnergy - uniqscore = {uniqscore}")
    
    maxeng = max([node.getscore() for node in nodelist])
    print(f"inside UBEnergy - eng: {[node.getscore() for node in nodelist]}")
    print(f"inside UBEnergy - maxeng: {[node.getscore() for node in nodelist]}")
    print('inside UBEnergy: END OF GETTING UNIQUENESS FACTOR')
    def UCT_Unique_Score(node, uniqval, doprint=False):
        parent = node.getparent()
        print(f"inside UCT_Unique_Score - parent = {parent}")
        energy = node.getscore()
        print(f"inside UCT_Unique_Score - energy = {energy}")
        if parent is None:
            return -1e30
        else:
            parenergy = parent.getscore()
            print(f"inside UCT_Unique_Score - parenergy = {parenergy}")
            parvisits = parent.getvisits()
            print(f"inside UCT_Unique_Score - parvisits = {parvisits}")
        
        visits = node.getvisits()
        print(f"inside UCT_Unique_Score - visits = {visits}")
        depth = node.getdepth()
        print(f"inside UCT_Unique_Score - depth = {depth}")
        playoutEList = node.getenergylist()
        print(f"inside UCT_Unique_Score - playoutEList = {playoutEList}")
        if doprint:
            print(str(node), len(playoutEList))
        nodeEnergy = node.getscore()
        print(f"inside UCT_Unique_Score - nodeEnergy = {nodeEnergy}")
        nodeweight = nodeEnergy
        print(f"inside UCT_Unique_Score - nodeweight = {nodeweight}")
        avgweight = nodeweight
        print(f"inside UCT_Unique_Score - avgweight = {avgweight}")
        if len(playoutEList) > 0:
            for energy in playoutEList:
                #print(f"inside UCT_Unique_Score - energy = {energy}")
                #print(f"inside UCT_Unique_Score - log(energy) = {log(energy)}")
#                avgweight += log(energy)
                avgweight = min(avgweight, log(energy))
            print(f"inside UCT_Unique_Score - avgweight = {avgweight}")
#            avgweight = avgweight/float(len(playoutEList)+1)
 
        nChildren = len(node.getchildren())
        print(f"inside UCT_Unique_Score - nChildren = {nChildren}")
        
        explore = -1
        if depth > 5:
            score = -1e20
        elif nChildren > 35/depth:
            score = -1e20
        else:
            try:
                explore = exploreconst*uniqval*sqrt(log(parvisits)/(visits*nParameters))
                print(f"inside UCT_Unique_Score - nParameters = {nParameters}")
                score = -avgweight + explore
            except (ValueError, ZeroDivisionError):
                explore = exploreconst*uniqval/sqrt(nParameters)
                score = -avgweight + explore
        if doprint:
            print("Node %s -- Depth:%s Exploit:%s Explore:%s Total-Score:%s "%(node.getid(), depth, -avgweight, explore, score))
        #quit()
        return score
    
    keylist = {}
    for i, node in enumerate(nodelist):
        keylist[str(node)] = uniqscore[i]
        print(type(node))
    print(f"inside UCT_Unique_Score - keylist = {keylist}")
    
    selection = sorted(nodelist, key=lambda x:x.getid())
    print(f"inside UCT_Unique_Score - selection = {selection}")
    selection = sorted(selection, key=lambda x:UCT_Unique_Score(x, keylist[str(x)], doprint=True))[-1]
    print(f"inside UCT_Unique_Score - selection = {selection}")
    print("Selecting Node %s with Score: %s"%(selection.getid(),  UCT_Unique_Score(selection, keylist[str(selection)], doprint=False)))
    #quit()
    return selection


#---Tree Class Test Code---
tree = Tree(seeddata=indata, playouts=5, selectfunction=UBEnergy, headexpansion=5)
#tree.loadtree('mctree.restart', seeddata=indata)
tree.expand(nExpansions=1, writeevery=5)
tree.setconstant(80)
for iLoop in range(1, 100000):
    print("Loop Number %s"%(iLoop))
    print('SIMULATE START')
    tree.simulate(nSimulations=5)
    print('SIMULATE END')
    print('PLAYEXPAND START')
    tree.playexpand(nExpansions=5)
    print('PLAYEXPAND END')
    print('EXPAND START')
    tree.expand(nExpansions=4)
    print('EXPAND END')
    if iLoop%15 == 0:
        tree.selectpath()
#    if iLoop%5 == 0:
    tree.savetree("mctree.restart")
    #quit()


tree.savetree("finaltree.restart")
del tree

print("Finished MCTS Loop")
