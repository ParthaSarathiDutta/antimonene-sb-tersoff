from random import choice, random, shuffle, seed, randint
from math import fabs, cos, pi, sqrt
from time import time
from copy import deepcopy
from datetime import datetime
#from OptFunctions import HolderTable
import sys
from hyperopt import hp, fmin, tpe, space_eval, Trials
import numpy as np

seed(datetime.now())

searchdx = 1.0
depthscale = [0.7, 0.5, 0.3, 0.2, 0.1, 0.05]

#----------------------------------------------
def loggen(x, upper, lower):
    logupper = log(upper/x)
    loglower = log(lower/x)
    dV = (upper-lower)*random()+lower
    x_new = x * exp(dV)
#----------------------------------------------
def lineargen(x, upper, lower):
    xnew = (upper-lower)*random() + lower
    return xnew
#----------------------------------------------
#=================================================
def maxlist(inlist, weightlist):
    maxweight = -1e50
    for i, weight in enumerate(weightlist):
        if weight > maxweight:
            maxweight = weight
            indx = i
    return indx, inlist[indx]
#================================================
class MissingParameters(Exception):
    pass
#================================================
class InvalidParameterBounds(Exception):
    pass
#================================================
class ParameterData(object):
    #------------------------------------------------
    def __init__(self, parameters=[], evaluator=None, lbounds=None, ubounds=None):
        self.parameters = parameters
        if len(self.parameters) < 1:
            raise MissingParameters("Please specifiy an initial parameter set")
        self.dimensions = len(self.parameters)
        self.energy = 0.0
        self.ntrials = 20
        self.lbounds = lbounds
        self.ubounds = ubounds

        # This portion checks the consistency of the number of bounds and the parameters
        # to make sure the user properly defines their optimization problem.
        if len(self.lbounds) != len(self.ubounds):
            raise InvalidParameterBounds("The number of upper and lower bounds do not match!")
        if len(self.lbounds) != self.dimensions:
            raise InvalidParameterBounds("The number of bounds given do not match the number of parameters!")
        for lower, upper in zip(self.lbounds, self.ubounds):
            if upper <= lower:
                print(lower, upper)
                raise InvalidParameterBounds("Lower Bound >= Upper Bound!")


        self.bay_trials = None
        
        if evaluator is None:
            self.evaluator = None
        else:
            self.evaluator = evaluator

    #----------------------------------------------------
    def __eq__(self, ising2):
        str2 = ising2.getparameters()
        return str2 == self.parameters
    #----------------------------------------------------
    def __str__(self):
        outstr = ' '.join([str(x) for x in self.parameters])
        return outstr

    #------------------------------------------------
    def newdataobject(self):
        newobj = ParameterData(parameters=self.parameters, lbounds=self.lbounds, ubounds=self.ubounds, evaluator=self.evaluator)
        return newobj

    #------------------------------------------------
    def perturbate(self, node=None, parOnly=False):
        weights = []
        structs = []

        for iTrial in range(self.ntrials):
            newlist = None
            while newlist is None:
                newlist = self.localshift(node=node)
                print('I AM perturbate')
                print(newlist)
            if newlist is None:
                print("Unable to perform move!")
                continue
            break
        selection = newlist
        print(f'IN PeRTURBATE, I SELECTED - {selection}')
        if not parOnly:
            newobj = ParameterData(parameters=selection, lbounds=self.lbounds, ubounds=self.ubounds, evaluator=self.evaluator)
            return newobj
        else:
            return selection
    #------------------------------------------------
    def localshift(self, node=None):
        if node is not None:
            depth = node.getdepth()
            if depth < len(depthscale):
                searchmax = depthscale[depth]
                ranfunc = lineargen
            else:
                searchmax = depthscale[-1]
                ranfunc = lineargen
        else:
            depth = 1
            searchmax = 1.0
            ranfunc = lineargen



        newlist = deepcopy(self.parameters)
        print('Hi! I AM LoCaLShIfT')
        print(f'newlist - {newlist}')
        print(f'self.parameters - {self.parameters}')
        for i, par in enumerate(self.parameters):
            upper = par + searchmax*(self.ubounds[i]-self.lbounds[i])
            if upper > self.ubounds[i]:
                upper = self.ubounds[i]
            lower = par - searchmax*(self.ubounds[i]-self.lbounds[i])
            if lower < self.lbounds[i]:
                lower = self.lbounds[i]
            xnew = (upper-lower)*random() + lower
            newlist[i] = xnew
        return newlist
    #----------------------------------------------------
    def runsim(self, playouts=1, moves=3, mass=True, node=None):
        tries = 0
        structlist = []
        energylist = []
        if node is not None:
            print("Node Depth: %s"%(node.getdepth()))
        for playout in range(playouts):
            newlist = self.parameters
            for iMove in range(moves):
                tempstr = None
                while tempstr is None:
                    
                    tempstr = self.perturbate(parOnly=True, node=node)
                    print(f'IN RUNSIM - {tempstr}')
                newlist = tempstr
            structlist.append(newlist)
            energy = self.evaluator(newlist)
            print(' I was here too')
            print("Playout %s Result: %s"%(playout, energy))
            energylist.append(energy)
        return energylist, structlist
     #----------------------------------------------------
    def runsim_bayesian(self, playouts=1, moves=3, mass=True, node=None):
        tries = 0
        structlist = []
        energylist = []

        if node is not None:
            depth = node.getdepth()
            if depth < len(depthscale):
                searchmax = depthscale[depth]
            else:
                searchmax = depthscale[-1]
            print("Node Depth: %s"%(node.getdepth()))

        space = []
        for i, par in enumerate(self.parameters):
            upper = par + searchmax*(self.ubounds[i]-self.lbounds[i])
            if upper > self.ubounds[i]:
                upper = self.ubounds[i]
            lower = par - searchmax*(self.ubounds[i]-self.lbounds[i])
            if lower < self.lbounds[i]:
                lower = self.lbounds[i]
            space.append(hp.uniform(str(i), lower, upper))

        if self.bay_trials is None:
            trials = Trials()
            nprev = 0
        else:
            trials = self.bay_trials
            nprev = len(trials.trials)


        try:
            results = fmin(self.evaluator, space=space, algo=tpe.suggest, max_evals=nprev+playouts, trials=trials, verbose=0)
        except:
            pass
        self.bay_trials = trials
        cnt = nprev
        print(len(trials.trials))
        print(len(trials.results))

        for trial, result in zip(trials.trials[nprev:], trials.results[nprev:]):
            newlist = [trial['misc']['vals'][str(i)][0] for i in range(len(self.parameters))]
            energy = result['loss']
            cnt += 1
            print("Playout %s Result: %s"%(cnt, energy))
            structlist.append(newlist)
            energylist.append(energy)
        newlist = [results[str(x)] for x in range(self.dimensions)]
        energy = min(trials.losses())
        return energylist, structlist     
    #----------------------------------------------------
    def computescore(self, node=None):
        self.energy = self.evaluator(self.parameters)
        return self.energy

    #----------------------------------------------------
    def getuniqueness(self, inlist=None, node=None, nodelist=None):
        if node is None:
            raise
        depth = node.getdepth()
        #print(f"Inside PO - depth: {depth}")
        if depth == 0:
            #print(f"Inside PO - depth: 0")
            return 0.0
        siblist = node.parent.getchildren()
        #print(f"Inside PO - siblist - {siblist}")
        if inlist is None:
            comparelist = self.parameters
            #print(f"Inside PO - comparelist - {comparelist}")
        else:
            comparelist = inlist
            #print(f"Inside PO - comparelist - {comparelist}")
        cnt = 0.0
        score = 0.0
        for othernode in siblist:
            if node.getid() == othernode.getid():
                continue
            libstruct = othernode.getdata().getstructure()
            #print(f"Inside PO - libstruct - {libstruct}")
            rsq = 0.0
            for i, x in enumerate(comparelist):
                rxi = comparelist[i]/(self.ubounds[i]-self.lbounds[i])
                rxj = libstruct[i]/(self.ubounds[i]-self.lbounds[i])
                rsq += (rxi-rxj)**2
            r = sqrt(rsq)
            #print(f"Inside PO - r = {r}")
            if depth > len(depthscale)-1:
                scale = depthscale[-1]
                print(f"Inside PO - depthscale = {depthscale}")
            else:
                scale = depthscale[depth]
            if r < 1e-5/scale:
                score -= 1.0
            score += r
            cnt += 1.0
        if cnt > 0.0:
            score = score/cnt
        else:
            score = 0.03
        return score


    #----------------------------------------------------
    def getlibuniqueness(self, inlist=None, node=None):
        if inlist is None:
            curlist = self.parameters
        else:
            curlist = inlist

        if node is not None:
            depth = node.getdepth()
        else:
            depth = 1
            
        cnt = 0.0
        score = 0.0
        with open("dumpfile.dat", "r") as infile:
            for line in infile:
                libstruct = line.split("|")[0].split()
                libstruct = [float(x) for x in libstruct]
                rsq = 0.0
                for pari, parj in zip(curlist, libstruct):
                    rsq += (pari-parj)**2
                r = sqrt(rsq)
                if r < 1e-10:
                    score -= 3.0
                score += r*sqrt(depth)
                cnt += 1.0
        score = score/cnt
        return score
    #----------------------------------------------------
    def getstructure(self):
        return self.parameters
    #----------------------------------------------------
    def setstructure(self, inlist):
        self.parameters = inlist
    #----------------------------------------------------
    def setevaluator(self, evalfunc):
        self.evaluator = evalfunc

    #----------------------------------------------------
    def convertstr(self, instr):
        par = [float(x) for x in instr.split()]
        return par

    #----------------------------------------------------
#================================================


def Rastrigin(parameters):   
    f = 10*len(parameters)
    for value in parameters:
        f += value**2 - 10*cos(2.0*pi*value)
    with open("dumpfile.dat", "a") as outfile:
        outstr = ' '.join([str(x) for x in parameters])
        outfile.write('%s | %s \n'%(outstr, f))
    return f


#================================================
def Rastrigin_Print(parameters):   
    f = 10*len(parameters)
    for value in parameters:
        f += value**2 - 10*cos(2.0*pi*value)
    with open("dumpfile.dat", "a") as outfile:
        outstr = ' '.join([str(x) for x in parameters])
        outfile.write('%s | %s \n'%(outstr, f))
    print("%s %s"%(parameters, f))
    return f
#================================================

