#!/usr/bin/pypy

import random as random
import math
from CAmodel import CAmodel
M =256  # number of interacting species
n = 100  # dimensions of (square) 2-D lattice
steps = 2 # number of steps
extProbs = dict()  # dictionary containing extinction probablities
# extProbs[16] = 0.01
# extProbs[90] = 0.001
model = CAmodel.CAmodel(M, n, steps, perturbimpact=10, hierarchy=True,
        replaceExtinct=False, hierarchyRange=5,
        neighbourhood_level =1, perturb=True, singleRuns=True)
# model.printMatrix()
model.run()
# import cProfile
# cProfile.run('model.run()', sort='tottime')
# print model.cmplxt[1:-1:50]
# print model.percentageAlive[1:-1:50]
# print model.numExtinct
model.printEntropy()
#print model.avalanceLengths
#model.printCSV(avalanche=True, entropy=True)
#print model.single_run_avalance
#print model.extinctions
# print model.dead
# model.printEntireMatrix()
