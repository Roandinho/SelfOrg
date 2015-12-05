#!/usr/bin/pypy

import random as random
# import numpy as np
import math
#import matplotlib.pyplot as plt


class Cell:
    """Store the information about species in cell"""

    def __init__(self, Level=0, State=True):
        """Initialize a new cell object

        Level -- identifier of the species
        State -- indicating wheither alive (True) or dead (False)
        """
        self.Level = Level
        self.State = State
        self.nextState = State
        self.Grid = None
        self.neighbours = None
        self.preys = []
        self.preds = []
        self.rest = []

    def setNeighbours(self, cells):
        self.neighbours = cells
        for cell in cells:
            if cell.Level == (self.Level + 1):
                self.preds.append(cell)
            elif cell.Level == (self.Level - 1):
                self.preys.append(cell) # changed preds to preys, right??
            else:
                self.rest.append(cell)

    def doConway(self):
        n_alive = 0
        if self.State:
            for cell in self.rest:
                if cell.State:
                    n_alive += 1
                if n_alive > 3:
                    self.nextState = False
                    return
            if n_alive < 2:
                self.nextState = False
                return
        else:
            for cell in self.rest:
                if cell.State:
                    n_alive += 1
                if n_alive > 3:
                    self.nextState = False
                    return
            if n_alive < 3:
                self.nextState = False
                return
        self.nextState = True

    def updateStep(self, forcedExtinct):
        # if cell is forced to extinct, do this at first
        if self.State in forcedExtinct:
            self.nextState = False
            return

        # a large portion of the cells do not have pred or prey neighbours, do
        # conway immediately
        if len(self.preds) + len(self.preys) == 0:
            self.doConway()
            return

        # calculate balance pred/prey
        n_alive = 0
        coef = 0
        for prey in self.preys:
            if prey.State:
                n_alive += 1
                coef += 1
        for preds in self.preds:
            if preds.State:
                n_alive += 1
                coef -= 1
        # no preds or prey alive is conway
        if not n_alive:
            self.doConway()
            return
        elif coef == 0:
            self.nextState = self.State
        elif coef < 0:
            self.nextState = False
        else:
            self.nextState = True

    # this makes that you can do 'print cell':
    def __str__(self):
        """Return string representation of a cell object"""
        return "Level %s (Alive: %s)" % (str(self.Level), str(self.State))

    def __float__(self):
        return float(self.Level)


class CAmodel:
    Grid = None

    def __init__(self, M, n, steps, extProbs):
        # store variables
        self.steps = steps
        self.n = n
        self.M = M
        self.extProbs = extProbs

        # variables that are needed for statistics
        self.cmplxt = []
        self.numExtinct = []

        # create grid
        CAmodel.Grid = self.initializeGrid(M, n)
        CAmodel.Grid = self.reorderGrid(self.Grid)
        self.setCellNeighbours()

    def initializeGrid(self, M, n):
        """Grid initizialed by funtions startLevel and startState
        Function takes the number of species M and dims of the
        2-D lattice (n by n)"""
        Grid = []
        for i in range(n):
            row = []
            for j in range(n):
                row.append(Cell(self.startLevel(M), self.startState()))
            Grid.append(row)
        return Grid

    def reorderGrid(self, Grid):
        numbers = []
        for i in range(self.M):
            numbers.append([i, 0])
        for row in Grid:
            for cell in row:
                numbers[cell.Level][1] += 1
        numbers = sorted(numbers, key=lambda x: x[1], reverse=True)

        new_Grid = []
        for row in Grid:
            new_row = []
            for cell in row:
                for i in range(len(numbers)):
                    if numbers[i][0] == cell.Level:
                        new_row.append(Cell(i, cell.State))
                        break
                else:
                    print ("error")
            new_Grid.append(new_row)
        return new_Grid

    def setCellNeighbours(self):
        for i in range(n):
            for j in range(n):
                CAmodel.Grid[i][j].setNeighbours(
                    [CAmodel.Grid[(i + 1) % n][(j - 1) % n],
                     CAmodel.Grid[(i + 1) % n][(j) % n],
                     CAmodel.Grid[(i + 1) % n][(j + 1) % n],
                     CAmodel.Grid[(i) % n][(j - 1) % n],
                     CAmodel.Grid[(i) % n][(j + 1) % n],
                     CAmodel.Grid[(i - 1) % n][(j - 1) % n],
                     CAmodel.Grid[(i - 1) % n][(j) % n],
                     CAmodel.Grid[(i - 1) % n][(j + 1) % n]])

    # done, completely random atm
    def startLevel(self, M):
        """Returns a completely random level (int) to start at"""
        return random.randint(0, M-1)

    # done, completely random atm
    def startState(self):
        """Returns True or False indicating dead or alive, randomly"""
        return bool(random.randint(0, 2))

    # done
    def updateGrid(self):

        # keep track of those species who are alive at this timestep
        alive = dict()

        # calculate if species is forced to extinct
        forcedExtinct = self.calcForced()
        for row in CAmodel.Grid:
            for cell in row:
                cell.updateStep(forcedExtinct)
        for row in CAmodel.Grid:
            for cell in row:
                cell.State = cell.nextState
                # check for extinction
                if cell.State:
                    alive[cell.Level] = 1

        # update statistics
        self.cmplxt.append(self.getComplexity())
        self.numExtinct.append(self.M - len(alive))


    # done
    def calcForced(self):
        """calculates (at each timestep) which species are forced to die out
        the function is written in such a way, that extinction is not
        necessarily irreversible"""
        forcedExtinct = []
        probs = self.extProbs
        for p in probs:
            if probs[p] > random.random():
                forcedExtinct.append(p)
        return forcedExtinct

    # done
    def getComplexity(self):
        p_i = [0]*self.M
        # loop over all cells to determine probabilities
        for row in CAmodel.Grid:
            for cell in row:
                if cell.State:
                    # if cell is alive, update corresp. probability
                    p_i[cell.Level] += 1
        p_i = [p_ii/float(self.n**2) for p_ii in p_i]

        # calculate entropy
        S = 0.
        for p in p_i:
            try:
                S += (p*math.log(p))
            except ValueError:
                pass
        S = -1.*S
        return S

    def printMatrix(self):
        new_grid = []
        for row in self.Grid:
            new_row = []
            for cell in row:
                new_row.append(float(cell))
            new_grid.append(new_row)
        fig, ax = plt.subplots()
        cax = ax.imshow(new_grid)
        ax.set_title('Randomly Initialized grid')
        cbar = fig.colorbar(cax, ticks=[0, 1, self.M/2.0, self.M])
        cbar.ax.set_yticklabels([0, self.M/2, self.M])
        plt.show()

    def printEntireMatrix(self):
        """prints Level plus +(Alive) or -(Dead) for each cell"""
        TheString = ""
        for row in self.Grid:
            for cell in row:
                TheString += str(cell.Level).zfill(3)
                if cell.State:
                    TheString += "+ "
                else:
                    TheString += "- "
            TheString += "\n\n"
        print TheString
        return

    def printEntropy(self):
        y = []
        for i in range(1, self.steps):
            y.append(i)
        plt.plot(y, model.cmplxt[1:])
        plt.show()

    # done
    def run(self):
        for t in range(self.steps):
            self.updateGrid()
        print ("model run completed")


M = 256  # number of interacting species
n = 100  # dimensions of (square) 2-D lattice
steps = 1000  # number of steps
extProbs = dict() # dictionary containing 
extProbs[16] = 0.01

model = CAmodel(M, n, steps, extProbs)
# model.printMatrix()
model.run()
# import cProfile
#cProfile.run('model.run()', sort='tottime')
#model.printEntropy()
print model.cmplxt[1:-1:50]
print model.numExtinct[1:-1:50]
model.printEntireMatrix()