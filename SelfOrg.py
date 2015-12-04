#!/usr/bin/python

import random as random
import numpy as np
import math
import matplotlib.pyplot as plt


class Cell:
    """Store the information about species in cell"""

    def __init__(self, Level="", State="" ):
        """Initialize a new cell object

        Level -- identifier of the species
        State -- indicating wheither alive (True) or dead (False)
        """
        self.Level = Level
        self.State = State

    # this makes that you can do 'print cell':
    def __str__(self):
        """Return string representation of a cell object"""
        return "Level %s (Alive: %s)" % (str(self.Level), str(self.State))
    def __float__(self):
        return float(self.Level)

class CAmodel:
    def __init__(self, M,n,steps):
        self.steps = steps
        self.n = n
        self.M = M

        self.cmplxt = []

        self.Grid = self.initializeGrid(M,n)
        self.Grid = self.reorderGrid(self.Grid)


    # done
    def initializeGrid(self,M,n):
        """Grid initizialed by funtions startLevel and startState
        Function takes the number of species M and dims of the
        2-D lattice (n by n)"""
        Grid = []
        for i in range(n):
            row = []
            for j in range(n):
                row.append(Cell(self.startLevel(M),self.startState()))
            Grid.append(row)
        return Grid

    def reorderGrid(self, Grid):
        numbers = []
        for i in range(self.M):
            numbers.append([i,0])
        for row in Grid:
            for cell in row:
                numbers[cell.Level-1][1] += 1
        numbers = sorted(numbers,key=lambda x:x[1],reverse=True)

        new_Grid = []
        for row in Grid:
            new_row = []
            for cell in row:
                for i in range(len(numbers)):
                    if numbers[i][0] == cell.Level-1:
                        new_row.append(Cell(i+1, cell.State))
                        break
                else:
                    print "error"
            new_Grid.append(new_row)
        return new_Grid


    # done, completely random atm
    def startLevel(self,M):
        """Returns a completely random level (int) to start at"""
        return random.randint(0,M-1)

    # done, completely random atm
    def startState(self):
        """Returns True or False indicating dead or alive, randomly"""
        return True#bool(random.getrandbits(1))

    # done
    def updateGrid(self):
        """calculate next step for every cell, and save in self.Grid"""
        self.cmplxt.append(self.getComplexity())
        new_grid = self.Grid
        for i in range(n):
            for j in range(n):
                newlevel,newstate = self.updateCell(i,j)
                new_grid[i][j] = Cell(newlevel,newstate)
        self.Grid = new_grid
        return

    # to do
    def updateCell(self,i,j):
        """Gets the coordinates and returns the new level and State"""
        num_pred,num_prey,num_alive = self.getNums(i,j)
        # case1
        if num_pred<num_prey and num_prey>0:
            return self.Grid[i][j].Level,True
        else:
            return self.Grid[i][j].Level,False

    # done
    def getNums(self,i,j):
        """Gets the coordinates and returns the number of preds and preys in neighbrhd"""
        nb = self.getNeighbors(i,j)
        act_level = self.Grid[i][j].Level
        num_prey = 0
        num_pred = 0
        num_alive = 0
        for k in nb:
            if k.Level == (act_level + 1) and k.State == True:
                num_pred += 1
            elif k.Level == (act_level -1) and k.State == True:
                num_prey += 1
            if k.State == True:
                num_alive += True
        return num_pred,num_prey,num_alive

    # done
    def getNeighbors(self,i,j):
        nb = []
        n = self.n
        for k in range(-1,2,1):
            for l in range(-1,2,1):
                if (k==0 and l==0):
                    continue
                else:
                    nb.append(self.Grid[(i+k)%n][(j+l)%n])
        return nb


    # working on it
    def getComplexity(self):
        p_i = np.zeros((self.M,1))
        n = self.n

        # loop over all cells to determine probabilities
        for i in range(n):
            for j in range(n):
                k = self.Grid[i][j]
                # if cell is alive, update corresp. probability
                if k.State == True:
                    p_i[k.Level] += 1./(n**2)

        # calculate entropy
        S = 0.
        for p in p_i:
            try:
                S += (p*math.log(p))[0]
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
        cbar = fig.colorbar(cax, ticks=[0,1,self.M/2.0,self.M])
        cbar.ax.set_yticklabels([0, self.M/2, self.M])  # vertically oriented colorbar
        plt.show()


    # done
    def run(self):
        for t in range(self.steps):
            self.updateGrid()
        print "model run completed"


M = 8 # number of interacting species
n = 100 # dimensions of (square) 2-D lattice
steps = 3 # number of steps

model = CAmodel(M,n,steps)
model.printMatrix()
model.run()
