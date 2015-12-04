#!/usr/bin/python

import random as random
import numpy as np

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

class CAmodel:
    def __init__(self, M,n,steps):
        self.Grid = self.initializeGrid(M,n)
        self.steps = steps
        self.n = n

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

    # done, completely random atm
    def startLevel(self,M):
        """Returns a completely random level (int) to start at"""
        return random.randint(1,M)

    # done, completely random atm
    def startState(self):
        """Returns True or False indicating dead or alive, randomly"""
        return bool(random.getrandbits(1))

    # done
    def updateGrid(self):
        """calculate next step for every cell, and save in self.Grid"""
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
        num_pred,num_prey, num_nb = self.getNums(i,j)
        # case 1 number of preys > predator
        if num_prey>num_pred:
            if self.Grid[i][j].State:
                # if the cell is alive, keep it alive
                return self.Grid[i][j].Level,True
            else:
                # if the cell is dead, return to life
                return self.Grid[i][j].Level,True

        #case 2 number of predators > prey
        elif num_pred>num_prey:
            #the cell dies
            return self.Grid[i][j].Level,False
        #case3, no predators or preys apply Conway's rules
        elif num_pred == 0 and  num_prey == 0:
            if self.Grid[i][j].State:
                # the cell was alive
                if num_nb == 2 or num_nb == 3:
                    # the cell stays alive when it has 2 or 3 living neighbors
                    return self.Grid[i][j].Level,True
                else:
                    # the cell dies otherwise
                    return self.Grid[i][j].Level,False
            else:
                # the cell was dead
                if num_nb == 3:
                    # return to life if it has 3 neighbors
                    return self.Grid[i][j].Level,True
                else:
                    # keep dead
                    return self.Grid[i][j].Level,False
        #case4, predators are in the same number as preys
        elif num_pred == num_prey:
            # do not change state
           return self.Grid[i][j].Level,self.Grid[i][j].State
    
        

    # to do
    def getNums(self,i,j):
        """Gets the coordinates and returns the number of preds and preys in neighbrhd"""
        nb = self.getNeighbors(i,j)
        return 1,2

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

    # done
    def run(self):
        for t in range(self.steps):
            self.updateGrid()
        print "model run completed"


M = 8 # number of interacting species
n = 100 # dimensions of (square) 2-D lattice
steps = 3 # number of steps

model = CAmodel(M,n,steps)
model.run()
# print model.Grid
