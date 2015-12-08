#!/usr/bin/pypy

import random as random
import math
plt = None


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
                self.preys.append(cell)  # changed preds to preys, right??
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
        if self.Level in forcedExtinct:
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

    # nice way to store previous states
    def __repr__(self):
        return "%s%s" % (str(self.Level), str(int(self.State)))

    # this makes that you can do 'print cell':
    def __str__(self):
        """Return string representation of a cell object"""
        return "Level %s (Alive: %s)" % (str(self.Level), str(self.State))

    def printNeighbors(self):
        print ("For cell of level", self.Level, "and alive:", self.State)
        print ("PREYS:")
        for prey in self.preys:
            print (prey)
        print ("PREDS:")
        for pred in self.preds:
            print (pred)
        print ("Rest:")
        for rest in self.rest:
            print (rest)

    def __float__(self):
        return float(self.Level)


class Grid(list):
    def __repr__(self):
        ret = ""
        for row in self:
            for cell in row:
                ret += repr(cell)
        return ret

    def __hash__(self):
        return hash(self.__repr__())

    def __str__(self):
        return self.__repr__()


class ConvergedException(Exception):
    pass


class CAmodel:
    Grid = None
    withExtinction = False

    def __init__(self, M, n, steps, extProbs=None, perturbimpact=1):
        # store variables
        self.calcConvergence = False
        self.previousStates = dict()
        self.perturbimpact = perturbimpact
        self.avalanceLengths = []
        self.steps = steps
        self.n = n
        self.M = M
        if extProbs:
            self.extProbs = extProbs
            self.withExtinction = True

        # variables that are needed for statistics
        self.cmplxt = []  # the complexity for all timestep
        # the number of species that are extinct for multiple timesteps
        self.numExtinct = []
        # the percentage of all cell's which are alive for every timestep
        self.percentageAlive = []
        # the number of species alive for every timestep
        self.nAlive = []

        # create grid
        CAmodel.Grid = self.initializeGrid(M, n)
        CAmodel.Grid = self.reorderGrid(self.Grid)
        self.setCellNeighbours()

    def initializeGrid(self, M, n):
        """Grid initizialed by funtions startLevel and startState
        Function takes the number of species M and dims of the
        2-D lattice (n by n)"""
        grid = Grid()
        for i in range(n):
            row = []
            for j in range(n):
                row.append(Cell(self.startLevel(M), self.startState()))
            grid.append(row)
        return grid

    def reorderGrid(self, grid):
        numbers = []
        for i in range(self.M):
            numbers.append([i, 0])
        for row in grid:
            for cell in row:
                numbers[cell.Level][1] += 1
        numbers = sorted(numbers, key=lambda x: x[1], reverse=True)

        new_Grid = Grid()
        for row in grid:
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

    def startLevel(self, M):
        """Returns a completely random level (int) to start at"""
        return random.randint(0, M-1)

    def startState(self):
        """Returns True or False indicating dead or alive, randomly"""
        return bool(random.randint(0, 2))

    def updateGrid(self, timestep):

        # keep track of those species who are alive at this timestep
        alive = [0]*self.M

        # calculate if species is forced to extinct
        forcedExtinct = []
        if self.withExtinction:
            forcedExtinct = self.calcForced()

        # percentage alive
        percAlive = 0

        # UPDATE, i.e.: calc next grid and then copy into grid
        for row in CAmodel.Grid:
            for cell in row:
                cell.updateStep(forcedExtinct)
                if cell.State:
                    percAlive += 1

        for row in CAmodel.Grid:
            for cell in row:
                cell.State = cell.nextState

                # check for extinction
                if cell.State:
                    alive[cell.Level] = 1

        #todo, append after update statistics?
        if self.calcConvergence:
            st = str(CAmodel.Grid)
            try:
                self.previousStates[st]
            except:
                self.previousStates[st] = 1
            else:
                raise ConvergedException()
        # update statistics
        self.cmplxt.append(self.getComplexity())
        self.percentageAlive.append((1.*percAlive)/self.n**2)
        self.nAlive.append(sum(alive))

        # if timestep % (self.steps/10) == 0:
        #    totalExtinct = sum(self.numExtinct)
        #    print (timestep, (self.steps/10), self.steps % (self.steps/10),
        #           self.M, len(alive), (self.M - len(alive)))
        #    self.numExtinct.append((self.M - len(alive))-totalExtinct)

    def calcForced(self):
        """calculates (at each timestep) which species are forced to die out
        the function is written in such a way, that extinction is not
        necessarily irreversible"""
        forcedExtinct = []
        probs = self.extProbs
        for p in probs:
            p_random = random.random()
            if probs[p] > p_random:
                forcedExtinct.append(p)
        return forcedExtinct

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
        import subprocess
        with open("entropy.py", "w") as f:
            f.write("import matplotlib.pyplot as plt\n" +
                    "y = []\n" +
                    "for i in range(1, " + str(len(self.cmplxt)) + "):\n " +
                    "  y.append(i)\n" +
                    "plt.plot(y, " + str(self.cmplxt[1:]) + ")\n" +
                    "plt.show()")
            subprocess.Popen(["python2", "entropy.py"])

    def perturb(self):
        for _ in range(self.perturbimpact):
            i = random.randint(0, self.n-1)
            j = random.randint(0, self.n-1)
            #toggle/perturb
            self.Grid[i][j].State = not self.Grid[i][j].State


    def run(self, perturb=False):
        if perturb:
            self.calcConvergence=True
            i = 0
            j = 0
            while True:
                try:
                    self.updateGrid(j)
                except ConvergedException:
                    self.perturb()
                    self.previousStates=dict()
                    i += 1
                    self.avalanceLengths.append(j)
                    j = 0
                    continue
                if i == self.steps:
                    break
                j += 1
        else:
            self.calcConvergence=False
            for t in range(self.steps):
                self.updateGrid(t)
            print ("model run completed")


M = 256  # number of interacting species
n = 100  # dimensions of (square) 2-D lattice
steps = 100  # number of steps
extProbs = dict()  # dictionary containing extinction probablities
# extProbs[16] = 0.01
# extProbs[90] = 0.001

model = CAmodel(M, n, steps, extProbs=None, perturbimpact=10)
# model.printMatrix()
model.run(perturb=True)
# import cProfile
# cProfile.run('model.run()', sort='tottime')
model.printEntropy()
# print model.cmplxt[1:-1:50]
# print model.percentageAlive[1:-1:50]
# print model.numExtinct
print model.avalanceLengths
# model.printEntireMatrix()
