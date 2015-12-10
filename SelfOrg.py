#!/usr/bin/pypy

import random as random
import math


class Cell:
    """Store the information about species in cell"""

    def __init__(self, model, Level, State=True):
        """Initialize a new cell object

        Level -- identifier of the species
        State -- indicating wheither alive (True) or dead (False)
        """
        self.model = model
        self.Level = Level
        self.State = State
        self.nextState = State
        self.neighbours = None
        self.preys = []
        self.preds = []
        self.rest = []

    def setNeighbours(self, cells = None):
        self.preys = []
        self.preds = []
        self.rest = []
        if cells:
            self.neighbours = cells
        if self.model.hierarchy:
            for cell in self.neighbours:
                if cell.Level <= self.Level + self.model.hierarchyRange and cell.Level > self.Level:
                    self.preds.append(cell)
                if cell.Level >= self.Level - self.model.hierarchyRange and cell.Level < self.Level:
                    self.preys.append(cell)
                else:
                    self.rest.append(cell)
        else:
            for cell in self.neighbours:
                a = self.model.interactions[cell.Level][self.Level]
                b = self.model.interactions[self.Level][cell.Level]
                if a < 0 and b > 0:
                    self.preys.append(cell)
                if a > 0 and b < 0:
                    self.preds.append(cell)
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


class Grid(list):

    _cells = None

    def __repr__(self):
        return "".join(map(repr, self.cells))

    @property
    def cells(self):
        if self._cells:
            return self._cells
        else:
            self._cells = []
            for row in self:
                for cell in row:
                    self._cells.append(cell)
            return self._cells

    def __hash__(self):
        return hash(self.__repr__())

    def __str__(self):
        return self.__repr__()


class ConvergedException(Exception):
    pass
class AllDeadException(Exception):
    pass


class CAmodel:

    def __init__(self, M, n, steps, extProbs=None, perturbimpact=10,
                 hierarchy=True, replaceExtinct=False, updateInteractions=False,
                 hierarchyRange = 1, neighbourhood_level=1):
        # store variables
        self.calcConvergence = False
        self.hierarchyRange = hierarchyRange
        self.neighbourhood_level = neighbourhood_level
        self.hierarchy = hierarchy
        self.replaceExtinct = replaceExtinct
        self.updateInteractions = updateInteractions
        self.withExtinction = False
        self.interactions = None
        self.previousStates = dict()
        self.perturbimpact = perturbimpact
        self.dead = [0]*M
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
        # the length of a avalanche when it occurs
        self.avalanceLengths = []
        # the number of extinctions per timestep (when they occur)
        self.extinctions = []

        # create grid
        self.grid = self.initializeGrid(M, n)
        self.grid = self.reorderGrid(self.grid)
        self.interactions = self.setInteractions()
        self.setCellNeighbours()

    def updateInteraction(self):
        for row in self.interactions:
            row[random.randint(0, len(row) -1)] = random.uniform(-1,1)

    def setInteractions(self):
        interactions = []
        for _ in range(self.M):
            row = []
            for _ in range(self.M):
                row.append(random.uniform(-1,1))
            interactions.append(row)
        return interactions
                
    def initializeGrid(self, M, n):
        """Grid initizialed by funtions startLevel and startState
        Function takes the number of species M and dims of the
        2-D lattice (n by n)"""
        grid = Grid()
        for i in range(n):
            row = []
            for j in range(n):
                row.append(Cell(self,self.startLevel(M), self.startState()))
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
                        new_row.append(Cell(self, i, cell.State))
                        break
                else:
                    print ("error")
            new_Grid.append(new_row)
        return new_Grid

    def setCellNeighbours(self):
        for i in range(n):
            for j in range(n):
                neighbours = []
                for x in range(-self.neighbourhood_level, self.neighbourhood_level+1):
                    for y in range(-self.neighbourhood_level, self.neighbourhood_level+1):
                        if not x and not y:
                            continue
                        neighbours.append(self.grid[(i + x) % n ][(j + y) % n])
                self.grid[i][j].setNeighbours(neighbours)

    def startLevel(self, M):
        """Returns a completely random level (int) to start at"""
        return random.randint(0, M-1)

    def startState(self):
        """Returns True or False indicating dead or alive, randomly"""
        return True
        return bool(random.randint(0, 2))

    def updateGrid(self, timestep):

        # keep track of those species who are alive at this timestep
        alive = self.dead[:]

        # calculate if species is forced to extinct
        forcedExtinct = []
        if self.withExtinction:
            forcedExtinct = self.calcForced()

        # percentage alive
        percAlive = 0

        # UPDATE, i.e.: calc next grid and then copy into grid
        for cell in self.grid.cells:
            cell.updateStep(forcedExtinct)
            if cell.State:
                percAlive += 1

        # Keep track of the size of a avalance
        if self.calcConvergence:
            temp = []
            for i in range(len(self.not_avalanched)):
                if self.not_avalanched[i].State == self.not_avalanched[i].nextState:
                    temp.append(self.not_avalanched[i])
            self.not_avalanched = temp


        for cell in self.grid.cells:
            cell.State = cell.nextState

            # check for extinction
            if cell.State:
                alive[cell.Level] = 1

        
        died = 0
        #if a Level does not exist anymore, it is extinct, remove it and
        #randomly assign a new level, also mark it as dead
        if self.replaceExtinct:
            for i in range(len(alive)):
                if not alive[i]:
                    died += 1
                    self.dead[i] = 1
                    # Set new level for all dead cells of this level
                    for cell in self.grid.cells:
                        if cell.Level == i:
                            #nasty but quick
                            newLevel = random.randint(0, self.M  - 1 - sum(self.dead))
                            j = 0
                            for i in range(newLevel):
                                while True:
                                    if not self.dead[j]:
                                        j += 1
                                        break
                                    else:
                                        j += 1
                            cell.Level = j
                            cell.State = True
                            cell.nextState = True
                    for cell in self.grid.cells:
                        cell.setNeighbours()
        else:
            died = len(alive)-sum(alive)

        if self.updateInteractions:
            self.updateInteraction()
            for cell in self.grid.cells:
                cell.setNeighbours()



        #todo, append after update statistics?
        if self.calcConvergence:
            st = str(self.grid)
            try:
                self.previousStates[st]
            except:
                self.previousStates[st] = 1
            else:
                raise ConvergedException()

        # update statistics
        #if died:
        #    self.extinctions.append(died)
        #    if sum(self.extinctions) == self.M:
        #        raise AllDeadException()
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
        for cell in self.grid.cells:
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
        for row in self.grid:
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


    def printCsvEntropy(self):
		import csv
		myfile = open("./data_analysis/data/output_entropy_test.csv", 'wb')
		wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
		wr.writerow(self.cmplxt)
    
    def printCsvAvalanche(self):
		import csv
		myfile = open("./data_analysis/data/output_ava_length_test.csv", 'wb')
		wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
		wr.writerow(self.avalanceLengths)

    def perturb(self):
        for _ in range(self.perturbimpact):
            i = random.randint(0, self.n-1)
            j = random.randint(0, self.n-1)
            s = random.randint(0, self.M-1)
            #toggle/perturb
            self.grid[i][j].State = not self.grid[i][j].State
            #self.grid[i][j].Level = s
        self.setCellNeighbours()


    def run(self, perturb=False, singleRuns=False):
        if perturb:
            self.calcConvergence=True
            self.single_run_avalance = []
            i = 0
            j = 0
            sr_steps = 0
            self.not_avalanched = self.grid.cells
            while True:
                try:
                    self.updateGrid(j)
                except ConvergedException:
                    self.perturb()
                    self.previousStates=dict()
                    i += 1
                    self.avalanceLengths.append((j, (n*n)-len(self.not_avalanched)))
                    print ((j, (n*n)-len(self.not_avalanched)))
                    if singleRuns and i == 2:
                        self.single_run_avalance.append((j, (n*n)-len(self.not_avalanched)))
                        #reinit grid, custom counter for steps
                        sr_steps += 1
                        i = 0
                        if sr_steps == self.steps:
                            i = self.steps
                        # create grid
                        self.grid = self.initializeGrid(M, n)
                        self.grid = self.reorderGrid(self.grid)
                        self.interactions = self.setInteractions()
                        self.setCellNeighbours()
                    j = 0
                    self.not_avalanched = self.grid.cells
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
steps = 1000 # number of steps
extProbs = dict()  # dictionary containing extinction probablities
# extProbs[16] = 0.01
# extProbs[90] = 0.001

model = CAmodel(M, n, steps, perturbimpact=100, hierarchy=True,
        replaceExtinct=False, updateInteractions=False, hierarchyRange=20,
        neighbourhood_level = 3)
# model.printMatrix()
model.run(perturb=False, singleRuns=False)
# import cProfile
# cProfile.run('model.run()', sort='tottime')
model.printEntropy()
model.printCsvEntropy()
model.printCsvAvalanche()
# print model.cmplxt[1:-1:50]
# print model.percentageAlive[1:-1:50]
# print model.numExtinct
#print model.avalanceLengths
#print model.single_run_avalance
#print model.extinctions
# print model.dead
# model.printEntireMatrix()
