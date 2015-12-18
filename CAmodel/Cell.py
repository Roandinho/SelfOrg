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
        self.conway_n = None
        self.preys = []
        self.preds = []
        self.rest = []

    def setConwayNeighbours(self,cells):
        self.conway_n = cells

    def setNeighbours(self, cells = None):
        """ Set the neighbours, when cells is None, only update them"""
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
        """Most used function, when no alive predators or preys are present, do
        Conway"""
        n_alive = 0
        if self.State:
            for cell in self.conway_n:
                if cell.State:
                    n_alive += 1
                if n_alive > 3:
                    self.nextState = False
                    return
            if n_alive < 2:
                self.nextState = False
                return
        else:
            for cell in self.conway_n:
                if cell.State:
                    n_alive += 1
                if n_alive > 3:
                    self.nextState = False
                    return
            if n_alive < 3:
                self.nextState = False
                return
        self.nextState = True

    def updateStep(self):
        """calculate self.nextState"""

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
