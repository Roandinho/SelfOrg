class Grid(list):
    """Special class for Grid, needed for optimization reasons"""

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
