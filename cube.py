#!/usr/bin/env python3
# Permutation approach to Rubik's cube class
# (c) 2018 youngsoo baek

"""
This file contains the class Cube that represents a Rubik's Cube in real life.

The Cube object has two attributes, a list of moves that has affected it and a list of strings that represent
each color for one of the 54 possible locations on a cube. This array of strings can be thought of as 
representing the current permutation of the cube. Each method swaps the strings in the list and thus changes
the permutation of the cube.
"""

__all__ = [ 'Cube' ]

def _map(l,lag,*pairs):
    """A function to implement permutation changes of an array l.
    Takes tuples of possible indices in l: (i,j,k,...). The ith item of l is mapped to the jth entry if lag=1.
    lag argument sets the lag between the entry and the output index. Used to control
    the number of cube rotations.
    Each tuple can be thought of as a 'cyclic notation' of the permutation of indices of l.

    >>> l = [1,2,3,4,5]
    >>> _map(l,1, (0,1), (2,3,4))
    >>> l
    [2, 1, 5, 3, 4]
    >>> _map(l,2, (0,1), (2,3,4))
    >>> l
    [2, 1, 3, 4, 5]
    """
    cp = list(l)
    for p in pairs:
        for index in range(len(p)):
            i = p[index]
            j = p[(index+lag)%(len(p))]
            l[j] = cp[i]

class Cube(object):
    """
    A class that represents the Rubik's cube. It represents a cube by a 54-length array of strings,
    which represents a permutation state of that cube.

    The Cube class has methods that manipulate each of its six faces. The methods for six basic moves 
    are named in Singmaster notation: U, D, F, B, L, and R. For more on this, refer to: 
    https://en.wikipedia.org/wiki/Rubik%27s_Cube#Move_notation

    The map of the locations on this cube:
    Edges:
    - WR (7,14), WB (6,23), WO (5,32), WG (8,41)
    - YR (52,16), YB (53,25), YO (50,34), YG (51,43)
    - RB (15,26), BO (24,35), OG (33,44), GR (42,17)
    Corners:
    - WRB (3,11,19), WBO (2,20,28), WOG (1,29,37), WGR (4,38,10)
    - YRB (49,12,22), YBO (46,21,31), YOG (47,30,40), YGR (48,39,13)
    Centers: (0,8,16,24,32,40)
    """

    __slots__ = [ '_m','_p' ]

    def __init__(self):
        """
        The Cube object has two attributes attached to it when initialized:
        - _m is a list of tuples that represent the moves made to the cube. 
        - _p is a list of strings that represent the colors for each position in the cube 
        at the current state.
        The colors of each face are fixed to be white, red, blue, orange, green, and yellow by convention.
        """
        self._m = [] 
        self._p = []
        self._p.extend( 'w' for _ in range(3*3) )
        self._p.extend( 'r' for _ in range(3*3) )
        self._p.extend( 'b' for _ in range(3*3) )
        self._p.extend( 'orange' for _ in range(3*3) )
        self._p.extend( 'g' for _ in range(3*3) )
        self._p.extend( 'y' for _ in range(3*3) )        
        
    def __str__(self):
        """Returns the array that represents the Cube in string."""
        return "{}".format(self._p)

    def __repr__(self):
        """The repr method."""
        return "Cube()"

    def _up(self,lag):
        """Corresponds to Singmaster U."""
        _map(self._p,lag, (1,2,3,4), (5,6,7,8), (10,37,28,19), (14,41,32,23), (11,38,29,20))

    def _down(self,lag):
        """Corresponds to Singmaster D."""
        _map(self._p,lag, (46,47,48,49), (50,51,52,53), (13,22,31,40), (16,25,34,43), (12,21,30,39))

    def _left(self,lag):
        """Corresponds to Singmaster L."""
        _map(self._p,lag, (37,38,39,40), (41,42,43,44), (1,10,48,30), (8,17,51,33), (4,13,47,29))

    def _right(self,lag):
        """Corresponds to Singmaster R."""
        _map(self._p,lag, (19,20,21,22), (23,24,25,26), (2,31,49,11), (6,35,53,15), (3,28,46,12))

    def _front(self,lag):
        """Corresponds to Singmaster F."""
        _map(self._p,lag, (10,11,12,13), (14,15,16,17), (4,19,49,39), (7,26,52,42), (3,22,48,38))

    def _back(self,lag):
        """Corresponds to Singmaster B."""
        _map(self._p,lag, (28,29,30,31), (32,33,34,35), (2,37,47,21), (5,44,50,24), (1,40,46,20))

    def U(self,direction='clockwise',turn=None):
        """Rotates the Up face. Direction can be either clockwise or counterclockwise.
        1 turn corresponds to a 90 degree rotation in that direction.
        For simplicity, a counterclockwise move is stored as a clockwise move with more rotations.
        """
        if turn is None:
            turn=1
        if direction=='counterclockwise':
            turn=4-turn
        self._up(turn)
        self._m.append( ('U', turn) )

    def D(self,direction='clockwise',turn=None):
        """Rotates the Down face. Direction can be either clockwise or counterclockwise.
        1 turn corresponds to a 90 degree rotation in that direction.
        For simplicity, a counterclockwise move is stored as a clockwise move with more rotations.
        """
        if turn is None:
            turn=1
        if direction=='counterclockwise':
            turn=4-turn
        self._down(turn)
        self._m.append( ('D', turn) )

    def L(self,direction='clockwise',turn=None):
        """Rotates the Left face. Direction can be either clockwise or counterclockwise.
        1 turn corresponds to a 90 degree rotation in that direction.
        For simplicity, a counterclockwise move is stored as a clockwise move with more rotations.
        """
        if turn is None:
            turn=1
        if direction=='counterclockwise':
            turn=4-turn
        self._left(turn)
        self._m.append( ('L', turn) ) 

    def R(self,direction='clockwise',turn=None):
        """Rotates the Right face. Direction can be either clockwise or counterclockwise.
        1 turn corresponds to a 90 degree rotation in that direction.
        For simplicity, a counterclockwise move is stored as a clockwise move with more rotations.
        """
        if turn is None:
            turn=1
        if direction=='counterclockwise':
            turn=4-turn
        self._right(turn)
        self._m.append( ('R', turn) )

    def F(self,direction='clockwise',turn=None):
        """Rotates the Front face. Direction can be either clockwise or counterclockwise.
        1 turn corresponds to a 90 degree rotation in that direction.
        For simplicity, a counterclockwise move is stored as a clockwise move with more rotations.
        """
        if turn is None:
            turn=1
        if direction=='counterclockwise':
            turn=4-turn
        self._front(turn)
        self._m.append( ('F', turn) ) 

    def B(self,direction='clockwise',turn=None):
        """Rotates the Back face. Direction can be either clockwise or counterclockwise.
        1 turn corresponds to a 90 degree rotation in that direction.
        For simplicity, a counterclockwise move is stored as a clockwise move with more rotations.
        """
        if turn is None:
            turn=1
        if direction=='counterclockwise':
            turn=4-turn
        self._back(turn)
        self._m.append( ('B', turn) ) 

    def shuffle(self):
        """Shuffles the cube by making a random number of possible moves on the cube.
        The possible number of moves for shuffling is arbitrarily set between 10 and 20.
        """
        from random import randint
        bound = randint(10,20)
        moves = (self.U, self.D, self.L, self.R, self.F, self.B)
        for i in range(bound):
            move = moves[randint(0, len(moves)-1)]
            turn = randint(1,3)
            move(turn=turn)

    def solve(self):
        """Returns a solution to the current state of the cube.
        In its current incarnation: simply returns the list of moves that reverse the order of moves made
        to the cube so far (in reverse directions). All moves are directed clockwise.
        """
        solution = []
        for face, turn in self._m[::-1]:
            solution.append( (face, 4-turn) )
        return solution

    def reset(self):
        """Resets the cube."""
        self._p = Cube()._p
        self._m.clear()
