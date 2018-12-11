#!/usr/bin/env python3
# Final project file.
# (c) 2018 youngsoo baek

"""
This file contains functions that are needed to run a simple Rubik's Cube game.
- plotCube: Plots the Cube class object.
- interact: Controls the user's interaction with the Cube through keyboard inputs.
- RubiksCube: The 'main method' that runs the Rubik's Cube solving game.
"""
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as a3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Rectangle
from itertools import product
from cube import Cube

# A constant tuple holding position numberings. It is used to color the corresponding coordinates.
POSITIONS = (39,43,40,42,36,44,38,41,37,22,25,21,26,18,24,19,23,20,13,16,12,17,9,15,10,14,11,30,34,31,33,27,35,29,32,28,48,52,49,51,45,53,47,50,46,4,7,3,8,0,6,1,5,2 )

def plotCube(cube):
    """A function that plots the Cube class object.
    A visual representation is produced with 6 rectangular planes in a 3D space.
    Each item in the array that represents the cube's permutation state is mapped as the color
    to its corresponding locations on the 6 sides of the cube.
    """
    ax = Axes3D(plt.figure())
    cols = cube._p
    xs = [ -3, -1, 1 ]
    ys = [ -3, -1, 1 ]
    zs = [ -3, 3 ]
    dirs = [ "x", "y", "z" ]
    locs = list(product(dirs,zs,ys,xs)) 
    # product conceals the four loops, nested in order
    for i in range(len(locs)):
        col = cols[POSITIONS[i]]
        x = locs[i][3]
        y = locs[i][2]
        z = locs[i][1]
        dir = locs[i][0]
        r = Rectangle( (x,y),2,2, edgecolor='k', facecolor=col, lw=2 )
        ax.add_patch(r)
        a3d.pathpatch_2d_to_3d(r,z=z,zdir=dir)
    ax.set_xlim3d(-3,3)
    ax.set_ylim3d(-3,3)
    ax.set_zlim3d(-3,3)
    ax.grid(False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    
def interact(cube):
    """Interacting with the cube.
    Each keyboard input of the user is passed as the face to rotate, the direction of the rotation,
    and the number of 90 degree rotations to be made. 
    """
    choice = input("Choose to continue, get solution, or reshuffle [c,s,r]: ")
    if choice == "c":
        move = input("Which face to move? ")
        while move not in ('U', 'D', 'L', 'R', 'F', 'B'):
            move = input("Which face to move? ")        
        direction = input("Which direction? ")
        while direction not in ('clockwise', 'counterclockwise'):
            direction = input("Which direction? ")
        turn = int(input("How many turns? "))
        while (turn < 1) or (turn > 3):
            turn = int(input("How many turns? "))
        move = getattr(cube, move) 
        # finds a method with name 'move'
        move(direction,turn)
    elif choice=='s':
        print("These are clockwise moves that can solve the current state of the cube: {}".format(cube.solve()))
    elif choice=='r':
        cube.reset()
        cube.shuffle()
    else:
        interact(cube)

def RubiksCube():
    """Runs the Rubik's Cube game."""
    solved = False
    while (input("Start game?: ") == "yes") and (not solved):
        c = Cube()
        c.shuffle()
        plotCube(c)
        plt.ion()
        plt.show()
        plt.pause(1)
        # pause needed to update matplotlib plots with keyboard inputs
        print("The possible faces to rotate are: U,D,L,R,F,B.\n The possible directions to rotate are: clockwise, counterclockwise.\n The possible numbers of rotations are: 1,2,3.\n")
        while not solved:
            interact(c)
            plt.close()
            plotCube(c)
            plt.show()
            plt.pause(1)
            if c._p == Cube()._p:
                solved = True
        if solved:
            print("You win!")

if __name__ == '__main__':
    RubiksCube()
