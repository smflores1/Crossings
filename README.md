# Crossings
This repository contains the simulations supporting predictions in our articles on crossing probabilities for critical lattice models inside polygons: [arxiv.org/abs/1608.00170](https://arxiv.org/abs/1608.00170) and [arxiv.org/abs/1407.8163](https://arxiv.org/abs/1407.8163).

The simulations section of our article explains the code that you will find in this repository. We have tried to name these files in a self-evident way:
* 'FkQ2' means the Q=2 critical FK (random cluster) model, 'SpinQ2' means the Q=2 Potts model, and 'Perc' means percolation.
* 'Rec' means rectangle, and 'Hex' means hexagon. 
* 'Indep' means independently wired boundary conditions, 'Mixed' means mixed boundary conditions, and 'Mutual' means mutually wired boundary conditions.

The aspect ratio controlling the shape of the rectangle or hexagon can be toggled between one of 33 (almost) evenly spaced values by changing the value of INTEGER, which must be an integer between 1 and 33 inclusive.
