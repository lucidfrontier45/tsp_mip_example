# Example program to solve TSP as mixed integer programming

## How to run

The project is managed by uv.
You can use the following command to automatically download dependencies and run the main program.

```sh
uv run main.py berlin52.csv
```

CP-SAT is used as the default solver. You can use `--solver` option to choose other solvers e.g. 'highs' and 'scip'.

## Problem structure


- $i$ : index of city
- $d_{ij}$ : distance between cities $i$ and $j$
- $x_{ij}$ : binary decision variable, indicates if edge between $i$ and $j$ is used or not.
- $u_i$ : integer auxiliary variable, indicates thata $i$ is visited at $u_i$th order. 
- $D(x) = \sum_{i,j}d_{ij}x_{ij}$ : objective function, total distance
- $X_{i.} = \sum_j x_{ij} = 1$ : constraint, there must be one edge from $i$
- $X_{.j} = \sum_i x_{ij} = 1$ : constraint, there must be one edge to $j$
- if $x_{ij} = 1$ then $u_i - u_j = -1$ otherwise $u_i - u_j \leq M$, wher $M$ is a sufficiently large number. In this case, $M$ can be set to $N-1$ where $N$ is total number of cities. -> $u_i - u_j \leq (1-x_{ij})N -1$. This constraint, known as the Miller–Tucker–Zemlin (MTZ) formulation, is used for sub tour elimination.