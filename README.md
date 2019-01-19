# The Max Clique Problem

This is an implementation of the [EWLS, EWCC](https://www.sciencedirect.com/science/article/pii/S0004370211000427), [NuMVC](https://arxiv.org/abs/1402.0584) algorithm to solve the max clique problem as a class final project.

Input example:

the first line is `N M` for N vertices and M edges; followed by M lines as an undirected edge `a b`.

```bash
4 4 
1 2 
3 4
2 3
1 3
```

Output:

The first line is the size of the found max clique.

The second line is the indexes of vertices in the found clique.

```bash
3
1 2 3
```

