# :cloud: :cloud: cellnet :cloud: :cloud:

cellnet is a python 2 package for transforming natural cellular networks.
Examples for the natural cellular networks addressed here are foams or biological cells, but not mobile phone infrastructure.
cellnet was developed for cellular networks of clouds and its application to this case it documented in the following publication:

F. Glassmeier & G. Feingold (2017): Network approach to patterns in stratocumulus clouds. (in press) 

## Examples

Browse the following ipython notebooks (no installation required) for an illustration of cellnet's functionality:

- [division of a cell in a cellular network](examples/illustrate_celldivision.ipynb)
- [merging of two cells in a cellular network](examples/illustrate_cellmerging.ipynb)
- [illustration of the network transformation governing open cloud cells](examples/illustrate_opencells.ipynb)
- [evolution of a cellular network under repeated flips of random edges](examples/random_edge_flip_evolution.ipynb)
- [evolution of a cellular network under repeated edge flips biased to represent cloud networks](examples/biased_edge_flip_evolution.ipynb)
- [evolution of a closed-cell cloud network](examples/closed_cell_Sc_cloud_evolution.ipynb)
- [evolution of an open-cell cloud network](examples/open_cell_Sc_cloud_evolution.ipynb)

## Installation

cellnet is a Python 2 package.

First install the following Python packages that cellnet depends on:

```
pip install jupyter numpy scipy "matplotlib<2" networkx pyvoro
```

Then download the zip file [here](https://github.com/fglassmeier/cellnet/archive/master.zip), unpack it, and run the following code to browse the example notebooks (which will open in your browser):

```
cd examples
jupyter notebook
```

:copyright: 2017 Franziska Glassmeier
