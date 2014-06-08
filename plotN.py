#!/usr/bin/python

from numpy import genfromtxt, arange, meshgrid, exp
import matplotlib.pyplot as plt
from igraph import read

def Nplot(filename):
  N = dict()
  my_data = genfromtxt(filename, delimiter=',')
  _x = my_data[:, 0]
  _y = my_data[:, 1]
  _z = my_data[:, 2]
  for (i, j, k) in zip(_x, _y, _z):
    N[(i, j)] = k

  x = arange(0,  51, 1)
  y = arange(0, _y.max() + 1, 1)
  X, Y = meshgrid(x, y)

  Z = [[0 if (i, j) not in N else N[(i, j)] for i in x] for j in y]
  plt.matshow(Z, cmap=plt.cm.Greys)
  plt.colorbar()
  plt.show()

def plotyN(filename):
  N = dict()
  my_data = genfromtxt(filename, delimiter=',')
  _x = my_data[:, 0]
  _y = my_data[:, 1]
  plt.hist(_y)
  plt.show()


def plotDegDistr(gml_file):
  g = read(gml_file)
  print g.summary()

  xs, ys = zip(*[(left, count) for left, _, count
    in g.degree_distribution().bins()])
  _sum = sum(ys)
  _novy = [float(ys[i]) / _sum for i in range(0, len(ys))]

  plt.stem(xs, _novy)
  plt.yscale("log")
  plt.xscale("log")
  plt.ylabel('p(k)')
  plt.xlabel('k')
  plt.show()


def sigmoid(z, a):
  return exp(a + z) / (1. + exp(a + z))

def plotSigmoid():
  x = arange(-6, 6, .01)
  S = sigmoid(x, 0)
  plt.plot(x, S)
  plt.xlabel('z')
  plt.ylabel('sigma(z)')
  plt.show()


def plotSigmoid2():
  x = arange(-6, 6, .01)
  plt.plot(x, sigmoid(x, 0))
  plt.plot(x, sigmoid(x, 1))
  plt.plot(x, sigmoid(x, 2))
  plt.plot(x, sigmoid(x, 3))
  plt.xlabel('z')
  plt.ylabel('sigma(a + z)')
  plt.show()

