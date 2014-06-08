#!/usr/bin/python

from __future__ import division
from numpy import log, exp, divide, random as r
from scipy.optimize import minimize
from scipy import array
import csv


def lognorm(x):
    return log(x + 1)

def P(params, xl_norm, yl_norm):
    """calculates sigmoid function"""
    alfa, beta, gama = params[0], params[1], params[2]
    z = alfa * xl_norm + beta * yl_norm + gama
    return divide(exp(z), 1 + exp(z))

def loglikeliPxy(params, A, N):
    """
    Returns the value of negative log likelihood
    and the derivatives
    with respect to alpha, beta and gama (in params)
    """
    dera, derb, derg = 0, 0, 0
    ret = 0

    for (x, y), val in A.iteritems():
        xl_norm = lognorm(x)
        yl_norm = lognorm(y)
        P_val = P(params, xl_norm, yl_norm)
        _val = 0
        if (x, y) in N:
            _val = N[x, y]
        ret += val * log(P_val) + _val * log(1 - P_val)
        delta = val * (1 - P_val) - _val * P_val
        dera += delta * xl_norm    # log(x + 1)
        derb += delta * yl_norm    # log(y + 1)
        derg += delta

    return ret * (-1),  \
        array([-dera, -derb, - derg])

def logreg(A, N):
    """
    Logistic regression estimator of alpha,beta, gama.
    Returns the array of estimated parameters
    """
    print "...running bfgs ..."
    initParams = [r.random_sample(), r.random_sample(), r.random_sample()]
    print initParams
    results = minimize(loglikeliPxy, initParams, args=(A, N),
      method = 'BFGS', jac = True, options={'maxIter': 5000, 'disp': True})
    print results.x
    return results.x

def min_referendum():
    """ Logistic regression estimator that reads
    dictionaries of A(x,y) and N(x, y) from file
    """

    A, N = dict(), dict()
    print "ucitavam N(x,y) i _N(x, y)"
    with open('Axy.csv', 'rb') as fp:
        csv_reader = csv.reader(fp)
        for row in csv_reader:
            A[int(row[0]), int(row[1])] = int(row[2])

    with open('Nxy.csv', 'rb') as fp:
        csv_reader = csv.reader(fp)
        for row in csv_reader:
            N[int(row[0]), int(row[1])] = int(row[2])

    return logreg(A, N)
