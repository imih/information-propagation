#!/usr/bin/python

import csv
from igraph import read
from numpy.random import permutation
from numpy import less, linspace
from logreg import logreg, min_referendum
from matplotlib import pyplot as plt


def getNs(unactive_ids, nodes, auth_time, window, _print=False, _small=True):
  """
  Calculate  A(x,y)  and N(x,y), where 
  A(x,y) = sum_t A(x, y, t);
  N(x,y) = sum_t N(x, y, t);

  A(x, y, t) - number of nodes that became active at time t
   when x of its neighbours were active and y authorities were active
  N(x, y, t) - number of unactive nodes that did not become active at time t
   when x of its neighbours were active and y authorities were active

  Keyword arguments:
  unactive_ids -- id's of nodes that never become active
  nodes -- dictionary that maps id of the peer node to its
       instance igraph.Vertex() for every not that eventually becomes active
  auth_time -- list of activation times for authority nodes
    (not present in the nodes dictionary)
  window -- number of units of time that are used as
    the length of the discrete time step
  _print -- if True, calculated N, _N will be written
    to files Axy.csv and Nxy.csv (default False)
  _small -- if True, nodes that will never become active
    will not be included in the calculations (default True)
  """

  node_list = [(_id, -1) for _id in unactive_ids]
  node_list += [(node["id"], int(node["regtime"]) / window)
      for node in nodes.values() if node["regtime"] != -1]
  node_list = sorted(node_list, key=lambda (x, _y): _y)
  assert(len(node_list) > 0)

  auth_time = map(lambda x: (x / window), auth_time)
  _time = set([t for (_id, t) in node_list if t != -1])
  _time |= set(auth_time)
  _time = sorted(_time)

  A, N = dict(), dict()
  bucket, which_b = dict(), dict()

  """
  bucket[x] contains id's of unactive nodes that have x active friends
  which_b[i]: bucket index of node with id = i

  Algorithm:
    initialize: add all nodes to bucket[0]
      (all nodes are unactive at time step t = 0)
    y = 0
    for every t:
     N(x,y) += len(bucket[x])  * (t - tprev)
     for every activation of node  n in t:
       x = which_b(n)
       N(x, y) -= 1
       A(x, y) + 1
     for every activation of node n in t:
       x = which_b(n)
       remove  n from bucket[which_b(n]
       for every neighbor i of n:
         remove i from bucket[which_b(i)]
         add i to bucket[which_b(i) + 1]
         which_b(i) += 1
     if t in auth_time:
       y = y + 1
  """

  bucket[0] = set([_id for (_id, t) in node_list if not _small or
    (_small and t != -1)])
  print "#nodes: " + str(len(bucket[0]))

  for _id in bucket[0]:
    which_b[_id] = 0

  lo , hi = 0, len(node_list)
  while lo < hi:
    mid = (lo + hi) // 2
    if node_list[mid][1] == -1:
      lo = mid + 1
    else:
      hi = mid
  i = lo

  y = 0
  tprev = -1
  for t in _time:
    if t > (7200 / window):  # pet dana
      break
    for x, binn in bucket.iteritems():
      N[x, y] = N.get((x, y), 0) + len(binn) * (t - tprev)

    iprev = i
    while len(node_list) > i and node_list[i][1] == t:
      _id = node_list[i][0]
      x = which_b[_id]
      N[x, y] = N.get((x, y), 0) - 1
      N[x, y] = A.get((x, y), 0) + 1
      bucket[x].remove(_id)
      i = i + 1

    i = iprev
    while len(node_list) > i and node_list[i][1] == t:
      neis = nodes[node_list[i][0]].neighbors()
      for node in neis:
        if node["id"] not in which_b:
          continue
        _x = which_b[node["id"]]
        if node["id"] in bucket[_x]:
          bucket[_x].remove(node["id"])
          if _x + 1 not in bucket:
            bucket[_x + 1] = set()
          bucket[_x + 1].add(node["id"])
          which_b[node["id"]] = _x + 1
      i = i + 1

    if t in auth_time:
      if y == 0:
        y = y + 1
    tprev = t

  ###############
  if _print:
    print "printing files..."
    with open('Axy.csv', 'wb') as fp:
      csv_writer = csv.writer(fp)
      assert(len(A.keys()) > 0)
      for (x, y), val in A.iteritems():
        csv_writer.writerow([x, y, val])

    with open('Nxy.csv', 'wb') as fp:
      csv_writer = csv.writer(fp)
      for (x, y), val in N.iteritems():
        csv_writer.writerow([x, y, val])
  return (A, N)


def calcExternalActivations(_time, bucket, which_b, node_list, nodes, window):
  """
  For every time t in set _time, the procedure calculates number of nodes that
  activate at time t and have 0 friends active at that time.
  The values are return as a dictionary with key being the timestamp.
  """
  y = dict()
  lo , hi = 0, len(node_list)
  while lo < hi:
    mid = (lo + hi) // 2
    if node_list[mid][1] == -1:
      lo = mid + 1
    else:
      hi = mid
  i = lo

  for t in _time:
    if t > (7200 / window):  # pet dana
      break
    yt = 0

    iprev = i
    while node_list[i][1] == t:
      _id = node_list[i][0]
      x = which_b[_id]
      if x == 0:
        yt = yt + 1
      i = i + 1
      bucket[x].remove(_id)

    i = iprev
    while node_list[i][1] == t:
      neis = nodes[node_list[i][0]].neighbors()
      for node in neis:
        if node["id"] not in which_b:
          continue
        _x = which_b[node["id"]]
        if node["id"] in bucket[_x]:
          bucket[_x].remove(node["id"])
          if _x + 1 not in bucket:
            bucket[_x + 1] = set()
          bucket[_x + 1].add(node["id"])
          which_b[node["id"]] = _x + 1
      i = i + 1
    y[t] = yt
  return y

def getNs2(unactive_ids, nodes, auth_time, window, _print=False, _small=True):
  """
  Calculate  N(x,y)  and _N(x,y), where N(x, y)
  N(x,y) = sum_t N(x, y, t);
  _N(x, y) = sum_t _N(x, y, t);

  N(x, y, t) - number of nodes that became active at time t
   when x of its neighbours were active and y nodes
   that have 0 friends became active
  _N(x, y, t) - number of unactive nodes
   that did not become active at time t
   when x of its neighbours were active and y nodes
   that have 0 friends became active

  Keyword arguments:
  nodes -- dictionary that maps id of the peer node to its
       instance igraph.Vertex()
  auth_time -- list of activation times for authority nodes
    (not present in the nodes dictionary)
  window -- number of units of time that are used as
    the length of the discrete time step
    (t in argument of N(x, y, t) and _N(x, y, t)
  _print -- if True, calculated N, _N will be written
    to files Axy.csv and Nxy.csv (default False)
  _small -- if True, nodes that will never become active
    will not be included in the calculations (default True)
  """

  node_list = [(_id, -1) for _id in unactive_ids]
  node_list += [(node["id"], int(node["regtime"]) / window)
      for node in nodes.values() if node["regtime"] != -1]
  node_list = sorted(node_list, key=lambda (x, _y): _y)

  times = dict()
  for node in nodes.values():
    times[node["id"]] = int(node["regtime"]) / window

  auth_time = map(lambda x: (x / window), auth_time)
  _time = set([t for (_id, t) in node_list if t != -1])
  _time |= set(auth_time)
  _time = sorted(_time)

  A, N = dict(), dict()
  bucket, which_b = dict(), dict()

  """
  bucket[x] contains id's of unactive nodes that have x active friends
  which_b[i]: bucket index of node with id = i

  Algorithm:
    initialize: add all nodes to bucket[0]
      (all nodes are unactive at time step t = 0)
    for every t:
     y = getY()
     N(x,y) += len(bucket[x])
     for every activation of node  n in t:
       x = which_b(n)
       N(x, y) -= 1
       A(x, y) + 1
     for every activation of node n in t:
       x = which_b(n)
       remove  n from bucket[which_b(n]
       for every neighbor i of n:
         remove i from bucket[which_b(i)]
         add i to bucket[which_b(i) + 1]
         which_b(i) += 1

  """

  bucket[0] = set([_id for (_id, t) in node_list if not _small or
    (_small and t != -1)])
  print "#nodes: " + str( len(bucket[0]))

  for _id in bucket[0]:
    which_b[_id] = 0

  print "racunam N-ove"

  y = calcExternalActivations(_time, bucket, which_b, node_list,  nodes, window)

  bucket[0] = set([_id for (_id, t) in node_list if not _small or
    (_small and t != -1)])

  for _id in bucket[0]:
    which_b[_id] = 0

  lo , hi = 0, len(node_list)
  while lo < hi:
    mid = (lo + hi) // 2
    if node_list[mid][1] == -1:
      lo = mid + 1
    else: 
      hi = mid
  i = lo

  tprev = -1
  for t in _time:
    if t > (7200 / window):  # pet dana
      break
    yt = y[t]

    for x, binn in bucket.iteritems():
      N[x, yt] = N.get((x, yt), 0) + len(binn) * (t - tprev)

    iprev = i
    while node_list[i][1] == t:
      _id = node_list[i][0]
      x = which_b[_id]
      N[x, yt] = N.get((x, yt), 0) - 1
      A[x, yt] = A.get((x, yt), 0) + 1
      bucket[x].remove(_id)
      i = i + 1

    i = iprev
    while node_list[i][1] == t:
      neis = nodes[node_list[i][0]].neighbors()
      for node in neis:
        if node["id"] not in which_b:
          continue
        _x = which_b[node["id"]]
        if node["id"] in bucket[_x]:
          bucket[_x].remove(node["id"])
          if _x + 1 not in bucket:
            bucket[_x + 1] = set()
          bucket[_x + 1].add(node["id"])
          which_b[node["id"]] = _x + 1
      i = i + 1
    tprev = t

  ###############
  if _print:
    print "printing files..."
    with open('Axy.csv', 'wb') as fp:
      csv_writer = csv.writer(fp)
      assert(len(A.keys()) > 0)
      for (x, y), val in A.iteritems():
        csv_writer.writerow([x, y, val])

    with open('Nxy.csv', 'wb') as fp:
      csv_writer = csv.writer(fp)
      for (x, y), val in N.iteritems():
        csv_writer.writerow([x, y, val])
  return (A, N)

def shuffle_test(gml_file, auth_time, N_getter, tempfile, _iter, window, _small):
  """
  The procedure generate a file with estimated coefficients 
  for instances of randomized dataset.
  Keyword arguments:
  gml_file --- the graph in  .gml format 
  auth_time --- list of activation times for authority  nodes
  N_getter --- function object that calculates the values of A(x, y)
    and N(x,y). Can be calc.getNs or calc.getNs2.
  tempfile --- name of the file the coefficients will be writen to
  _iter --- number of iterations - runs of the randomization and estimation 
    steps
  window --- length of the window for one time step
  _small -- if True, the estimation will be run on restricted, only 
    with nodes that eventually become active
  """
  print "...reading graph..."
  g = read(gml_file)
  nodes = dict()
  unactive_ids = []
  for node in g.vs:
    if int(node["regtime"]) == -1:
      unactive_ids.append(node["id"])
    else:
      nodes[node["id"]] = node

  f = open(tempfile, 'a')
  (A, N) = N_getter(unactive_ids, nodes, auth_time, window, False, _small)
  f.write(" ".join(str(num) for num in logreg(A, N)) + "\n")

  all_time = [node["regtime"] for node in g.vs
      if node["regtime"] != -1] + auth_time

  for j in range(0, _iter):
    print j

    all_time = permutation(all_time)
    auth_time = []

    i = 0
    for _id, node in nodes.iteritems():
      while len(auth_time) != 11 and all_time[i] != -1:
        auth_time.append(all_time[i])
        i = i + 1
      if node["regtime"] == -1:
        continue
      nodes[_id]["regtime"] = all_time[i]
      i = i + 1

    (A, N) = N_getter(unactive_ids, nodes, auth_time, window, False, _small)
    f.write(" ".join(str(num) for num in logreg(A, N)) + "\n")
  f.close()

def computeStrength(tempfile):
  """
  This procedure reads parameter estimations from file,
  calculates the strength of peer and authority influence
  and plots the frequnecies of values for alpha and beta.
  """
  org = []
  a, b, g = [], [], []
  alfaS = 0.0
  betaS = 0.0
  _iter = -1
  f = open(tempfile, 'r')
  for line in f.readlines():
    _iter = _iter + 1
    nums = line.split()
    if len(org) == 0:
      org.append(float(nums[0]))
      org.append(float(nums[1]))
      org.append(float(nums[2]))
    else:
      if float(nums[0]) < 100:
        a.append(float(nums[0]))
      if float(nums[1]) < 100 and float(nums[1]) > -100:
        b.append(float(nums[1]))
      if float(nums[2]) < 0 and float(nums[2]) > -50:
        g.append(float(nums[2]))
      if less(float(nums[0]), org[0]):
        alfaS += 1
      if less(float(nums[1]), org[1]):
        betaS += 1

  f.close()
  alfaS = alfaS / _iter
  betaS = betaS / _iter
  print alfaS, betaS
  bins = linspace(0.6, 1.4, 90)

  n, bins, patches = plt.hist(a, bins, alpha=0.5, color='b', label='alpha')
  plt.hist(b, bins, alpha=0.5,  color='g', label='beta')
  plt.legend(loc='upper right')
  plt.show()

