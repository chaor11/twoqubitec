# Implement the fault-tolerant error correction of [[15,7,3]] code using only two ancilla qubit. 

from utility import *

# Perform weight-1 Pauli correction according to the syndromes of eight stabilizers.
def correctErrorsUsingSyndromes(errors, syndromes):
  xsyndrome = (syndromes[0]<<3) + (syndromes[1]<<2) + (syndromes[2]<<1) + syndromes[3]
  if xsyndrome:
    errors.z ^= 1<<(xsyndrome-1)
  zsyndrome = (syndromes[4]<<3) + (syndromes[5]<<2) + (syndromes[6]<<1) + syndromes[7]
  if zsyndrome:
    errors.x ^= 1<<(zsyndrome-1)

# The matrix indexing the non-identity positions of the stabilizers 
cnotWires = [[k for k in range(15) if ((k+1)>>(3-i))&1==1] for i in range(4)]

# Extract the syndromes of X stabilizers using one qubit a time.
# For CSS codes we sometimes only have to measure X or Z stabilizers alone.
def extractXSyndromes(errors, errorRates):
  syndromes = [0 for i in range(8)]
  for i in range(4):
    prepX(15, errors, errorRates)
    for j in range(8):
      cnot(15, cnotWires[i][j], errors, errorRates)
    syndromes[i] = measX(15, errors, errorRates)
  return syndromes

# Extract the syndromes of Z stabilizers using one qubit a time.
def extractZSyndromes(errors, errorRates):
  syndromes = [0 for i in range(8)] 
  for i in range(4):
    prepZ(15, errors, errorRates)
    for j in range(8):
      cnot(cnotWires[i][j], 15, errors, errorRates)
    syndromes[i+4] = measZ(15, errors, errorRates)
  return syndromes

def extractSyndromes(errors, errorRates):
  xsyn = extractXSyndromes(errors, errorRates)
  zsyn = extractZSyndromes(errors, errorRates)
  return [xsyn[i]+zsyn[i] for i in range(8)]

# Permutation for the CNOT wires extracting a single stabilizer syndrome. Note that it ranges from 1 to 8.
perm = [1,2,3,5,4,7,6,8]

# flagWire contains four matrices. Each matrix contains 8 distinct syndromes of correlated error when measuring a stabilizer.
flagWire = [[[0 for k in range(8)] for j in range(4)] for i in range(4)]
f = [[0 for k in range(8)] for j in range(4)] 
temp = [0 for j in range(8)]
for i in range(7):
  temp[6-i] = (7+perm[7-i]) ^ temp[7-i]
for i in range(4):
  for j in range(8):
    f[i][j] = (temp[j]>>(3-i))&1
flagWire = [f,  [f[1],f[0],f[2],f[3]],  [f[1],f[2],f[0],f[3]],  [f[1],f[2],f[3],f[0]]]

# Implement the error correction procedure in Section III in the paper. For example, the circuit for measurement of IIIIIIIZZZZZZZZ follows FIG.3 (c).
def correctErrors(errors, errorRates, verbose=False):
  for i in range(4):
    if verbose: print "starting syndrome%d" % i
    prepX(15, errors, errorRates)
    prepZ(16, errors, errorRates)
    cnot(15, cnotWires[i][perm[0]-1], errors, errorRates)
    cnot(15, 16, errors, errorRates)
    for j in range(6):
      cnot(15, cnotWires[i][perm[j+1]-1], errors, errorRates)
    cnot(15, 16, errors, errorRates)
    cnot(15, cnotWires[i][perm[7]-1], errors, errorRates)
    syndrome = measX(15, errors, errorRates)
    flag = measZ(16, errors, errorRates)
    if flag:
      if verbose: print "flag%d"% i
      syndromes = extractZSyndromes(errors, errorRates)
      if verbose: print "corrX:", syndromes
      for j in range(8):
        if syndromes == [0,0,0,0,flagWire[i][0][j],flagWire[i][1][j],flagWire[i][2][j],flagWire[i][3][j]]:
          for k in range(j):
            errors.x ^= 1<<cnotWires[i][perm[k]-1]
      syndromes = extractXSyndromes(errors, errorRates)
      if verbose: print "Z:", syndromes
      correctErrorsUsingSyndromes(errors, syndromes)
      return
    elif syndrome:
      if verbose: print "syndrome%d"% i
      syndromes = extractSyndromes(errors, errorRates)
      if verbose: print syndromes
      correctErrorsUsingSyndromes(errors, syndromes)
      return

  for i in range(4):
    if verbose: print "starting syndrome%d" % (i+4)
    prepZ(15, errors, errorRates)
    prepX(16, errors, errorRates)
    cnot(cnotWires[i][perm[0]-1], 15, errors, errorRates)
    cnot(16, 15, errors, errorRates)
    for j in range(6):
      cnot(cnotWires[i][perm[j+1]-1], 15, errors, errorRates)
    cnot(16, 15, errors, errorRates)
    cnot(cnotWires[i][perm[7]-1], 15, errors, errorRates)
    syndrome = measZ(15, errors, errorRates)
    flag = measX(16, errors, errorRates)
    if flag:
      if verbose: print "flag%d"% (i+4)
      syndromes = extractXSyndromes(errors, errorRates)
      if verbose: print "corrZ:", syndromes
      for j in range(8):
        if syndromes == [flagWire[i][0][j],flagWire[i][1][j],flagWire[i][2][j],flagWire[i][3][j],0,0,0,0]:
          for k in range(j):
            errors.z ^= 1<<cnotWires[i][perm[k]-1]
      syndromes = extractZSyndromes(errors, errorRates)
      if verbose: print "X:", syndromes
      correctErrorsUsingSyndromes(errors, syndromes)
      return
    elif syndrome:
      if verbose: print "syndrome%d"% (i+4)
      syndromes = extractSyndromes(errors, errorRates)
      if verbose: print syndromes
      correctErrorsUsingSyndromes(errors, syndromes)
      return

# Find least weight representation modulo stabilizers.
def weight(errors):
  return bin((errors.x | errors.z) & ((1 << 15) - 1)).count("1")

stabilizers = [[0,0] for i in range(8)]
for i in range(4):
  for j in range(8):
    stabilizers[i][0] += 1<<cnotWires[i][j]
for i in range(4):
  for j in range(8):
    stabilizers[i+4][1] += 1<<cnotWires[i][j]

def reduceError(errors): 
  bestErrors = Errors(errors.x, errors.z)
  bestWeight = weight(bestErrors)
  trialErrors = Errors(0, 0)
  for k in range(1, 1<<(len(stabilizers))):
    trialErrors.x = errors.x
    trialErrors.z = errors.z
    for digit in range(len(stabilizers)):
      if (k>>digit)&1: 
        trialErrors.x ^= stabilizers[digit][0]
        trialErrors.z ^= stabilizers[digit][1]
    if weight(trialErrors) < bestWeight: 
      bestErrors.x = trialErrors.x
      bestErrors.z = trialErrors.z
      bestWeight = weight(bestErrors)
  return bestErrors

# Run consecutive trials of error correction with physical error rate of gamma, and count the number of failures, i.e., when the trialing error is not correctable by perfect error correction.
# The logical error rate is calculated as the ratio of failures over trials. 
def simulateErrorCorrection(gamma, trials): 
  errors = Errors(0, 0)
  errorsCopy = Errors(0, 0)
  
  errorRates0 = ErrorRates(0, 0, 0)
  errorRates = ErrorRates((4/15.)*gamma, gamma, (4/15.)*gamma)
  
  failures = 0
  for k in xrange(trials): 
    correctErrors(errors, errorRates)
    errorsCopy.x = errors.x
    errorsCopy.z = errors.z
    correctErrors(errorsCopy, errorRates0)
    errorsCopy = reduceError(errorsCopy)
    if (errorsCopy.x & ((1<<15)-1)) or (errorsCopy.z & ((1<<15)-1)): 
      failures += 1
      errors.x = 0
      errors.z = 0
  print failures

# Wrapper function for the plot. More trials are needed for small gammas due to the confidence interval.
gammas = [10**(i/10.-4) for i in range(21)]
for i in range(10):
  print "gamma=10^(%d/10-4), trials=10^7"% i
  simulateErrorCorrection(gammas[i], 10**7)
for i in range(11):
  print "gamma=10^(%d/10-4), trials=10^6"% (i+10)
  simulateErrorCorrection(gammas[i+10], 10**6)

