# Implement the fault-tolerant error correction of [[15,7,3]] code using Shor's scheme. 

from utility import *

# Perform weight-1 Pauli correction according to the syndromes of eight stabilizers.
def correctErrorsUsingSyndromes(errors, syndromes):
  xsyndrome = (syndromes[0]<<3) + (syndromes[1]<<2) + (syndromes[2]<<1) + syndromes[3]
  if xsyndrome:
    errors.z ^= 1<<(xsyndrome-1)
  zsyndrome = (syndromes[4]<<3) + (syndromes[5]<<2) + (syndromes[6]<<1) + syndromes[7]
  if zsyndrome:
    errors.x ^= 1<<(zsyndrome-1)

# This is the matrix indexing the non-identity positions of the stabilizers 
cnotWires = [[k for k in range(15) if ((k+1)>>(3-i))&1==1] for i in range(4)]

# Extract the syndromes of eight stabilizers using one qubit a time.
def extractXSyndromes(errors, errorRates):
  syndromes = [0 for i in range(8)]
  for i in range(4):
    prepX(15, errors, errorRates)
    for j in range(8):
      cnot(15, cnotWires[i][j], errors, errorRates)
    syndromes[i] = measX(15, errors, errorRates)
  return syndromes

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

# Prepare eight-qubit cat state, following circuit in Section III A. Postselect on measuring trivial verification qubit. 
def prepCat(errors, errorRates, verbose):
  z = 1
  while(z):
    prepX(15, errors, errorRates)
    prepZ(16, errors, errorRates)
    prepZ(17, errors, errorRates)
    prepZ(18, errors, errorRates)
    prepZ(19, errors, errorRates)
    prepZ(20, errors, errorRates)
    prepZ(21, errors, errorRates)
    prepZ(22, errors, errorRates)
    cnot(15, 16, errors, errorRates)
    cnot(15, 17, errors, errorRates)
    cnot(15, 18, errors, errorRates)
    cnot(15, 19, errors, errorRates)
    cnot(15, 20, errors, errorRates)
    cnot(15, 21, errors, errorRates)
    cnot(15, 22, errors, errorRates)
    prepZ(23, errors, errorRates)
    cnot(15, 23, errors, errorRates)
    cnot(16, 23, errors, errorRates)
    z = measZ(23, errors, errorRates)
    if verbose&z: print "cat fail"

# Measure eight stabilizers in turn, using cat states. Whenever have non-trivial syndrome, measure all stabilizers again with bare ancilla, because there is no fault anymore.
def correctErrors(errors, errorRates, verbose=False):
  for i in range(4):
    prepCat(errors, errorRates, verbose)
    for j in range(8):
      cnot(15+j, cnotWires[i][j], errors, errorRates)
    syndrome = 0
    for j in range(8):
      syndrome ^= measX(15+j, errors, errorRates)
    if syndrome:
      if verbose: print "syndrome%d"%i
      syndromes = extractSyndromes(errors, errorRates)
      if verbose: print syndromes
      correctErrorsUsingSyndromes(errors, syndromes)
      return 1
  for i in range(4):
    prepCat(errors, errorRates, verbose)
    for j in range(8):
      cz(15+j, cnotWires[i][j], errors, errorRates)
    syndrome = 0
    for j in range(8):
      syndrome ^= measX(15+j, errors, errorRates)
    if syndrome:
      if verbose: print "syndrome%d"%(i+4)
      syndromes = extractSyndromes(errors, errorRates)
      if verbose: print syndromes
      correctErrorsUsingSyndromes(errors, syndromes)
      return 1
  return 0

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

