# Implement the fault-tolerant error correction of [[7,1,3]] code using Shor's scheme. 

from utility import *

# Perform weight-1 Pauli correction according to the syndromes of six stabilizers.
def correctErrorsUsingSyndromes(errors, syndromes): 
  xsyndrome = (syndromes[0]<<2) + (syndromes[1]<<1) + syndromes[2]
  if xsyndrome:
    errors.z ^= 1<<(xsyndrome-1)
  zsyndrome = (syndromes[3]<<2) + (syndromes[4]<<1) + syndromes[5]
  if zsyndrome:
    errors.x ^= 1<<(zsyndrome-1)

# Extract the syndromes of six stabilizers using one qubit a time.
def extractSyndromes(errors, errorRates): 
  syndromes = [0 for i in range(6)]
  prepX(7, errors, errorRates)
  cnot(7, 6, errors, errorRates)
  cnot(7, 5, errors, errorRates)
  cnot(7, 4, errors, errorRates)
  cnot(7, 3, errors, errorRates)
  syndromes[0] = measX(7, errors, errorRates)
  prepX(7, errors, errorRates)
  cnot(7, 6, errors, errorRates)
  cnot(7, 5, errors, errorRates)
  cnot(7, 2, errors, errorRates)
  cnot(7, 1, errors, errorRates)
  syndromes[1] = measX(7, errors, errorRates)
  prepX(7, errors, errorRates)
  cnot(7, 6, errors, errorRates)
  cnot(7, 4, errors, errorRates)
  cnot(7, 2, errors, errorRates)
  cnot(7, 0, errors, errorRates)
  syndromes[2] = measX(7, errors, errorRates)
  prepZ(7, errors, errorRates)
  cnot(6, 7, errors, errorRates)
  cnot(5, 7, errors, errorRates)
  cnot(4, 7, errors, errorRates)
  cnot(3, 7, errors, errorRates)
  syndromes[3] = measZ(7, errors, errorRates)
  prepZ(7, errors, errorRates)
  cnot(6, 7, errors, errorRates)
  cnot(5, 7, errors, errorRates)
  cnot(2, 7, errors, errorRates)
  cnot(1, 7, errors, errorRates)
  syndromes[4] = measZ(7, errors, errorRates)
  prepZ(7, errors, errorRates)
  cnot(6, 7, errors, errorRates)
  cnot(4, 7, errors, errorRates)
  cnot(2, 7, errors, errorRates)
  cnot(0, 7, errors, errorRates)
  syndromes[5] = measZ(7, errors, errorRates)
  return syndromes

# Prepare four-qubit cat state, following circuit in Section III A. Postselect on measuring trivial verification qubit. 
def prepCat(errors, errorRates, verbose):
  z = 1
  while(z):
    prepX(7, errors, errorRates)
    prepZ(8, errors, errorRates)
    prepZ(9, errors, errorRates)
    prepZ(10, errors, errorRates)
    cnot(7, 8, errors, errorRates)
    cnot(7, 9, errors, errorRates)
    cnot(7, 10, errors, errorRates)
    prepZ(11, errors, errorRates)
    cnot(7, 11, errors, errorRates)
    cnot(8, 11, errors, errorRates)
    z = measZ(11, errors, errorRates)
    if z&verbose: print "cat fail"

# Measure six stabilizers in turn, using cat states. Whenever have non-trivial syndrome, measure all stabilizers again with bare ancilla, because there is no fault anymore.
def correctErrors(errors, errorRates, verbose=False):
  if verbose: print "starting syndrome0"
  prepCat(errors, errorRates, verbose)
  cnot(7, 3, errors, errorRates)
  cnot(8, 4, errors, errorRates)
  cnot(9, 5, errors, errorRates)
  cnot(10, 6, errors, errorRates)
  if (measX(7, errors, errorRates)^measX(8, errors, errorRates)^measX(9, errors, errorRates)^measX(10, errors, errorRates))==1:
    if verbose: print "syndrome0"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    correctErrorsUsingSyndromes(errors, syndromes)
    return 1
  if verbose: print "starting syndrome1"
  prepCat(errors, errorRates, verbose)
  cnot(7, 1, errors, errorRates)
  cnot(8, 2, errors, errorRates)
  cnot(9, 5, errors, errorRates)
  cnot(10, 6, errors, errorRates)
  if (measX(7, errors, errorRates)^measX(8, errors, errorRates)^measX(9, errors, errorRates)^measX(10, errors, errorRates))==1:
    if verbose: print "syndrome1"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    correctErrorsUsingSyndromes(errors, syndromes)
    return 1
  if verbose: print "starting syndrome2"
  prepCat(errors, errorRates, verbose)
  cnot(7, 0, errors, errorRates)
  cnot(8, 2, errors, errorRates)
  cnot(9, 4, errors, errorRates)
  cnot(10, 6, errors, errorRates)
  if (measX(7, errors, errorRates)^measX(8, errors, errorRates)^measX(9, errors, errorRates)^measX(10, errors, errorRates))==1:
    if verbose: print "syndrome2"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    correctErrorsUsingSyndromes(errors, syndromes)
    return 1
  if verbose: print "starting syndrome3"
  prepCat(errors, errorRates, verbose)
  cz(7, 3, errors, errorRates)
  cz(8, 4, errors, errorRates)
  cz(9, 5, errors, errorRates)
  cz(10, 6, errors, errorRates)
  if (measX(7, errors, errorRates)^measX(8, errors, errorRates)^measX(9, errors, errorRates)^measX(10, errors, errorRates))==1:
    if verbose: print "syndrome3"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    correctErrorsUsingSyndromes(errors, syndromes)
    return 1
  if verbose: print "starting syndrome4"
  prepCat(errors, errorRates, verbose)
  cz(7, 1, errors, errorRates)
  cz(8, 2, errors, errorRates)
  cz(9, 5, errors, errorRates)
  cz(10, 6, errors, errorRates)
  if (measX(7, errors, errorRates)^measX(8, errors, errorRates)^measX(9, errors, errorRates)^measX(10, errors, errorRates))==1:
    if verbose: print "syndrome4"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    correctErrorsUsingSyndromes(errors, syndromes)
    return 1
  if verbose: print "starting syndrome5"
  prepCat(errors, errorRates, verbose)
  cz(7, 0, errors, errorRates)
  cz(8, 2, errors, errorRates)
  cz(9, 4, errors, errorRates)
  cz(10, 6, errors, errorRates)
  if (measX(7, errors, errorRates)^measX(8, errors, errorRates)^measX(9, errors, errorRates)^measX(10, errors, errorRates))==1:
    if verbose: print "syndrome5"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    correctErrorsUsingSyndromes(errors, syndromes)
    return 1
  return 0

# Find least weight representation modulo stabilizers.
def weight(errors):
  return bin((errors.x | errors.z) & ((1 << 7) - 1)).count("1")

def reduceError(errors): 
  stabilizers = \
  [[(1<<6)+(1<<5)+(1<<4)+(1<<3),0], \
  [(1<<6)+(1<<5)+(1<<2)+(1<<1),0], \
  [(1<<6)+(1<<4)+(1<<2)+(1<<0),0], \
  [0,(1<<6)+(1<<5)+(1<<4)+(1<<3)], \
  [0,(1<<6)+(1<<5)+(1<<2)+(1<<1)], \
  [0,(1<<6)+(1<<4)+(1<<2)+(1<<0)], \
  ]
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
    if (errorsCopy.x & ((1<<7)-1)) or (errorsCopy.z & ((1<<7)-1)): 
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

