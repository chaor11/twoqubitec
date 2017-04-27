# Implement the fault-tolerant error correction of [[5,1,3]] code using only two ancilla qubit. 

from utility import *

# Perform weight-1 Pauli correction according to the syndromes of four stabilizers.
def correctErrorsUsingSyndromes(errors, syndromes): 
  if syndromes == [0,0,0,0]: 
    pass
  elif syndromes == [0,0,0,1]: 
    errors.x ^= 1<<0
  elif syndromes == [1,0,1,1]: 
    errors.x ^= 1<<0
    errors.z ^= 1<<0
  elif syndromes == [1,0,1,0]: 
    errors.z ^= 1<<0
  elif syndromes == [1,0,0,0]: 
    errors.x ^= 1<<1
  elif syndromes == [1,1,0,1]: 
    errors.x ^= 1<<1
    errors.z ^= 1<<1
  elif syndromes == [0,1,0,1]: 
    errors.z ^= 1<<1
  elif syndromes == [1,1,0,0]: 
    errors.x ^= 1<<2
  elif syndromes == [1,1,1,0]: 
    errors.x ^= 1<<2
    errors.z ^= 1<<2
  elif syndromes == [0,0,1,0]: 
    errors.z ^= 1<<2
  elif syndromes == [0,1,1,0]: 
    errors.x ^= 1<<3
  elif syndromes == [1,1,1,1]: 
    errors.x ^= 1<<3
    errors.z ^= 1<<3
  elif syndromes == [1,0,0,1]: 
    errors.z ^= 1<<3
  elif syndromes == [0,0,1,1]: 
    errors.x ^= 1<<4
  elif syndromes == [0,1,1,1]: 
    errors.x ^= 1<<4
    errors.z ^= 1<<4
  elif syndromes == [0,1,0,0]: 
    errors.z ^= 1<<4

# Extract the syndromes of four stabilizers using one qubit a time.
def extractSyndromes(errors, errorRates): 
  syndromes = [0 for i in range(4)]
  for j in range(4): 
    prepZ(5, errors, errorRates)
    dualcz(j%5, 5, errors, errorRates)
    cnot((j+1)%5, 5, errors, errorRates)
    cnot((j+2)%5, 5, errors, errorRates)
    dualcz((j+3)%5, 5, errors, errorRates)
    syndromes[j] = measZ(5, errors, errorRates)
  return syndromes

# Implement the error correction procedure in Section II in the paper. For example, the circuit for measurement of XZZXI follows FIG.2 (c).
def correctErrors(errors, errorRates, verbose=False): 
  if verbose: print "starting syndrome0"
  prepZ(5, errors, errorRates)
  prepX(6, errors, errorRates)
  dualcz(0, 5, errors, errorRates)
  cnot(6, 5, errors, errorRates)
  cnot(1, 5, errors, errorRates)
  cnot(2, 5, errors, errorRates)
  cnot(6, 5, errors, errorRates)
  dualcz(3, 5, errors, errorRates)
  syndrome0 = measZ(5, errors, errorRates)
  flag0 = measX(6, errors, errorRates)
  if flag0: 
    if verbose: print "flag0"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    if syndromes == [0,0,0,1]: 
      errors.x ^= 1<<0
    elif syndromes == [0,1,0,0]: 
      errors.x ^= 1<<3
      errors.z ^= 1<<2
    elif syndromes == [0,1,1,0]: 
      errors.x ^= 1<<3
    elif syndromes == [1,0,0,0]: 
      errors.x ^= (1<<2) ^ (1<<3)
      errors.z ^= 1<<2
    elif syndromes == [1,0,0,1]: 
      errors.x ^= (1<<0) ^ (1<<1)
    elif syndromes == [1,0,1,0]: 
      errors.x ^= (1<<2) ^ (1<<3)
    elif syndromes == [1,1,0,0]: 
      errors.x ^= (1<<3) + (1<<4)
      errors.z ^= 1<<3
    return 1
  elif syndrome0: 
    if verbose: print "syndrome0"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    correctErrorsUsingSyndromes(errors, syndromes)
    return 1
  
  if verbose: print "starting syndrome1"
  prepZ(5, errors, errorRates)
  prepX(6, errors, errorRates)
  dualcz(1, 5, errors, errorRates)
  cnot(6, 5, errors, errorRates)
  cnot(2, 5, errors, errorRates)
  cnot(3, 5, errors, errorRates)
  cnot(6, 5, errors, errorRates)
  dualcz(4, 5, errors, errorRates)
  syndrome1 = measZ(5, errors, errorRates)
  flag1 = measX(6, errors, errorRates)
  if flag1: 
    if verbose: print "flag1"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    if syndromes == [0,0,1,1]: 
      errors.x ^= 1<<4
    elif syndromes == [0,1,0,0]: 
      errors.x ^= (1<<1) ^ (1<<2)
    elif syndromes == [0,1,0,1]: 
      errors.x ^= (1<<3) ^ (1<<4)
    elif syndromes == [0,1,1,0]: 
      errors.x ^= (1<<0) ^ (1<<4)
      errors.z ^= 1<<4
    elif syndromes == [1,0,0,0]: 
      errors.x ^= 1<<1
    elif syndromes == [1,0,1,0]: 
      errors.x ^= 1<<4
      errors.z ^= 1<<3
    elif syndromes == [1,1,0,0]: 
      errors.x ^= (1<<3) + (1<<4)
      errors.z ^= 1<<3
    return 1
  elif syndrome1: 
    if verbose: print "syndrome1"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    correctErrorsUsingSyndromes(errors, syndromes)
    return 1
  
  if verbose: print "starting syndrome2"
  prepZ(5, errors, errorRates)
  prepX(6, errors, errorRates)
  dualcz(2, 5, errors, errorRates)
  cnot(6, 5, errors, errorRates)
  cnot(3, 5, errors, errorRates)
  cnot(4, 5, errors, errorRates)
  cnot(6, 5, errors, errorRates)
  dualcz(0, 5, errors, errorRates)
  syndrome2 = measZ(5, errors, errorRates)
  flag2 = measX(6, errors, errorRates)
  if flag2: 
    if verbose: print "flag2"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    if syndromes == [0,0,0,1]: 
      errors.x ^= 1<<0
    elif syndromes == [0,0,1,0]: 
      errors.x ^= (1<<0) ^ (1<<4)
    elif syndromes == [0,0,1,1]: 
      errors.x ^= (1<<0) ^ (1<<1)
      errors.z ^= 1<<0
    elif syndromes == [0,1,0,1]: 
      errors.x ^= 1<<0
      errors.z ^= 1<<4
    elif syndromes == [0,1,1,0]: 
      errors.x ^= (1<<0) ^ (1<<4)
      errors.z ^= 1<<4
    elif syndromes == [1,0,1,0]: 
      errors.x ^= (1<<2) ^ (1<<3)
    elif syndromes == [1,1,0,0]: 
      errors.x ^= 1<<2
    return 1
  elif syndrome2: 
    if verbose: print "syndrome2"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    correctErrorsUsingSyndromes(errors, syndromes)
    return 1
  
  if verbose: print "starting syndrome3"
  prepZ(5, errors, errorRates)
  prepX(6, errors, errorRates)
  dualcz(3, 5, errors, errorRates)
  cnot(6, 5, errors, errorRates)
  cnot(4, 5, errors, errorRates)
  cnot(0, 5, errors, errorRates)
  cnot(6, 5, errors, errorRates)
  dualcz(1, 5, errors, errorRates)
  syndrome3 = measZ(5, errors, errorRates)
  flag3 = measX(6, errors, errorRates)
  if flag3: 
    if verbose: print "flag3"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    if syndromes == [0,0,0,1]: 
      errors.x ^= (1<<3) ^ (1<<4)
      errors.z ^= 1<<4
    elif syndromes == [0,0,1,0]: 
      errors.x ^= 1<<1
      errors.z ^= 1<<0
    elif syndromes == [0,0,1,1]: 
      errors.x ^= (1<<0) ^ (1<<1)
      errors.z ^= 1<<0
    elif syndromes == [0,1,0,1]: 
      errors.x ^= (1<<3) ^ (1<<4)
    elif syndromes == [0,1,1,0]: 
      errors.x ^= 1<<3
    elif syndromes == [1,0,0,0]: 
      errors.x ^= 1<<1
    elif syndromes == [1,0,0,1]: 
      errors.x ^= (1<<0) ^ (1<<1)
    return 1
  elif syndrome3: 
    if verbose: print "syndrome3"
    syndromes = extractSyndromes(errors, errorRates)
    if verbose: print syndromes
    correctErrorsUsingSyndromes(errors, syndromes)
    return 1
  return 0

# Find least weight representation modulo stabilizers.
def weight(errors):
  return bin((errors.x | errors.z) & ((1 << 5) - 1)).count("1")

def reduceError(errors): 
  stabilizers = [[(1<<0)+(1<<3),(1<<1)+(1<<2)], [(1<<1)+(1<<4),(1<<2)+(1<<3)], [(1<<2)+(1<<0),(1<<3)+(1<<4)], [(1<<3)+(1<<1),(1<<4)+(1<<0)]]
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
    if (errorsCopy.x & ((1<<5)-1)) or (errorsCopy.z & ((1<<5)-1)): 
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

