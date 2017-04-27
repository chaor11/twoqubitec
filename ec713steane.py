# Implement the fault-tolerant error correction of [[7,1,3]] code using Steane's scheme. For the preparation of logical states, we adopt the simpler two-ancialla scheme in our paper.

from utility import *

# Prepare logical |0> and |+> states for syndrome measurement, following the circuit in Appendix C FIG. 9, which uses bare flag ancilla.
def prepLogical0(errors, errorRates, verbose):
  z = 1
  while(z):
    prepX(7, errors, errorRates)
    prepX(8, errors, errorRates)
    prepX(10, errors, errorRates)
    prepZ(9, errors, errorRates)
    prepZ(11, errors, errorRates)
    prepZ(12, errors, errorRates)
    prepZ(13, errors, errorRates)
    prepZ(14, errors, errorRates)
    cnot(7, 14, errors, errorRates)
    cnot(8,  14, errors, errorRates)
    cnot(10, 14, errors, errorRates)
    cnot(7, 9, errors, errorRates)
    cnot(7, 11, errors, errorRates)
    cnot(7, 13, errors, errorRates)
    cnot(8, 9, errors, errorRates)
    cnot(8, 12, errors, errorRates)
    cnot(8, 13, errors, errorRates)
    cnot(10, 11, errors, errorRates)
    cnot(10, 12, errors, errorRates)
    cnot(10, 13, errors, errorRates)
    cnot(10, 14, errors, errorRates)
    cnot(8, 14, errors, errorRates)
    cnot(7, 14, errors, errorRates)
    z = measZ(14, errors, errorRates)
    if verbose&z:
        print "ancilla zero fail"

def prepLogicalPlus(errors, errorRates, verbose):
  z = 1
  while(z):
    prepZ(7, errors, errorRates)
    prepZ(8, errors, errorRates)
    prepZ(10, errors, errorRates)
    prepX(9, errors, errorRates)
    prepX(11, errors, errorRates)
    prepX(12, errors, errorRates)
    prepX(13, errors, errorRates)
    prepX(14, errors, errorRates)
    cnot(14, 7, errors, errorRates)
    cnot(14, 8, errors, errorRates)
    cnot(14, 10, errors, errorRates)
    cnot(9, 7, errors, errorRates)
    cnot(11, 7, errors, errorRates)
    cnot(13, 7, errors, errorRates)
    cnot(9, 8, errors, errorRates)
    cnot(12, 8, errors, errorRates)
    cnot(13, 8, errors, errorRates)
    cnot(11, 10, errors, errorRates)
    cnot(12, 10, errors, errorRates)
    cnot(13, 10, errors, errorRates)
    cnot(14, 10, errors, errorRates)
    cnot(14, 8, errors, errorRates)
    cnot(14, 7, errors, errorRates)
    z = measX(14, errors, errorRates)
    if verbose&z:
        print "ancilla zero fail"

# Perform transversal cnot between data block and logical |0> and |+>. Correct the error through parity checks.
def correctErrors(errors, errorRates, verbose=False):
  if verbose: print "starting Xsyndrome"
  prepLogical0(errors, errorRates, verbose)
  Xsyndrome = 0
  for i in range(7):
    cnot(7+i, i, errors, errorRates)
    if measX(7+i, errors, errorRates):
      Xsyndrome ^= (1+i)
  if Xsyndrome:
    errors.z ^= 1<<(Xsyndrome-1)
  if verbose: print "starting Zsyndrome"
  prepLogicalPlus(errors, errorRates, verbose)
  Zsyndrome = 0
  for i in range(7):
    cnot(i, 7+i, errors, errorRates)
    if measZ(7+i, errors, errorRates):
      Zsyndrome ^= (1+i)
  if Zsyndrome:
    errors.x ^= 1<<(Zsyndrome-1)

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

