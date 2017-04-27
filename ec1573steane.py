# Implement the fault-tolerant error correction of [[15,7,3]] code using Steane's scheme. For the preparation of logical states, we adopt the simpler two-ancialla scheme in our paper.

from utility import *

# Prepare logical |0> and |+> states for syndrome measurement, following the circuit in Appendix C FIG. 9, which uses bare flag ancilla.
def prepLogical0(errors, errorRates, verbose):
  z = 1
  while(z):
    prepX(15, errors, errorRates)
    prepX(16, errors, errorRates)
    prepZ(17, errors, errorRates)
    prepX(18, errors, errorRates)
    prepZ(19, errors, errorRates)
    prepZ(20, errors, errorRates)
    prepZ(21, errors, errorRates)
    prepX(22, errors, errorRates)
    prepZ(23, errors, errorRates)
    prepZ(24, errors, errorRates)
    prepZ(25, errors, errorRates)
    prepZ(26, errors, errorRates)
    prepZ(27, errors, errorRates)
    prepZ(28, errors, errorRates)
    prepZ(29, errors, errorRates)
    prepZ(30, errors, errorRates)
    cnot(22, 30, errors, errorRates)
    cnot(18, 30, errors, errorRates)
    cnot(16, 30, errors, errorRates)
    cnot(15, 30, errors, errorRates)
    for i in range(4):
      for j in range(7):
        cnot(15+cnotWires[i][0], 15+cnotWires[i][j+1], errors, errorRates)
    cnot(15, 30, errors, errorRates)
    cnot(16, 30, errors, errorRates)
    cnot(18, 30, errors, errorRates)
    cnot(22, 30, errors, errorRates)
    z = measZ(30, errors, errorRates)
    if verbose&z:
      print "ancilla zero fail"

def prepLogicalPlus(errors, errorRates, verbose):
  z = 1
  while(z):
    prepZ(15, errors, errorRates)
    prepZ(16, errors, errorRates)
    prepX(17, errors, errorRates)
    prepZ(18, errors, errorRates)
    prepX(19, errors, errorRates)
    prepX(20, errors, errorRates)
    prepX(21, errors, errorRates)
    prepZ(22, errors, errorRates)
    prepX(23, errors, errorRates)
    prepX(24, errors, errorRates)
    prepX(25, errors, errorRates)
    prepX(26, errors, errorRates)
    prepX(27, errors, errorRates)
    prepX(28, errors, errorRates)
    prepX(29, errors, errorRates)
    prepX(30, errors, errorRates)
    cnot(30, 22, errors, errorRates)
    cnot(30, 18, errors, errorRates)
    cnot(30, 16, errors, errorRates)
    cnot(30, 15, errors, errorRates)
    for i in range(4):
      for j in range(7):
        cnot(15+cnotWires[i][j+1], 15+cnotWires[i][0], errors, errorRates)
    cnot(30, 15, errors, errorRates)
    cnot(30, 16, errors, errorRates)
    cnot(30, 18, errors, errorRates)
    cnot(30, 22, errors, errorRates)
    z = measX(30, errors, errorRates)
    if verbose&z:
      print "ancilla zero fail"

# Perform transversal cnot between data block and logical |0> and |+>. Correct the error through parity checks.
def correctErrors(errors, errorRates, verbose=False):
  if verbose: print "starting Xsyndrome"
  prepLogical0(errors, errorRates, verbose)
  Xsyndrome = 0
  for i in range(15):
    cnot(15+i, i, errors, errorRates)
    if measX(15+i, errors, errorRates):
      Xsyndrome ^= (1+i)
  if Xsyndrome:
    errors.z ^= 1<<(Xsyndrome-1)
  if verbose: print "starting Zsyndrome"
  prepLogicalPlus(errors, errorRates, verbose)
  Zsyndrome = 0
  for i in range(15):
    cnot(i, 15+i, errors, errorRates)
    if measZ(15+i, errors, errorRates):
      Zsyndrome ^= (1+i)
  if Zsyndrome:
    errors.x ^= 1<<(Zsyndrome-1)


# Find least weight representation modulo stabilizers.
def weight(errors):
  return bin((errors.x | errors.z) & ((1 << 15) - 1)).count("1")

# This is the matrix indexing the non-identity positions of the stabilizers 
cnotWires = [[k for k in range(15) if ((k+1)>>(3-i))&1==1] for i in range(4)]

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
# for i in range(10):
#   print "gamma=10^(%d/10-4), trials=10^7"% i
#   simulateErrorCorrection(gammas[i], 10**7)
for i in range(11):
  print "gamma=10^(%d/10-4), trials=10^6"% (i+10)
  simulateErrorCorrection(gammas[i+10], 10**6)

