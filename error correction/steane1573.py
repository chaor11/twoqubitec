import random

def etostr(e, n):
  bits = ['0'] * n
  for i in range(n):
    if (e>>i)&1: 
      bits[i] = '1'
  return ''.join(bits)

class ErrorRates():
  def __str__(self):
    return "cnot= %.6f, prep= %.6f, meas= %.6f" \
      % (self.cnot, self.prep, self.meas)

class Errors(object):
  def __str__(self): 
    return "x= %s, z= %s" % (etostr(self.x, 15), etostr(self.z, 15))

def errorsX1(bit, errors, rate): 
  if random.random() < rate: 
    errors.x ^= 1<<bit

def errorsZ1(bit, errors, rate): 
  if random.random() < rate: 
    errors.z ^= 1<<bit

def errors1(bit, errors, rate): 
  r = random.random()
  if r < rate: 
    r = 1 + int(3 * r / rate)
    if r&1: 
      errors.x ^= 1<<bit
    if r>>1&1: 
      errors.z ^= 1<<bit

def errors2(bit1, bit2, errors, rate): 
  r = random.random()    
  if r < rate: 
    r = 1 + int(15 * r / rate)  
    if r&1: 
      errors.x ^= 1<<bit1
    if r>>1&1:
      errors.z ^= 1<<bit1
    if r>>2&1:
      errors.x ^= 1<<bit2
    if r>>3&1:
      errors.z ^= 1<<bit2

def cnot0(bit1, bit2, errors): 
  errors.x ^= (errors.x>>bit1&1)<<bit2
  errors.z ^= (errors.z>>bit2&1)<<bit1

def hadamard0(bit, errors): 
  errxz = (errors.x & (1<<bit)) ^ (errors.z & (1<<bit))
  errors.x ^= errxz
  errors.z ^= errxz

def cz0(bit1, bit2, errors): 
  hadamard0(bit2, errors)
  cnot0(bit1, bit2, errors)
  hadamard0(bit2, errors)

def dualcz0(bit1, bit2, errors): 
  hadamard0(bit1, errors)
  cnot0(bit1, bit2, errors)
  hadamard0(bit1, errors)

def prep0(bit, errors): 
  errors.x ^= errors.x & (1<<bit)
  errors.z ^= errors.z & (1<<bit)

def prepZ(bit, errors, errorRates): 
  prep0(bit, errors)
  errorsX1(bit, errors, errorRates.prep)

def prepX(bit, errors, errorRates): 
  prep0(bit, errors)
  errorsZ1(bit, errors, errorRates.prep)

def measZ(bit, errors, errorRates): 
  errorsX1(bit, errors, errorRates.meas)
  return (errors.x>>bit)&1

def measX(bit, errors, errorRates): 
  errorsZ1(bit, errors, errorRates.meas)
  return (errors.z>>bit)&1

def cnot(bit1, bit2, errors, errorRates): 
  cnot0(bit1, bit2, errors)  
  errors2(bit1, bit2, errors, errorRates.cnot)

def cz(bit1, bit2, errors, errorRates): 
  cz0(bit1, bit2, errors)
  errors2(bit1, bit2, errors, errorRates.cnot)

def dualcz(bit1, bit2, errors, errorRates): 
  dualcz0(bit1, bit2, errors)
  errors2(bit1, bit2, errors, errorRates.cnot)

def prep(bit, errors, errorRates): 
  prep0(bit, errors)
  errors1(bit, errors, errorRates.prep)

def correctErrorsUsingSyndromes(errors, syndromes):
  xsyndrome = (syndromes[0]<<3) + (syndromes[1]<<2) + (syndromes[2]<<1) + syndromes[3]
  if xsyndrome:
    errors.z ^= 1<<(xsyndrome-1)
  zsyndrome = (syndromes[4]<<3) + (syndromes[5]<<2) + (syndromes[6]<<1) + syndromes[7]
  if zsyndrome:
    errors.x ^= 1<<(zsyndrome-1)

# This is the matrix indexing the non-identity positions of the stabilizers 
cnotWires = [[k for k in range(15) if ((k+1)>>(3-i))&1==1] for i in range(4)]

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

def prepLogic0(errors, errorRates, verbose):
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

def prepLogicPlus(errors, errorRates, verbose):
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

def correctErrors(errors, errorRates, verbose=False):
  if verbose: print "starting Xsyndrome"
  prepLogic0(errors, errorRates, verbose)
  Xsyndrome = 0
  for i in range(15):
    cnot(15+i, i, errors, errorRates)
    if measX(15+i, errors, errorRates):
      Xsyndrome ^= (1+i)
  if Xsyndrome:
    errors.z ^= 1<<(Xsyndrome-1)
  if verbose: print "starting Zsyndrome"
  prepLogicPlus(errors, errorRates, verbose)
  Zsyndrome = 0
  for i in range(15):
    cnot(i, 15+i, errors, errorRates)
    if measZ(15+i, errors, errorRates):
      Zsyndrome ^= (1+i)
  if Zsyndrome:
    errors.x ^= 1<<(Zsyndrome-1)

def weight(errors, n=15):
  return bin((errors.x | errors.z) & ((1 << n) - 1)).count("1")

stabilizers = [[0,0] for i in range(8)]
for i in range(4):
  for j in range(8):
    stabilizers[i][0] += 1<<cnotWires[i][j]
for i in range(4):
  for j in range(8):
    stabilizers[i+4][1] += 1<<cnotWires[i][j]

def reduceError(errors): 
  bestErrors = Errors()
  bestErrors.x = errors.x
  bestErrors.z = errors.z
  bestWeight = weight(bestErrors,15)
  trialErrors = Errors()
  for k in range(1, 1<<(len(stabilizers))):
    trialErrors.x = errors.x
    trialErrors.z = errors.z
    for digit in range(len(stabilizers)):
      if (k>>digit)&1==1: 
        trialErrors.x ^= stabilizers[digit][0]
        trialErrors.z ^= stabilizers[digit][1]
    if weight(trialErrors,15) < bestWeight: 
      bestErrors.x = trialErrors.x
      bestErrors.z = trialErrors.z
      bestWeight = weight(bestErrors,15)
  return bestErrors

def simulateErrorCorrection(gamma, trials): 
  errors = Errors()
  errorsCopy = Errors()
  
  errorRates0 = ErrorRates()
  errorRates0.cnot = 0
  errorRates0.prep = 0
  errorRates0.meas = 0

  errorRates = ErrorRates()
  errorRates.cnot = gamma
  errorRates.prep = (4/15.) * errorRates.cnot
  errorRates.meas = (4/15.) * errorRates.cnot
  
  failures = 0
  errors.x = 0
  errors.z = 0
  for k in xrange(trials): 
    correctErrors(errors, errorRates, verbose=False)
    errorsCopy.x = errors.x
    errorsCopy.z = errors.z
    correctErrors(errorsCopy, errorRates0, verbose=False)
    errorsCopy = reduceError(errorsCopy)
    if (errorsCopy.x & ((1<<15)-1)) or (errorsCopy.z & ((1<<15)-1)): 
      failures += 1
      errors.x = 0
      errors.z = 0
  print failures

gammas = [10**(i/10.-4) for i in range(21)]
print "steane1573"
for i in range(10):
  print "gamma=10^(%d/10-4), trials=10^6"% i
  simulateErrorCorrection(gammas[i], 10**6)
for i in range(11):
  print "gamma=10^(%d/10-4), trials=10^6"% (i+10)
  simulateErrorCorrection(gammas[i+10], 10**6)

  