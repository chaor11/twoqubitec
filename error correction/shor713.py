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
    return "x= %s, z= %s" % (etostr(self.x, 7), etostr(self.z, 7))

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
  xsyndrome = (syndromes[0]<<2) + (syndromes[1]<<1) + syndromes[2]
  if xsyndrome:
    errors.z ^= 1<<(xsyndrome-1)
  zsyndrome = (syndromes[3]<<2) + (syndromes[4]<<1) + syndromes[5]
  if zsyndrome:
    errors.x ^= 1<<(zsyndrome-1)

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

def extractXSyndromes(errors, errorRates): 
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
  return syndromes

def extractZSyndromes(errors, errorRates): 
  syndromes = [0 for i in range(6)]
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

def prepCat(errors, errorRates, verbose):
  z = 1
  while(z):
    prepX(7, errors, errorRates)
    prepZ(8, errors, errorRates)
    prepZ(9, errors, errorRates)
    prepZ(10, errors, errorRates)
    cnot(7, 8, errors, errorRates)
    cnot(7, 9, errors, errorRates)
    cnot(9, 10, errors, errorRates)
    prepZ(11, errors, errorRates)
    cnot(8, 11, errors, errorRates)
    cnot(10, 11, errors, errorRates)
    z = measZ(11, errors, errorRates)
    if z&verbose: print "cat fail"

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

def weight(errors, n=7):
  return bin((errors.x | errors.z) & ((1 << n) - 1)).count("1")

def reduceError(errors): 
  stabilizers = \
  [[(1<<6)+(1<<5)+(1<<4)+(1<<3),0], \
  [(1<<6)+(1<<5)+(1<<2)+(1<<1),0], \
  [(1<<6)+(1<<4)+(1<<2)+(1<<0),0], \
  [0,(1<<6)+(1<<5)+(1<<4)+(1<<3)], \
  [0,(1<<6)+(1<<5)+(1<<2)+(1<<1)], \
  [0,(1<<6)+(1<<4)+(1<<2)+(1<<0)], \
  ]
  bestErrors = Errors()
  bestErrors.x = errors.x
  bestErrors.z = errors.z
  bestWeight = weight(bestErrors,7)
  trialErrors = Errors()
  for k in range(1, 1<<(len(stabilizers))):
    trialErrors.x = errors.x
    trialErrors.z = errors.z
    for digit in range(len(stabilizers)):
      if (k>>digit)&1: 
        trialErrors.x ^= stabilizers[digit][0]
        trialErrors.z ^= stabilizers[digit][1]
    if weight(trialErrors,7) < bestWeight: 
      bestErrors.x = trialErrors.x
      bestErrors.z = trialErrors.z
      bestWeight = weight(bestErrors,7)
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
    if (errorsCopy.x & ((1<<7)-1)) or (errorsCopy.z & ((1<<7)-1)): 
      failures += 1
      errors.x = 0
      errors.z = 0
  print failures

gammas = [10**(i/10.-4) for i in range(21)]
print "shor713"
for i in range(10):
  print "gamma=10^(%d/10-4), trials=10^6"% i
  simulateErrorCorrection(gammas[i], 10**6)
for i in range(11):
  print "gamma=10^(%d/10-4), trials=10^6"% (i+10)
  simulateErrorCorrection(gammas[i+10], 10**6)

