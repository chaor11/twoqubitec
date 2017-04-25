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
    return "x= %s, z= %s" % (etostr(self.x, 5), etostr(self.z, 5))

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

def weight(errors, n=5):
  return bin((errors.x | errors.z) & ((1 << n) - 1)).count("1")

def reduceError(errors): 
  stabilizers = [[(1<<0)+(1<<3),(1<<1)+(1<<2)], [(1<<1)+(1<<4),(1<<2)+(1<<3)], [(1<<2)+(1<<0),(1<<3)+(1<<4)], [(1<<3)+(1<<1),(1<<4)+(1<<0)]]
  bestErrors = Errors()
  bestErrors.x = errors.x
  bestErrors.z = errors.z
  bestWeight = weight(bestErrors)
  trialErrors = Errors()
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
    if (errorsCopy.x & ((1<<5)-1)) or (errorsCopy.z & ((1<<5)-1)): 
      failures += 1
      errors.x = 0
      errors.z = 0
  print failures

gammas = [10**(i/10.-4) for i in range(21)]
print "flag513"
for i in range(10):
  print "gamma=10^(%d/10-4), trials=10^7"% i
  simulateErrorCorrection(gammas[i], 10**7)
for i in range(11):
  print "gamma=10^(%d/10-4), trials=10^6"% (i+10)
  simulateErrorCorrection(gammas[i+10], 10**6)

