# Implement the preparation and measurement of logical X states for [[5,1,3]] code. Each module uses only two ancilla qubits. The circuits follow Appendix B FIG. 8.

from utility import *

# Follow the circuit in FIG. 8 (a). The circuit actually prepare logical |->. Can perform transversal Pauli to get |+>, in which case the error is unchanged.
def prepLogicalX(errors, errorRates, verbose=False):
  z=1
  while(z):
    for i in range(5):
      prepX(i, errors, errorRates)
    z = 0
    cz(0, 1, errors, errorRates)
    cz(2, 3, errors, errorRates)
    cz(1, 2, errors, errorRates)
    cz(3, 4, errors, errorRates)
    cz(0, 4, errors, errorRates)
    prepX(5, errors, errorRates)
    prepZ(6, errors, errorRates)
    cnot(5, 0, errors, errorRates)
    cnot(5, 6, errors, errorRates)
    cz(5, 1, errors, errorRates)
    cnot(5, 6, errors, errorRates)
    cz(5, 4, errors, errorRates)
    z += measX(5, errors, errorRates)
    z += measZ(6, errors, errorRates)
    prepX(5, errors, errorRates)
    prepZ(6, errors, errorRates)
    cz(5, 0, errors, errorRates)
    cnot(5, 6, errors, errorRates)
    cnot(5, 1, errors, errorRates)
    cnot(5, 6, errors, errorRates)
    cz(5, 2, errors, errorRates)
    z += measX(5, errors, errorRates)
    z += measZ(6, errors, errorRates)
    prepX(5, errors, errorRates)
    prepZ(6, errors, errorRates)
    cnot(5, 3, errors, errorRates)
    cnot(5, 6, errors, errorRates)
    cz(5, 2, errors, errorRates)
    cnot(5, 6, errors, errorRates)
    cz(5, 4, errors, errorRates)
    z += measZ(6, errors, errorRates)
    z += measX(5, errors, errorRates)

# The decoder in FIG. 8 (b) Part IV.
def decodeLogicalX(errors, rates, verbose=False): 
  cz(0, 4, errors, rates)
  cz(3, 4, errors, rates)
  cz(1, 2, errors, rates)
  cz(2, 3, errors, rates)
  cz(0, 1, errors, rates)
  measurements = [measX(j, errors, rates) for j in range(5)]
  if measurements == [0] * 5 or measurements == [1] * 5: 
    return measurements[0]
  for qubit in range(5): 
    if measurements == [0,1,0,0,1] or measurements == [1,0,1,1,0]: # X error on first qubit
      return measurements[0]
    if measurements == [1,1,0,0,1] or measurements == [0,0,1,1,0]: # Y error on first qubit
      return measurements[2]
    if measurements == [1,0,0,0,0] or measurements == [0,1,1,1,1]: # Z error on first qubit
      return measurements[1]
    measurements = measurements[1:] + measurements[:1] 
  return measurements[0]  

# Non-deterministic procedure measuring logical X. Follow the circuit in FIG. 8 (b).
def measureLogicalX(errors, rates, verbose=False): 
  prepX(5,   errors, rates)
  prepZ(6,   errors, rates)
  cz(  5, 4, errors, rates)
  cnot(5, 6, errors, rates)
  cnot(5, 0, errors, rates)
  cnot(5, 6, errors, rates)
  cz(  5, 1, errors, rates)
  s1 = measX(5, errors, rates)
  f1 = measZ(6, errors, rates)
  if f1: 
    cz(0, 4, errors, rates)
    cz(3, 4, errors, rates)
    cz(1, 2, errors, rates)
    cz(2, 3, errors, rates)
    cz(0, 1, errors, rates)
    measurements = [measX(j, errors, rates) for j in range(5)]
    if measurements == [0] * 5 or measurements == [1] * 5: # error is IIIII
      return measurements[0]
    if measurements == [0,1,0,0,0] or measurements == [1,0,1,1,1]: # error is IZIII
      return measurements[0]
    if measurements == [0,0,0,0,1] or measurements == [1,1,1,1,0]: # error is XIIIZ
      return measurements[0]
    if measurements == [1,0,0,0,1] or measurements == [0,1,1,1,0]: # error is YIIIZ
      return measurements[1]
    if measurements == [1,1,0,0,0] or measurements == [0,0,1,1,1]: # error is ZZIII
      return measurements[2]
    return measurements[0]  # this should never happen
  prepX(5,   errors, rates)
  prepZ(6,   errors, rates)
  cz(  5, 0, errors, rates)
  cnot(5, 6, errors, rates)
  cnot(5, 1, errors, rates)
  cnot(5, 6, errors, rates)
  cz(  5, 2, errors, rates)
  s2 = measX(5, errors, rates)
  f2 = measZ(6, errors, rates)
  if f2: 
    return s1
  if s2 != s1: 
    return decodeLogicalX(errors, rates, verbose)
  prepX(5,   errors, rates)
  prepZ(6,   errors, rates)
  cz(  5, 2, errors, rates)
  cnot(5, 6, errors, rates)
  cnot(5, 3, errors, rates)
  cnot(5, 6, errors, rates)
  cz(  5, 4, errors, rates)
  s4 = measX(5, errors, rates)
  f4 = measZ(6, errors, rates)
  if f4: 
    return s1
  if s4 == s2: 
    return s4
  # flag is not triggered, s4 != s1 = s2
  return decodeLogicalX(errors, rates, verbose)

# Run independent trials of preparation and measurement of logical X states, with pysical error rate being gamma. 
# Count the number of wrong measurements and calculate its ratio over trials.
def simulatePrepMeasLogicalX(gamma, trials):
  errors = Errors(0, 0)
  errorRates = ErrorRates((4/15.)*gamma, gamma, (4/15.)*gamma)
  
  failures = 0
  for k in xrange(trials): 
    errors.x = 0
    errors.z = 0
    prepLogicalX(errors, errorRates)
    failures += measureLogicalX(errors, errorRates)
  print failures

# Wrapper function for the plot. More trials are needed for small gammas due to the confidence interval.
gammas = [10**(i/10.-4) for i in range(21)]
# for i in range(10):
#   print "gamma=10^(%d/10-4), trials=10^7"% i
#   simulatePrepMeasLogicalX(gammas[i], 10**7)
for i in range(11):
  print "gamma=10^(%d/10-4), trials=10^6"% (i+10)
  simulatePrepMeasLogicalX(gammas[i+10], 10**6)


