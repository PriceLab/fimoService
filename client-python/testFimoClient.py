from FimoClient import *
from datetime import datetime
#------------------------------------------------------------------------------------------------------------------------
if len(sys.argv) != 3:
   print("required args:  hostname port")
   sys.exit(0)

hostname = sys.argv[1]
port = sys.argv[2]

fimoClient = FimoClient(hostname,port)
print(fimoClient.getHost())
sequences = {"tert_wt1": "CCCGGAGGGGG", "tert_wt2": "CCCGGGAGGGG", "tert_mut": "CCCCTTCCGGG"}

for i in range (0, 100):
   start = datetime.now()
   tbl = fimoClient.request(sequences)
   end = datetime.now()
   delta = end - start
   msecs = delta.seconds + (delta.microseconds / 1000)
   print("testFimoClient.py, loop %4d: %6f.1" % (i, msecs))
   assert(str(type(tbl)) == "<class 'pandas.core.frame.DataFrame'>")
   assert(tbl.shape == (4,9))
   assert(list(tbl.columns.values) == ['#pattern name', 'matched sequence', 'p-value', 'q-value', 'score',
                                       'sequence name', 'start', 'stop', 'strand'])
