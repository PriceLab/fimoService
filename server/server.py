import zmq
import time
import sys
import json
import subprocess
import os
import os.path
import tempfile
import pandas
#------------------------------------------------------------------------------------------------------------------------
def runTests():

   test_writeSequencesToTemporaryFastaFile()
   test_runFimo()

#------------------------------------------------------------------------------------------------------------------------
def runFimo(sequences, fimoExecutable, motifsFile):

   sequencesFile = writeSequencesToTemporaryFastaFile(sequences);
   outputDirectory = tempfile.mktemp()
   outputDirectorySwitch = "--oc %s" % outputDirectory
   args = [fimoExecutable, "--oc", outputDirectory, motifsFile, sequencesFile]
   devnull = open(os.devnull, 'w')
   processStatus = subprocess.check_call(args, stdout=devnull, stderr=devnull)

   tbl = pandas.DataFrame()

   if(processStatus == 0):
     filename = "%s/%s" % (outputDirectory, "fimo.txt")
     tbl = pandas.read_csv(filename, delimiter="\t")

   return(tbl)

#------------------------------------------------------------------------------------------------------------------------
def writeSequencesToTemporaryFastaFile(sequences):

   f = tempfile.NamedTemporaryFile(mode='w', suffix=".fa", delete=False)
   filename = f.name

   for key in sequences.keys():
      sequence = sequences[key]
      f.write("> %s\n" % key)
      f.write("%s\n" % sequence)

   f.close()

   return(filename)

#------------------------------------------------------------------------------------------------------------------------
def sampleSequences():

   sequences = {"tert_wt1": "CCCGGAGGGGG", "tert_wt2": "CCCGGGAGGGG", "tert_mut": "CCCCTTCCGGG"}
   return(sequences)

#------------------------------------------------------------------------------------------------------------------------
def test_writeSequencesToTemporaryFastaFile():

   print("--- test_writeSequencesToTemporaryFastaFile")

   filename = writeSequencesToTemporaryFastaFile(sampleSequences())
   assert(os.path.exists(filename))
   lines = open(filename).read().split("\n")
   assert(len(lines) == 7)   # an empty token on the end

#------------------------------------------------------------------------------------------------------------------------
def test_runFimo():

   print("--- test_runFimo")

   tbl = runFimo(sampleSequences())
   assert(tbl.shape == (4,9))

#------------------------------------------------------------------------------------------------------------------------
if len(sys.argv) != 4:
   print("required args:  port   path.to.fimo.executable   path.to.motifs.file")
   sys.exit(0)

port = sys.argv[1]
fimoExecutable = sys.argv[2]  #  "/Users/paul/meme/bin/fimo"
assert(os.path.exists(fimoExecutable))
motifsFile = sys.argv[3]       # = "/Users/paul/s/work/priceLab/cory/footprintWorkflow/JASPAR_CORE_plus_seth.meme"
assert(os.path.exists(motifsFile))

context = zmq.Context()
socket = context.socket(zmq.REP)
socket.bind("tcp://*:%s" % port)

while True:
    request = json.loads(socket.recv_string())
    sequences = request['sequences']
    print("calling runFimo on %d sequences" % len(sequences))
    tbl = runFimo(sequences, fimoExecutable, motifsFile)
    obj = pandas.DataFrame.to_json(tbl)
    socket.send_string(obj)
