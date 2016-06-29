import zmq
import time
import sys
import json
import subprocess
import os
import os.path
import tempfile
import shutil
import pandas
#------------------------------------------------------------------------------------------------------------------------
class FimoServer:
   "run the Meme Suite's Fimo tool in response to a ZeroMQ message"

   host = "localhost"
   port = 6123
   fimoExecutable = ""
   motifsVile = ""
   socketContext = None
   socket = None

   #--------------------------------------------------------------------------------
   def __init__(self, host="localhost", port=6123, fimoExecutable=None, motifsFile=None):

      self.host = host
      self.port = int(port)
      assert(os.path.exists(fimoExecutable))
      self.fimoExecutable = fimoExecutable
      assert(os.path.exists(motifsFile))
      self.motifsFile = motifsFile
      self.socketContext = zmq.Context()
      self.socket = self.socketContext.socket(zmq.REP)
      self.socket.bind("tcp://*:%s" % self.port)
      print("started FimoServer on port %d" % self.port)


   #--------------------------------------------------------------------------------
   def getHost(self):
      return(self.host)

   #--------------------------------------------------------------------------------
   def getPort(self):
      return(self.port)


   #--------------------------------------------------------------------------------
   def runFimo(self, sequences):

      sequencesFile = self.writeSequencesToTemporaryFastaFile(sequences);
      outputDirectory = tempfile.mktemp()
      outputDirectorySwitch = "--oc %s" % outputDirectory
      print("writing fimo files to %s" % outputDirectory)
      args = [self.fimoExecutable, "--oc", outputDirectory, self.motifsFile, sequencesFile]
      devnull = open(os.devnull, 'w')
      processStatus = subprocess.check_call(args, stdout=devnull, stderr=devnull)
      tbl = pandas.DataFrame()
      if(processStatus == 0):
         filename = "%s/%s" % (outputDirectory, "fimo.txt")
         tbl = pandas.read_csv(filename, delimiter="\t")
      shutil.rmtree(outputDirectory)
      print("deleted %s" % outputDirectory)
      return(tbl)

   #--------------------------------------------------------------------------------
   def writeSequencesToTemporaryFastaFile(self, sequences):

      f = tempfile.NamedTemporaryFile(mode='w', suffix=".fa", delete=False)
      filename = f.name

      for key in sequences.keys():
         sequence = sequences[key]
         f.write("> %s\n" % key)
         f.write("%s\n" % sequence)

      f.close()

      return(filename)


   #--------------------------------------------------------------------------------
   def handleRequest(self):

      request = json.loads(self.socket.recv_string())
      sequences = request;
      #sequences = request['sequences']
      print("%s  calling runFimo on %d sequences" % \
            (time.strftime("%Y-%m-%d %H:%M:%S"), len(sequences)))
      tbl = self.runFimo(sequences)
      obj = pandas.DataFrame.to_json(tbl)
      self.socket.send_string(obj)


#-----------------------------------------------------------------------------------


