import zmq
import sys
import json
import pandas
#------------------------------------------------------------------------------------------------------------------------
class FimoClient:
   'Access to motif/sequence match on a (remote) server'

   host = None
   port = None
   context = None
   socket = None

   def __init__(self, host, port, quiet=True):
      self.host = host
      self.port = int(port)
      self.context = zmq.Context()
      self.quiet = quiet
      if(not self.quiet):
         print("Connecting to server %s:%s" % (self.host, self.port))
      self.socket = self.context.socket(zmq.REQ)
      self.socket.connect("tcp://%s:%s" % (host, port))

   def getHost(self):
      return('%s:%d' % (self.host, self.port))

   def request(self, sequences):
      if(not self.quiet):
         print("Sending request, length %d" % len(sequences))
      self.socket.send_string(json.dumps(sequences))
      #self.socket.send_string(json.dumps({'sequences': sequences}))
      responseString = self.socket.recv_string()
      if(not self.quiet):
         print(responseString)
      tbl = pandas.read_json(responseString)
      if(not self.quiet):
         print("Received table with dimensions: ",  tbl.shape)
      return(tbl)

