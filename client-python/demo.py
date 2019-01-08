import zmq
import sys
import json
import pandas
#------------------------------------------------------------------------------------------------------------------------
if len(sys.argv) != 3:
   print("required args:  hostname port")
   sys.exit(0)

pandas.set_option('display.width', 1000)
hostname = sys.argv[1]
port = sys.argv[2]

context = zmq.Context()
print("Connecting to server...")
socket = context.socket(zmq.REQ)
socket.connect("tcp://%s:%s" % (hostname, port))

sequences = {"tert_wt1": "CCCGGAGGGGG", "tert_wt2": "CCCGGGAGGGG", "tert_mut": "CCCCTTCCGGG"}

msg = {"sequences": sequences, "pvalThreshold": 0.0001}

for request in range (1,3):
    print("Sending request ", request,"...")
    socket.send_string(json.dumps(msg))
    responseString = socket.recv_string()
    #print(responseString)
    tbl = pandas.read_json(responseString)
    print("Received table with dimensions: ",  tbl.shape)
    print(tbl)
