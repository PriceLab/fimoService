from FimoServer import *

if len(sys.argv) != 4:
   print("required args:  port   path.to.fimo.executable   path.to.motifs.file")
   sys.exit(0)

hostname = "localhost"
port = sys.argv[1]
fimoExecutable = sys.argv[2]  #  "/Users/paul/meme/bin/fimo"
motifsFile = sys.argv[3]       # = "/Users/paul/s/work/priceLab/cory/footprintWorkflow/JASPAR_CORE_plus_seth.meme"

server = FimoServer(hostname, port, fimoExecutable, motifsFile)

while True:
    server.handleRequest()
