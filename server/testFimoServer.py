from FimoServer import *
import unittest

#hostname = None
#port = None

class TestFimoServer(unittest.TestCase):

  def testConstructor(self):
     print("--- testConstructor")
     server = FimoServer(hostname, port)
     self.assertEqual(server.getHost(), hostname)
     self.assertEqual(server.getPort(), port)

  def testRequest(self):
     print("--- testRequest")
     server = FimoServer(hostname, port)
     self.assertEqual(server.getHost(), hostname)
     self.assertEqual(server.getPort(), port)


if __name__ == '__main__':
   if len(sys.argv) != 3:
      print("required args:  port   path.to.fimo.executable   path.to.motifs.file")
      sys.exit(0)

   hostname = sys.argv[1]
   port = int(sys.argv[2])
   print("starting unittests, '%s' and '%d'" % (hostname, port))
   unittest.main()

