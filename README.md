# fimoService
The meme suite's fimo (Find Individual Motif Occurences) as a microservice.


William Noble's lab at the University of Washington publishes the [MEME
Suite](http://meme-suite.org/index.html) - Motif-based sequence analysis tools.

We here provide R and python convenience functions for running
[FIMO](http://meme-suite.org/doc/fimo.html?man_type=web) on specified DNA sequence/s.

To demonstrate the FIMO service, we use an example from the 2013 paper,
["Highly recurrent TERT promoter mutations in human melanoma"](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4423787),
which reports the creation of a transcript factor binding motif due to a single-base mutation.

This microservice has a simple design and implementation:

  * The MEME suite and FIMO are installed on your localhost, or (at the ISB) on whovian.
  * A motif library file is avaialable at a known location on that host.
  * The FimoServer, written Python, opens a [ZMQ](http://zeromq.org/) REPLY socket on
    a specified port, and listens for ZMQ REQUEST messages
  * These messages contain a JSON named list of DNA sequences
  * The python FimoServer executes FIMO as a separate process and collects the results.
  * The results are captured into a dataframe (in R, or in python pandas) and is returned (as shown below)


Start the FimoServer.  In many cases, those wishing to use the FimoService will not need
to do this because (as at the ISB, behind the firewall, on whovian) the service is
already running.
```
PORT=5558
FIMO="/Users/paul/meme/bin/fimo"
MOTIF_FILE="/Users/paul/s/work/priceLab/cory/footprintWorkflow/JASPAR_CORE_plus_seth.meme"
python -i runServer.py $PORT $FIMO $MOTIF_FILE
```

Test this out with the python FimoClient module:

```
python testFimoClient.py localhost 5558
```

At the ISB, with no need to start your own server (and thus no need to install the MEME suite):
```
python testFimoClient.py whovian 5558
```

In R, after installing two prerequisite packages:

````
source("http://bioconductor.org/biocLite.R")
biocLite("pbdZMQ")
biocLite("devtools")
library(devtools)
install_github("PriceLab/fimoService", subdir="client-R/FimoClient")
``` 

With this installation complete, try it out:

```
library(FimoClient)
fimo <- FimoClient("whovian", 5558, quiet=TRUE)
sequences <- list(tert_wt1="CCCGGAGGGGG", tert_wt2="CCCGGGAGGGG", tert_mut="CCCCTTCCGGG")
tbl <- requestMatch(fimo, sequences)
```
Producing this output.  Notice that only the tert_mut sequence matches any PWM, as suggested
by the Tert Promoter Mutations paper:

|X.pattern.name | sequence.name | start | stop | strand | score | p.value | q.value | matched.sequence|
|---------------|:--------------|:------|:-----|:-------|:------|:--------|:--------|-----------------|
|             MA0076.2  |    tert_mut  |   1  | 11   |   + | 13.1818 | 2.22e-05 | 0.000133 |    CCCCTTCCGGG
| ELK1,4_GABP{A,B1}.p3  |    tert_mut  |   2  | 11   |   - | 11.7000 | 4.50e-05 | 0.000539 |     CCCGGAAGGG
|          ETV6_full_2  |    tert_mut  |   2  | 11   |   - | 10.7245 | 9.30e-05 | 0.001120 |      CCCGGAAGGG
|             MA0645.1  |    tert_mut  |   2  | 11   |   - | 10.7308 | 9.54e-05 | 0.001150 |     CCCGGAAGGG



