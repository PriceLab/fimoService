PORT=60000
FIMO="$(HOME)/meme/bin/fimo"
MOTIFS="../pfms/human-jaspar2018-hocomoco-swissregulon.meme"
TEAD1_MOTIF="../pfms/human-tead1-jaspar2018-ma0090.2.meme"
FOXP1_MOTIF="../pfms/human-foxp1-jaspar2018-ma0481.1.meme"
PYTHON_SCRIPT=$(HOME)/github/fimoService/server/runFimoServer.py
PYTHON=$(HOME)/anaconda3/bin/python
PYTHON_KHALEESI=$(HOME)/anaconda3/bin/python3
start:
	$(PYTHON) -i $(PYTHON_SCRIPT) $(PORT) $(FIMO) $(MOTIFS)&

tead1:
	$(PYTHON) -i $(PYTHON_SCRIPT) 60010 $(FIMO) $(TEAD1_MOTIF)&

tead1.khaleesi:
	nohup $(PYTHON_KHALEESI) -i $(PYTHON_SCRIPT) 60010 $(FIMO) $(TEAD1_MOTIF) &

foxp1.khaleesi:
	nohup $(PYTHON_KHALEESI) -i $(PYTHON_SCRIPT) 60011 $(FIMO) $(FOXP1_MOTIF) &



