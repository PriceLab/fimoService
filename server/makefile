PORT=60000
FIMO="$(HOME)/meme/bin/fimo"
MOTIFS="../pfms/human-jaspar2018-hocomoco-swissregulon.meme"
PYTHON_SCRIPT=$(HOME)/github/fimoService/server/runFimoServer.py

start:
	python -i $(PYTHON_SCRIPT) $(PORT) $(FIMO) $(MOTIFS)&

