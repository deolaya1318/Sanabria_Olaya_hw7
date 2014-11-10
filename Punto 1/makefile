plot.py: string_rho.dat
	python plot.py

string_rho.dat: string.x
	./string.x 0.01

string.x: string.c
	cc string.c -o string.x -lm

clean:
	rm *.x *.dat *.pdf