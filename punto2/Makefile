graficas.py: estado_0.1.dat
		python graficas.py estado_0.1.dat

estado_0.1.dat: sod_test.x
		./sod_test.x 0.1

sod_test.x: tubo.c
		cc tubo.c -o sod_test.x -lm

clean:
	rm *.pdf *.dat *.x
