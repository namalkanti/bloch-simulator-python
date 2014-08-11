bloch: bloch/bloch_simulator.c bloch/bloch.c
	python setup.up build_ext
	mv build/lib.linux-x86_64-3.4/bloch/bloch_simulator.cpython-34m.so bloch/
