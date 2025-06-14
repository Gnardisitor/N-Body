clang: nbody.c
	clang -lm -o nbody nbody.c

gcc: nbody.c
	gcc -lm -o nbody nbody.c

test:
	echo `time ./nbody 8 1000 0.1 euler`

clean:
	rm *.csv
