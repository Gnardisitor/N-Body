clang: nbody.c
	clang -lm -lcurl -o nbody nbody.c

gcc: nbody.c
	gcc -lm -lcurl -o nbody nbody.c

clean:
	rm *.json nbody
