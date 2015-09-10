all: V1EXE

# von Tatjana
run: V1EXE
	./V1EXE -S 18 -B 3 -T 4 -d -7 -Y 4 -x 0.0 -M 0 -L 1 -R 500.0 

V1EXE: main.o topology.o mignicheweb.o getargs.o holling2.o evolveweb.o robustness.o gillespie.o
	g++ -o V1EXE main.o topology.o mignicheweb.o getargs.o holling2.o robustness.o gillespie.o evolveweb.o -lm -lgsl -lgslcblas -Wall

main.o: main.cpp
	g++ -c main.cpp
	
getargs.o: getargs.c
	g++ -c getargs.c

topology.o: topology.c
	g++ -c topology.c

mignicheweb.o: mignicheweb.c
	g++ -c mignicheweb.c

holling2.o: holling2.c
	g++ -c holling2.c

evolveweb.o: evolveweb.c
	g++ -c evolveweb.c

robustness.o: robustness.cpp
	g++ -c robustness.cpp

gillespie.o: gillespie.c
	g++ -c gillespie.c
	
# mit der Eingabe: 	make clean	können die .o Dateien und V1EXE gelöscht werden
clean:
	rm -rf *o V1EXE
