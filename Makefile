CC = g++

parallel_bmesolver: 
	g++ -std=c++11 -D_REENTRANT -w -O3 -Xpreprocessor -fopenmp -lomp -I${XPRESSDIR}/include/ -L${XPRESSDIR}/lib Main.cpp -o BME-Solver -lxprs -lxprb
	cp BME-Solver ../Tests/
	
clean:	
	rm -f BME-Solver;
	rm -f ../Tests/BME-Solver
	clear;
