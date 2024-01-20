CC=gcc
OFLAGS=-O3
CFLAGS=-march=native -g3
LFLAGS=-fopenmp -mfma -mavx2
SRC_DIR=nbody


all: nbody0 nbody1 nbody2 nbody3 nbody4 nbody5 nbody6 nbody7 nbody8 #nbody9

nbody0: $(SRC_DIR)/nbody.c
	$(CC) $(CFLAGS) $(OFLAGS) $(LFLAGS) $< -o $@ -lm

nbody1: $(SRC_DIR)/nbody1.c
	$(CC) $(CFLAGS) $(OFLAGS) $(LFLAGS) $< -o $@ -lm

nbody2: $(SRC_DIR)/nbody2.c
	$(CC) $(CFLAGS) $(OFLAGS) $(LFLAGS) $< -o $@ -lm

nbody3: $(SRC_DIR)/nbody3.c
	$(CC) $(CFLAGS) $(OFLAGS) $(LFLAGS) $< -o $@ -lm

nbody4: $(SRC_DIR)/nbody4.c
	$(CC) $(CFLAGS) $(OFLAGS) $(LFLAGS) $< -o $@ -lm

nbody5: $(SRC_DIR)/nbody5.c
	$(CC) $(CFLAGS) $(OFLAGS) $(LFLAGS) $< -o $@ -lm

nbody6: $(SRC_DIR)/nbody6.c
	$(CC) $(CFLAGS) $(OFLAGS) $(LFLAGS) $< -o $@ -lm

nbody7: $(SRC_DIR)/nbody7.c
	$(CC) $(CFLAGS) $(OFLAGS) $(LFLAGS) $< -o $@ -lm

nbody8: $(SRC_DIR)/nbody8.c
	$(CC) $(CFLAGS) $(OFLAGS) $(LFLAGS) $< -o $@ -lm

#nbody9: $(SRC_DIR)/nbody9.c
#	$(CC) $(CFLAGS) $(OFLAGS) $(LFLAGS) $< -o $@ -lm

clean:
	rm -Rf *~ nbody0 nbody1 nbody2 nbody3 nbody4 nbody5 nbody6 nbody7 nbody8 nbody9 *.optrpt

result_del:
	cd results/result_gcc/O1 && rm -Rf *~ *.dat
	cd results/result_gcc/O2 && rm -Rf *~ *.dat
	cd results/result_gcc/O3 && rm -Rf *~ *.dat
	cd results/result_gcc/OFast && rm -Rf *~ *.dat

	cd results/result_clang/O1 && rm -Rf *~ *.dat
	cd results/result_clang/O2 && rm -Rf *~ *.dat
	cd results/result_clang/O3 && rm -Rf *~ *.dat
	cd results/result_clang/OFast && rm -Rf *~ *.dat

plot_del:
	cd results/perf_plot && rm -Rf *~ *.png