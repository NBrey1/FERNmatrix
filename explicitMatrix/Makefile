CC=g++
CCDBGFLAGS=-g
CCOFLAGS=-O2
CFLAGS=$(CCDBGFLAGS) $(CCOFLAGS) -I.

#compiler flags
LIBS=-lgsl -lgslcblas -lgmp
all: explicitMatrix.cpp #EMATS.cpp # TEST.cpp #NATS.cpp CEM.cpp #dtrange.cpp CATS.cpp
#	$(CC) -o explicitMatrix explicitMatrix.cpp $(CFLAGS) $(LIBS)
#	$(CC) -o EM.out explicitMatrix.cpp $(CFLAGS) $(LIBS)
	$(CC) -o postFERN postFERN.cpp $(CFLAGS) $(LIBS)
#	$(CC) -o analysis.out postFERN.cpp $(CFLAGS) $(LIBS)
#	$(CC) -o test testread.cpp $(CFLAGS) $(LIBS)
#	$(CC) -o test.out testread.cpp $(CFLAGS) $(LIBS)
#	$(CC) -o FERNA FERNA.cpp $(CFLAGS) $(LIBS)
	$(CC) -o FERN FERN.cpp $(CFLAGS) $(LIBS)

vectorOps: vectorOps.cpp
	$(CC) -c vectorOps.cpp

# make sure boost-test is installed before calling this target
# dnf install boost-test boost-devel
BOOST_TEST_LIBS=-lboost_unit_test_framework
vectorTest: vectorTest.cpp vectorOps
	$(CC) -o vectorTest vectorTest.cpp vectorOps.o $(CFLAGS) $(LIBS) $(BOOST_TEST_LIBS)

.PHONY: clean

clean:
	rm -f test *.0
	rm -f explicitMatrix vectorTest *.o
	rm -f postFERN *.o
	rm -f FERNA *.o
