.PHONY: all clean

CFLAGS = -O3 -Wall -std=c++11
LDFLAGS+=-lm     # link to math library


TARGET=test_hmm

all: test_hmm.c model.h
	g++ $(LDFLAGS) test_hmm.c model.h -O3 -o test_hmm
# type make/make all to compile test_hmm

clean:
	$(RM) $(TARGET) train test_acc


train: main.cpp model.h
	@g++ main.cpp hmm.h model.h $(CFLAGS) -o train

test: test.cpp
	@g++ test.cpp $(CFLAGS) -o test_acc
	
