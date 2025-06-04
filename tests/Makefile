F90=gfortran
FLAGS=-I../src -DDOUBLE -ffree-line-length-none -O2

BUILD=build
SRC=../src/kdtree2.F90

.PHONY: all test clean

all: test

$(BUILD)/test_pq: test_pq.f90 $(SRC)
	mkdir -p $(BUILD)
	$(F90) $(FLAGS) -c $(SRC) -o $(BUILD)/kdtree2.o
	$(F90) $(FLAGS) test_pq.f90 $(BUILD)/kdtree2.o -o $(BUILD)/test_pq

test: $(BUILD)/test_pq
	./$(BUILD)/test_pq

clean:
	rm -rf $(BUILD)
