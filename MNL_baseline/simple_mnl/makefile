LIBS = 	-L/Users/johngao/work/thesis/CoinIpopt/build/lib		\
		-L/usr/local/lib										\
		-lcppad_ipopt -lipopt -llapack -lblas -ldl -lcoinhsl	\
		-L/usr/local/lib/gcc/5 									\
		-lgfortran -lm -lquadmath

INCS = 	-I/Users/johngao/work/thesis/CoinIpopt/build/include/coin/ThirdParty	\
		-I/Users/johngao/work/thesis/CoinIpopt/build/include 					\
		-I/Users/johngao/work/thesis/fast-cpp-csv-parser

simple_mnl.o: simple_mnl.cpp
	g++ -std=c++0x -c simple_mnl.cpp ${INCS}

simple_mnl: simple_mnl.o
	g++ -std=c++0x -o simple_mnl simple_mnl.o ${LIBS}
