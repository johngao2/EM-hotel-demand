LIBS = 	-L/Users/johngao/work/thesis/CoinIpopt/build/lib	\
		-L/usr/local/lib									\
		-lcppad_ipopt -lipopt -llapack -lblas -ldl -lcoinhsl	\
		-L/usr/local/lib/gcc/5 \
		-lgfortran -lm -lquadmath

INCS = 	-I/home/mnediak/CoinIpopt/build/include/coin/ThirdParty	\
		-I/home/mnediak/CoinIpopt/build/include 				\
		-I/home/mnediak/include

get_started.o: get_started.cpp
	g++ -c get_started.cpp ${INCS}

get_started: get_started.o
	g++ -o get_started get_started.o ${LIBS}
