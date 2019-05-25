LIBS = 	-L/Users/johngao/work/thesis/CoinIpopt/build/lib		\
		-L/usr/local/lib										\
		-lcppad_ipopt -lipopt -llapack -lblas -ldl -lcoinhsl	\
		-L/usr/local/lib/gcc/5 									\
		-lgfortran -lm -lquadmath

INCS = 	-I/Users/johngao/work/thesis/CoinIpopt/build/include/coin/ThirdParty	\
		-I/Users/johngao/work/thesis/CoinIpopt/build/include 					\
		-I/Users/johngao/work/thesis/fast-cpp-csv-parser

censored_EM.o: censored_EM.cpp
	g++ -std=c++0x -c censored_EM.cpp ${INCS}

censored_EM: censored_EM.o
	g++ -std=c++0x -o censored_EM censored_EM.o ${LIBS}
