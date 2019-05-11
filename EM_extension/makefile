LIBS = 	-L/Users/johngao/work/thesis/CoinIpopt/build/lib		\
		-L/usr/local/lib										\
		-lcppad_ipopt -lipopt -llapack -lblas -ldl -lcoinhsl	\
		-L/usr/local/lib/gcc/5 									\
		-lgfortran -lm -lquadmath

INCS = 	-I/Users/johngao/work/thesis/CoinIpopt/build/include/coin/ThirdParty	\
		-I/Users/johngao/work/thesis/CoinIpopt/build/include 					\
		-I/Users/johngao/work/thesis/fast-cpp-csv-parser

EM_cabot.o: EM_cabot.cpp
	g++ -std=c++0x -c EM_cabot.cpp ${INCS}

EM_cabot: EM_cabot.o
	g++ -std=c++0x -o EM_cabot EM_cabot.o ${LIBS}
