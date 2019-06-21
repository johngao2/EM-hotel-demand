LIBS = 	-L/Users/johngao/work/thesis/CoinIpopt/build/lib		\
		-L/usr/local/lib										\
		-lcppad_ipopt -lipopt -llapack -lblas -ldl -lcoinhsl	\
		-L/usr/local/lib/gcc/5 									\
		-lgfortran -lm -lquadmath

INCS = 	-I/Users/johngao/work/thesis/CoinIpopt/build/include/coin/ThirdParty	\
		-I/Users/johngao/work/thesis/CoinIpopt/build/include 					\
		-I/Users/johngao/work/thesis/fast-cpp-csv-parser						\
		-I/Users/johngao/work/thesis/boost_1_62_0								\

EM_cabot_s4.o: EM_cabot_s4.cpp
	g++ -std=c++0x -c EM_cabot_s4.cpp ${INCS}

EM_cabot_s4: EM_cabot_s4.o
	g++ -std=c++0x -o EM_cabot_s4 EM_cabot_s4.o ${LIBS}
