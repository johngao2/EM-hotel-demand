LIBS = -L/Users/johngao/work/thesis/cppad-20180000.0/build -L/Users/johngao/work/thesis/CoinIpopt/build/lib -lcppad_ipopt -lipopt -llapack -lblas -ldl -lcoinhsl -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L /lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath

INCS = -I/Users/johngao/work/thesis/CoinIpopt/build/include/coin/ThirdParty -I/Users/johngao/work/thesis/CoinIpopt/build/include -I/usr/local/include

get_started.o: get_started.cpp
	g++ -c get_started.cpp ${INCS}

get_started: get_started.o
	g++ -o get_started get_started.o ${LIBS}
