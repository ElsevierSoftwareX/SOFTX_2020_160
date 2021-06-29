ROOTINC="${ROOTSYS}/include"
ROOTLIB="${ROOTSYS}/lib"

ifdef _USE_ICC
CPP         = icc
else
CPP         = g++
endif


ifdef _USE_CPP11
CPPFLAGS += -std=c++11
endif


SOURCES=Function.cc StdFunction.cc NuFunction.cc CompositeFunction.cc CosFourier.cc SaveFilters.cc

main :  
	$(CPP) -I${ROOTINC}  ${SOURCES} -L${ROOTLIB} $(CPPFLAGS) -lMathCore -lstdc++ -lm -o ComputeWDMFourier

clean:
	rm -f ComputeWDMFourier
