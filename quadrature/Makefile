-include make.inc

default: mwrap mex

# Rebuild mylegendre.m and gateway.cpp via MWrap
mwrap: mylegendre.m gateway.cpp
mylegendre.m: mylegendre.mw
	$(MWRAP) -mex gateway -m mylegendre.m mylegendre.mw
gateway.cpp: mylegendre.mw
	$(MWRAP) -mex gateway -c gateway.cpp mylegendre.mw

# Rebuild gateway.mex* via MEX
mex: gateway.mex*
gateway.mex*: gateway.cpp
	$(MEX) gateway.cpp $(INCDIR) $(MEXFLAGS) $(LIBDIR) $(LIBS) -output gateway

# Remove the MEX interface
clean:
	rm -f gateway.mex*

# Remove the MEX interface, MATLAB caller, and MEX file
# Note: You will need MWrap to rebuild the deleted files!
mwrapclean: clean
	rm -f mylegendre.m gateway.cpp

.PHONY: mwrap mex clean mwrapclean
