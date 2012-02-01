
include conf/mara.conf

MARA_E = bin/mara
MARA_A = lib/libmara.a

default : $(MARA_E)

$(MARA_A) : FORCE
	@make -C src

$(MARA_E) : $(MARA_A)
	$(CXX) -o $@ -Llib -lmara -llua -llunum $(CLIBS) $(HDF5_L) $(FFTW_L)

clean : 
	@make -C src clean
	@rm -f bin/mara

realclean : clean
	@rm -rf bin lib include

FORCE : 
