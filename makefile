# Compile with root the two files in the make directory
BUILDDIR = make
ROOTCMD = root -l -b -q

all: $(BUILDDIR)/mk_DijetHistosFill.C $(BUILDDIR)/mk_CondFormats.C
	$(ROOTCMD) $(BUILDDIR)/mk_CondFormats.C
	$(ROOTCMD) $(BUILDDIR)/mk_DijetHistosFill.C
	
clean:
	rm -f *.so *.d *.pcm
	rm -f $(BUILDDIR)/*.so $(BUILDDIR)/*.d $(BUILDDIR)/*.pcm
	rm -f src/*.so src/*.d src/*.pcm
	rm -f CondFormats/JetMETObjects/src/*.so CondFormats/JetMETObjects/src/*.d CondFormats/JetMETObjects/src/*.pcm
