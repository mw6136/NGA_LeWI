include ../Makefile.in

BINFILES = init_flow data2plot3d config2plot3d chucks2arts \
	   config2ensight data2ensight createChemtable \
	   dns2les interpolatedata \
	   biinflow2inflow chemtable2plot3d \
	   editData editInflow statdata2plot3d statconfig2plot3d spectrum \
	   rescaler compare mergeData bl2ml channel2jet \
	   create1DChemtable createChemtable mergeInflow \
	   createMMMChemtable createPremtable createUnsteadyPremtable \
	   combineInflow createUnsteadyChemtable viewTable \
	   InflowStat extractData \
	   planeStats volumeStats readPlaneData comparePlaneData \
	   dmd_p1 dmd dmd_p3 post_dmd pre_dmd \
	   volume2ensight volume_vort dmd_vort \
	   volume_diff dmd_diff \
	   volume_react_rate dmd_react_rate \
	   FlameRegimes mmsTool rescaleChannel \
	   FlameRegimes mmsTool rescaleChannel dnsbox_corr \
	   average_data stat2ascii dataslice2csv editVolumeData

all:	Makefile
	@for i in $(BINFILES); do make -C $$i; done

clean: 
	@for i in $(BINFILES); do make clean -C $$i; done

# premtable2plot3d
