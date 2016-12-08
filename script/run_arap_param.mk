EXE  =../build/bin/test_arap_param
MODEL=../dat/lilium.obj

arap_param/lscm_param.obj: $(EXE) $(MODEL)
	@$(EXE) $(MODEL)

clean:
	@rm -rf arap_param/*
