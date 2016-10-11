PROJ_DIR=/home/jerry/usr/workspace/riemann

EXE  = $(PROJ_DIR)/build/bin/test_vol_frame
MESH = $(PROJ_DIR)/dat/tets/sculpture.c59k.vtk
WS   = 1e0
WA   = 1e3
EPS  = 1e-8
MAXITS = 3000

MESH_PREFIX = $(notdir $(MESH))
IDENT = $(basename $(MESH_PREFIX))
OUTDIR = $(PROJ_DIR)/result/volume_frame/$(IDENT)-ws$(WS)-wa$(WA)/


$(OUTDIR)/zyz.txt: $(EXE) $(MESH)
	mkdir -p $(OUTDIR)
	$(EXE) mesh=$(MESH) out_dir=$(OUTDIR) weight.smooth=$(WS) weight.align=$(WA) lbfgs.epsf=$(EPS) lbfgs.maxits=$(MAXITS) | tee $(OUTDIR)/log.txt


vis_sing:
	$(FF_EXE) prog=draw_3d_frame_sing tet=$(OUTDIR)/tet.vtk zyz=$(OUTDIR)/zyz.txt out=$(OUTDIR)/sing.vtk


clean:
	-rm -rf $(OUTDIR)/*
