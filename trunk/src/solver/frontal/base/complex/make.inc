#
#  frontal solver routines (complex version)
#
frontsolz = $(SRC_DIR)/frontsolz/bckwrd.F \
	$(SRC_DIR)/frontsolz/complt.F \
	$(SRC_DIR)/frontsolz/dest.F \
	$(SRC_DIR)/frontsolz/desvec.F \
	$(SRC_DIR)/frontsolz/elmrhs.F \
	$(SRC_DIR)/frontsolz/elmsol.F \
	$(SRC_DIR)/frontsolz/frwcp.F \
	$(SRC_DIR)/frontsolz/frwrs.F \
	$(SRC_DIR)/frontsolz/locls.F \
	$(SRC_DIR)/frontsolz/loclu.F \
	$(SRC_DIR)/frontsolz/locr.F \
	$(SRC_DIR)/frontsolz/resol.F \
	$(SRC_DIR)/frontsolz/semrhs.F \
	$(SRC_DIR)/frontsolz/surfsp.F \
	$(SRC_DIR)/frontsolz/surfss.F \
	$(SRC_DIR)/frontsolz/symasm.F \
	$(SRC_DIR)/frontsolz/symelm.F \
	$(SRC_DIR)/frontsolz/unsasm.F \
	$(SRC_DIR)/frontsolz/unselm.F \
	$(SRC_DIR)/frontsolz/zdirio.F \
	$(SRC_DIR)/frontsolz/zdiriodum.F \
	$(SRC_DIR)/frontsolz/zecond.F \
	$(SRC_DIR)/frontsolz/zeror.F \
	$(SRC_DIR)/frontsolz/zeros.F \
	$(SRC_DIR)/frontsolz/zerou.F 
