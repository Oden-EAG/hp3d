#
#  frontal solver routines (real version)
#
frontsol = $(SRC_DIR)/frontsol/bckwrdR.F \
	$(SRC_DIR)/frontsol/compltR.F \
	$(SRC_DIR)/frontsol/destR.F \
	$(SRC_DIR)/frontsol/desvecR.F \
	$(SRC_DIR)/frontsol/elmrhsR.F \
	$(SRC_DIR)/frontsol/elmsolR.F \
	$(SRC_DIR)/frontsol/frwcpR.F \
	$(SRC_DIR)/frontsol/frwrsR.F \
	$(SRC_DIR)/frontsol/loclsR.F \
	$(SRC_DIR)/frontsol/locluR.F \
	$(SRC_DIR)/frontsol/locrR.F \
	$(SRC_DIR)/frontsol/resolR.F \
	$(SRC_DIR)/frontsol/semrhsR.F \
	$(SRC_DIR)/frontsol/surfspR.F \
	$(SRC_DIR)/frontsol/surfssR.F \
	$(SRC_DIR)/frontsol/symasmR.F \
	$(SRC_DIR)/frontsol/symelmR.F \
	$(SRC_DIR)/frontsol/unsasmR.F \
	$(SRC_DIR)/frontsol/unselmR.F \
	$(SRC_DIR)/frontsol/zdirioR.F \
	$(SRC_DIR)/frontsol/zecondR.F \
	$(SRC_DIR)/frontsol/zerorR.F \
	$(SRC_DIR)/frontsol/zerosR.F \
	$(SRC_DIR)/frontsol/zerouR.F 
