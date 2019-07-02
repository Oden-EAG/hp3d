# add module
SRC0 += $(SRC_DIR)/constrs_util/module/constraints.F
#
constrs_util = $(SRC_DIR)/constrs_util/setcnstr_edge.F \
	       $(SRC_DIR)/constrs_util/setcnstr_trian.F \
	       $(SRC_DIR)/constrs_util/setcnstr_trian_aniso.F \
	       $(SRC_DIR)/constrs_util/setcnstr_trian_iso.F \
	       $(SRC_DIR)/constrs_util/test_trian_aniso.F 
  



