SRC0 += $(SRC_DIR)/constrs/module/constrained_nodes.F
#
ifeq ($(REF),1)
constrs =   $(SRC_DIR)/constrs/celem.F \
	$(SRC_DIR)/constrs/get_connect_info.F \
	$(SRC_DIR)/constrs/logic_new.F \
	$(SRC_DIR)/constrs/logic_nodes.F \
	$(SRC_DIR)/constrs/nodcor.F \
	$(SRC_DIR)/constrs/nodcor_vert.F \
	$(SRC_DIR)/constrs/solelm.F 
  
else
constrs =   $(SRC_DIR)/constrs/celem.F \
	$(SRC_DIR)/constrs/get_connect_info.F \
	$(SRC_DIR)/constrs/logic.F \
	$(SRC_DIR)/constrs/logic_nodes.F \
	$(SRC_DIR)/constrs/nodcor.F \
	$(SRC_DIR)/constrs/nodcor_vert.F \
	$(SRC_DIR)/constrs/solelm.F 
endif



