#ifdef _INT_NODE_TYPE
#define NODE_TYPE_TYPE integer
#define NODE_TYPE_EMPTY -1
#define NODE_TYPE_VERTEX 0
#define NODE_TYPE_EDGE 1
#define NODE_TYPE_TRI 2
#define NODE_TYPE_QUAD 3
#define NODE_TYPE_TETRA 4
#define NODE_TYPE_PRISM 5
#define NODE_TYPE_PYRAM 6
#define NODE_TYPE_HEXA 7
#define ELEM_TYPE_TETRA 4
#define ELEM_TYPE_PRISM 5
#define ELEM_TYPE_PYRAM 6
#define ELEM_TYPE_HEXA 7

#else
#define NODE_TYPE_TYPE character(len=4)
#define NODE_TYPE_EMPTY "none"
#define NODE_TYPE_VERTEX "vert"
#define NODE_TYPE_EDGE "medg"
#define NODE_TYPE_TRI "mdlt"
#define NODE_TYPE_QUAD "mdlq"
#define NODE_TYPE_TETRA "mdln"
#define NODE_TYPE_PRISM "mdlp"
#define NODE_TYPE_PYRAM "mdld"
#define NODE_TYPE_HEXA "mdlb"
#define ELEM_TYPE_TETRA "tetr"
#define ELEM_TYPE_PRISM "pris"
#define ELEM_TYPE_PYRAM "pyra"
#define ELEM_TYPE_HEXA "bric"
#endif
