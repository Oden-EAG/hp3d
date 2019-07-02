module node_types
  implicit none
  integer, parameter :: NODE_TYPE_EMPTY = -1, &
                        NODE_TYPE_VERTEX = 0, &
                        NODE_TYPE_EDGE = 1, &
                        NODE_TYPE_TRI = 2, &
                        NODE_TYPE_QUAD = 3, &
                        NODE_TYPE_TETRA = 4, &
                        NODE_TYPE_PRISM = 5, &
                        NODE_TYPE_PYRAMID = 6, &
                        NODE_TYPE_HEXA = 7
end module node_types
