
#include "geometry.h"
#include "util.h"

#include "vtkCellType.h"

int geometry::ReadHexa(ifstream &ifs){
  if (IPRINT) cout << "geometry::ReadHexa"<<endl;

  // Read hexahedron number and number of subdivisions
  int number,nsub;
  ifs >> number >> nsub;
  if (IPRINT){
    cout << "number = " << number << endl;
    cout << "nsub   = " << nsub   << endl;
  }

  // Read points
  int nrpt = nrpt_hexa(nsub);
  vtkIdType *ip = new vtkIdType[nrpt];
  ReadPoints(ip,nrpt,ifs);

  // Add cells
  HexaCells(ip,nsub,number);

  delete [] ip;

  return nrpt;
}

void geometry::HexaCells(vtkIdType *ip,int nsub,int number){
  if (IPRINT) cout << "geometry::HexaCells"<<endl;

  int i,j,k;
  vtkIdType cell[8];

  //--------------------------------------------------------------------
  // Hexahedral cells for element volume
  //--------------------------------------------------------------------

  for (k=0; k<nsub; k++){
  for (j=0; j<nsub; j++){
  for (i=0; i<nsub; i++){
    cell[0] = ip[hexa(nsub,i,  j,  k  )];
    cell[1] = ip[hexa(nsub,i+1,j,  k  )];
    cell[2] = ip[hexa(nsub,i+1,j+1,k  )];
    cell[3] = ip[hexa(nsub,i,  j+1,k  )];
    cell[4] = ip[hexa(nsub,i,  j,  k+1)];
    cell[5] = ip[hexa(nsub,i+1,j,  k+1)];
    cell[6] = ip[hexa(nsub,i+1,j+1,k+1)];
    cell[7] = ip[hexa(nsub,i,  j+1,k+1)];
    vol->InsertNextCell(8,cell);
    volTyp->InsertNextType(VTK_HEXAHEDRON);
    volElt->InsertNextValue(number);
    volVis->InsertNextValue(1);
  }}}

  //--------------------------------------------------------------------
  // Quad cells for element faces
  //--------------------------------------------------------------------

  // side 1 (k=0)
  for (j=0; j<nsub; j++){
  for (i=0; i<nsub; i++){
    cell[0] = ip[hexa(nsub,i,  j,  0)];
    cell[1] = ip[hexa(nsub,i+1,j,  0)];
    cell[2] = ip[hexa(nsub,i+1,j+1,0)];
    cell[3] = ip[hexa(nsub,i,  j+1,0)];
    sur->InsertNextCell(4,cell);
    surTyp->InsertNextType(VTK_QUAD);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);
  }}

  // side 2 (k=nsub)
  for (j=0; j<nsub; j++){
  for (i=0; i<nsub; i++){
    cell[0] = ip[hexa(nsub,i,  j,  nsub)];
    cell[1] = ip[hexa(nsub,i+1,j,  nsub)];
    cell[2] = ip[hexa(nsub,i+1,j+1,nsub)];
    cell[3] = ip[hexa(nsub,i,  j+1,nsub)];
    sur->InsertNextCell(4,cell);
    surTyp->InsertNextType(VTK_QUAD);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);
  }}

  // side 3 (j=0)
  for (k=0; k<nsub; k++){
  for (i=0; i<nsub; i++){
    cell[0] = ip[hexa(nsub,i,  0,k  )];
    cell[1] = ip[hexa(nsub,i+1,0,k  )];
    cell[2] = ip[hexa(nsub,i+1,0,k+1)];
    cell[3] = ip[hexa(nsub,i,  0,k+1)];
    sur->InsertNextCell(4,cell);
    surTyp->InsertNextType(VTK_QUAD);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);
  }}

  // side 4 (i=nsub)
  for (k=0; k<nsub; k++){
  for (j=0; j<nsub; j++){
    cell[0] = ip[hexa(nsub,nsub,j,  k  )];
    cell[1] = ip[hexa(nsub,nsub,j+1,k  )];
    cell[2] = ip[hexa(nsub,nsub,j+1,k+1)];
    cell[3] = ip[hexa(nsub,nsub,j,  k+1)];
    sur->InsertNextCell(4,cell);
    surTyp->InsertNextType(VTK_QUAD);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);
  }}

  // side 5 (j=nsub)
  for (k=0; k<nsub; k++){
  for (i=0; i<nsub; i++){
    cell[0] = ip[hexa(nsub,i,  nsub,k  )];
    cell[1] = ip[hexa(nsub,i+1,nsub,k  )];
    cell[2] = ip[hexa(nsub,i+1,nsub,k+1)];
    cell[3] = ip[hexa(nsub,i,  nsub,k+1)];
    sur->InsertNextCell(4,cell);
    surTyp->InsertNextType(VTK_QUAD);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);
  }}

  // side 6 (i=0)
  for (k=0; k<nsub; k++){
  for (j=0; j<nsub; j++){
    cell[0] = ip[hexa(nsub,0,j,  k  )];
    cell[1] = ip[hexa(nsub,0,j+1,k  )];
    cell[2] = ip[hexa(nsub,0,j+1,k+1)];
    cell[3] = ip[hexa(nsub,0,j,  k+1)];
    sur->InsertNextCell(4,cell);
    surTyp->InsertNextType(VTK_QUAD);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);
  }}

  //--------------------------------------------------------------------
  // Linear cells for element edges
  //--------------------------------------------------------------------

  // edge 1 (j=0, k=0)
  for (i=0; i<nsub; i++){
    cell[0] = ip[hexa(nsub,i,0,0)];
    cell[1] = ip[hexa(nsub,i+1,0,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 2 (i=nsub, k=0)
  for (j=0; j<nsub; j++){
    cell[0] = ip[hexa(nsub,nsub,j,0)];
    cell[1] = ip[hexa(nsub,nsub,j+1,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 3 (j=nsub, k=0)
  for (i=0; i<nsub; i++){
    cell[0] = ip[hexa(nsub,i,nsub,0)];
    cell[1] = ip[hexa(nsub,i+1,nsub,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 4 (i=0, k=0)
  for (j=0; j<nsub; j++){
    cell[0] = ip[hexa(nsub,0,j,0)];
    cell[1] = ip[hexa(nsub,0,j+1,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 5 (j=0, k=nsub)
  for (i=0; i<nsub; i++){
    cell[0] = ip[hexa(nsub,i,0,nsub)];
    cell[1] = ip[hexa(nsub,i+1,0,nsub)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 6 (i=nsub, k=nsub)
  for (j=0; j<nsub; j++){
    cell[0] = ip[hexa(nsub,nsub,j,nsub)];
    cell[1] = ip[hexa(nsub,nsub,j+1,nsub)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 7 (j=nsub, k=nsub)
  for (i=0; i<nsub; i++){
    cell[0] = ip[hexa(nsub,i,nsub,nsub)];
    cell[1] = ip[hexa(nsub,i+1,nsub,nsub)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 8 (i=0, k=nsub)
  for (j=0; j<nsub; j++){
    cell[0] = ip[hexa(nsub,0,j,nsub)];
    cell[1] = ip[hexa(nsub,0,j+1,nsub)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 9 (i=0, j=0)
  for (k=0; k<nsub; k++){
    cell[0] = ip[hexa(nsub,0,0,k)];
    cell[1] = ip[hexa(nsub,0,0,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 10 (i=nsub, j=0)
  for (k=0; k<nsub; k++){
    cell[0] = ip[hexa(nsub,nsub,0,k)];
    cell[1] = ip[hexa(nsub,nsub,0,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 11 (i=nsub, j=nsub)
  for (k=0; k<nsub; k++){
    cell[0] = ip[hexa(nsub,nsub,nsub,k)];
    cell[1] = ip[hexa(nsub,nsub,nsub,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 12 (i=0, j=nsub)
  for (k=0; k<nsub; k++){
    cell[0] = ip[hexa(nsub,0,nsub,k)];
    cell[1] = ip[hexa(nsub,0,nsub,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }
}
