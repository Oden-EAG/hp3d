
#include "geometry.h"
#include "util.h"

#include "vtkCellType.h"

int geometry::ReadPrism(ifstream &ifs){
  if (IPRINT) cout << "geometry::ReadPrism"<<endl;

  // Read prism number and number of subdivisions
  int number,nsub;
  ifs >> number >> nsub;
  if (IPRINT){
    cout << "number = " << number << endl;
    cout << "nsub   = " << nsub   << endl;
  }

  // Read points
  int nrpt = nrpt_prism(nsub);
  vtkIdType *ip = new vtkIdType[nrpt];
  ReadPoints(ip,nrpt,ifs);

  // Add cells
  PrismCells(ip,nsub,number);

  delete [] ip;

  return nrpt;
}

void geometry::PrismCells(vtkIdType *ip,int nsub,int number){
  if (IPRINT) cout << "geometry::PrismCells"<<endl;

  int i,j,k;
  vtkIdType cell[6];

  //--------------------------------------------------------------------
  // Wedge cells for element volume
  //--------------------------------------------------------------------

  for (k=0; k<nsub; k++){
  for (j=0; j<nsub; j++){
  for (i=0; i<nsub-j; i++){

    cell[0] = ip[prism(nsub,i,  j,  k)];
    cell[1] = ip[prism(nsub,i,  j+1,k)];
    cell[2] = ip[prism(nsub,i+1,j,  k)];
    cell[3] = ip[prism(nsub,i,  j,  k+1)];
    cell[4] = ip[prism(nsub,i,  j+1,k+1)];
    cell[5] = ip[prism(nsub,i+1,j,  k+1)];
    vol->InsertNextCell(6,cell);
    volTyp->InsertNextType(VTK_WEDGE);
    volElt->InsertNextValue(number);
    volVis->InsertNextValue(1);

    if (i<nsub-j-1){
      cell[0] = ip[prism(nsub,i+1,j+1,k)];
      cell[1] = ip[prism(nsub,i+1,j,  k)];
      cell[2] = ip[prism(nsub,i,  j+1,k)];
      cell[3] = ip[prism(nsub,i+1,j+1,k+1)];
      cell[4] = ip[prism(nsub,i+1,j,  k+1)];
      cell[5] = ip[prism(nsub,i,  j+1,k+1)];
      vol->InsertNextCell(6,cell);
      volTyp->InsertNextType(VTK_WEDGE);
      volElt->InsertNextValue(number);
      volVis->InsertNextValue(1);
    }
  }}}

  //--------------------------------------------------------------------
  // Triangle and quad cells for element faces
  //--------------------------------------------------------------------

  // face 1 (k=0)
  for (j=0; j<nsub; j++){
  for (i=0; i<nsub-j; i++){

    cell[0] = ip[prism(nsub,i,  j,  0)];
    cell[1] = ip[prism(nsub,i+1,j,  0)];
    cell[2] = ip[prism(nsub,i,  j+1,0)];
    sur->InsertNextCell(3,cell);
    surTyp->InsertNextType(VTK_TRIANGLE);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);

    if (i<nsub-j-1){
      cell[0] = ip[prism(nsub,i+1,j+1,0)];
      cell[1] = ip[prism(nsub,i,  j+1,0)];
      cell[2] = ip[prism(nsub,i+1,j,  0)];
      sur->InsertNextCell(3,cell);
      surTyp->InsertNextType(VTK_TRIANGLE);
      surElt->InsertNextValue(number);
      surVis->InsertNextValue(1);
    }
  }}

  // face 2 (k=nsub)
  for (j=0; j<nsub; j++){
  for (i=0; i<nsub-j; i++){

    cell[0] = ip[prism(nsub,i,  j,  nsub)];
    cell[1] = ip[prism(nsub,i+1,j,  nsub)];
    cell[2] = ip[prism(nsub,i,  j+1,nsub)];
    sur->InsertNextCell(3,cell);
    surTyp->InsertNextType(VTK_TRIANGLE);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);

    if (i<nsub-j-1){
      cell[0] = ip[prism(nsub,i+1,j+1,nsub)];
      cell[1] = ip[prism(nsub,i,  j+1,nsub)];
      cell[2] = ip[prism(nsub,i+1,j,  nsub)];
      sur->InsertNextCell(3,cell);
      surTyp->InsertNextType(VTK_TRIANGLE);
      surElt->InsertNextValue(number);
      surVis->InsertNextValue(1);
    }
  }}

  // face 3 (j=0)
  for (k=0; k<nsub; k++){
  for (i=0; i<nsub; i++){
    cell[0] = ip[prism(nsub,i,  0,k  )];
    cell[1] = ip[prism(nsub,i+1,0,k  )];
    cell[2] = ip[prism(nsub,i+1,0,k+1)];
    cell[3] = ip[prism(nsub,i,  0,k+1)];
    sur->InsertNextCell(4,cell);
    surTyp->InsertNextType(VTK_QUAD);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);
  }}

  // face 4 (i=nsub-j)
  for (k=0; k<nsub; k++){
  for (j=0; j<nsub; j++){
    cell[0] = ip[prism(nsub,nsub-j,  j,  k  )];
    cell[1] = ip[prism(nsub,nsub-j-1,j+1,k  )];
    cell[2] = ip[prism(nsub,nsub-j-1,j+1,k+1)];
    cell[3] = ip[prism(nsub,nsub-j,  j,  k+1)];
    sur->InsertNextCell(4,cell);
    surTyp->InsertNextType(VTK_QUAD);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);
  }}

  // face 5 (i=0)
  for (k=0; k<nsub; k++){
  for (j=0; j<nsub; j++){
    cell[0] = ip[prism(nsub,0,j+1,k  )];
    cell[1] = ip[prism(nsub,0,j,  k  )];
    cell[2] = ip[prism(nsub,0,j,  k+1)];
    cell[3] = ip[prism(nsub,0,j+1,k+1)];
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
    cell[0] = ip[prism(nsub,i,0,0)];
    cell[1] = ip[prism(nsub,i+1,0,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 2 (i=nsub-j, k=0)
  for (j=0; j<nsub; j++){
    cell[0] = ip[prism(nsub,nsub-j,j,0)];
    cell[1] = ip[prism(nsub,nsub-j-1,j+1,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 3 (i=0, k=0)
  for (j=0; j<nsub; j++){
    cell[0] = ip[prism(nsub,0,j,0)];
    cell[1] = ip[prism(nsub,0,j+1,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 4 (j=0, k=nsub)
  for (i=0; i<nsub; i++){
    cell[0] = ip[prism(nsub,i,0,nsub)];
    cell[1] = ip[prism(nsub,i+1,0,nsub)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 5 (i=nsub-j, k=nsub)
  for (j=0; j<nsub; j++){
    cell[0] = ip[prism(nsub,nsub-j,j,nsub)];
    cell[1] = ip[prism(nsub,nsub-j-1,j+1,nsub)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 6 (i=0, k=nsub)
  for (j=0; j<nsub; j++){
    cell[0] = ip[prism(nsub,0,j,nsub)];
    cell[1] = ip[prism(nsub,0,j+1,nsub)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 7 (i=0, j=0)
  for (k=0; k<nsub; k++){
    cell[0] = ip[prism(nsub,0,0,k)];
    cell[1] = ip[prism(nsub,0,0,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 8 (i=nsub, j=0)
  for (k=0; k<nsub; k++){
    cell[0] = ip[prism(nsub,nsub,0,k)];
    cell[1] = ip[prism(nsub,nsub,0,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 9 (i=0, j=nsub)
  for (k=0; k<nsub; k++){
    cell[0] = ip[prism(nsub,0,nsub,k)];
    cell[1] = ip[prism(nsub,0,nsub,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }
}


