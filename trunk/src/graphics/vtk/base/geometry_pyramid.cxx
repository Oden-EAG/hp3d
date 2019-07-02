
#include "geometry.h"
#include "util.h"

#include "vtkCellType.h"

int geometry::ReadPyramid(ifstream &ifs){
  if (IPRINT) cout << "geometry::ReadPyramid"<<endl;

  // Read pyramid number and number of subdivisions
  int number,nsub;
  ifs >> number >> nsub;
  if (IPRINT){
    cout << "number = " << number << endl;
    cout << "nsub   = " << nsub   << endl;
  }

  // Read points
  int nrpt = nrpt_pyramid(nsub);
  vtkIdType *ip = new vtkIdType[nrpt];
  ReadPoints(ip,nrpt,ifs);

  // Add cells
  PyramidCells(ip,nsub,number);
  delete [] ip;

  return nrpt;
}

void geometry::PyramidCells(vtkIdType *ip,int nsub,int number){
  if (IPRINT) cout << "geometry::PyramidCells"<<endl;

  int i,j,k;
  vtkIdType cell[5];

  //--------------------------------------------------------------------
  // Pyramid and tetra cells for element volume
  //--------------------------------------------------------------------

  for (k=0; k<nsub; k++){
  for (j=0; j<nsub-k; j++){
  for (i=0; i<nsub-k; i++){

    // base pyramid
    cell[0] = ip[pyramid(nsub,i,  j,  k  )];
    cell[1] = ip[pyramid(nsub,i+1,j,  k  )];
    cell[2] = ip[pyramid(nsub,i+1,j+1,k  )];
    cell[3] = ip[pyramid(nsub,i,  j+1,k  )];
    cell[4] = ip[pyramid(nsub,i,  j,  k+1)];
    vol->InsertNextCell(5,cell);
    volTyp->InsertNextType(VTK_PYRAMID);
    volElt->InsertNextValue(number);
    volVis->InsertNextValue(1);

    // first tetra
    if (i<nsub-k-1){
      cell[0] = ip[pyramid(nsub,i+1,j,  k  )];
      cell[1] = ip[pyramid(nsub,i+1,j+1,k  )];
      cell[2] = ip[pyramid(nsub,i,  j,  k+1)];
      cell[3] = ip[pyramid(nsub,i+1,j,  k+1)];
      vol->InsertNextCell(4,cell);
      volTyp->InsertNextType(VTK_TETRA);
      volElt->InsertNextValue(number);
      volVis->InsertNextValue(1);
    }

    // second tetra
    if (j<nsub-k-1){
      cell[0] = ip[pyramid(nsub,i+1,j+1,k  )];
      cell[1] = ip[pyramid(nsub,i,  j+1,k  )];
      cell[2] = ip[pyramid(nsub,i,  j,  k+1)];
      cell[3] = ip[pyramid(nsub,i,  j+1,k+1)];
      vol->InsertNextCell(4,cell);
      volTyp->InsertNextType(VTK_TETRA);
      volElt->InsertNextValue(number);
      volVis->InsertNextValue(1);
    }

    // top pyramid
    if (i<nsub-k-1 && j<nsub-k-1){
      cell[0] = ip[pyramid(nsub,i+1,j+1,k+1)];
      cell[1] = ip[pyramid(nsub,i+1,j,  k+1)];
      cell[2] = ip[pyramid(nsub,i,  j,  k+1)];
      cell[3] = ip[pyramid(nsub,i,  j+1,k+1)];
      cell[4] = ip[pyramid(nsub,i+1,j+1,k  )];
      vol->InsertNextCell(5,cell);
      volTyp->InsertNextType(VTK_PYRAMID);
      volElt->InsertNextValue(number);
      volVis->InsertNextValue(1);
    }
  }}}

  //--------------------------------------------------------------------
  // Quad and triangle cells for element faces
  //--------------------------------------------------------------------

  // face 1 (k=0)
  for (j=0; j<nsub; j++){
  for (i=0; i<nsub; i++){
    cell[0] = ip[pyramid(nsub,i,  j,  0)];
    cell[1] = ip[pyramid(nsub,i,  j+1,0)];
    cell[2] = ip[pyramid(nsub,i+1,j+1,0)];
    cell[3] = ip[pyramid(nsub,i+1,j,  0)];
    sur->InsertNextCell(4,cell);
    surTyp->InsertNextType(VTK_QUAD);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);
  }}

  // face 2 (j=0)
  for (k=0; k<nsub; k++){
  for (i=0; i<nsub-k; i++){

    // lower triangle
    cell[0] = ip[pyramid(nsub,i,  0,k  )];
    cell[1] = ip[pyramid(nsub,i+1,0,k  )];
    cell[2] = ip[pyramid(nsub,i  ,0,k+1)];
    sur->InsertNextCell(3,cell);
    surTyp->InsertNextType(VTK_TRIANGLE);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);

    // upper triangle
    if (i<nsub-k-1){
      cell[0] = ip[pyramid(nsub,i,  0,k+1)];
      cell[1] = ip[pyramid(nsub,i+1,0,k  )];
      cell[2] = ip[pyramid(nsub,i+1,0,k+1)];
      sur->InsertNextCell(3,cell);
      surTyp->InsertNextType(VTK_TRIANGLE);
      surElt->InsertNextValue(number);
      surVis->InsertNextValue(1);
    }
  }}

  // face 3 (i=nsub-k)
  for (k=0; k<nsub; k++){
  for (j=0; j<nsub-k; j++){

    // lower triangle
    cell[0] = ip[pyramid(nsub,nsub-k,  j,  k  )];
    cell[1] = ip[pyramid(nsub,nsub-k,  j+1,k  )];
    cell[2] = ip[pyramid(nsub,nsub-k-1,j,  k+1)];
    sur->InsertNextCell(3,cell);
    surTyp->InsertNextType(VTK_TRIANGLE);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);

    // upper triangle
    if (j<nsub-k-1){
      cell[0] = ip[pyramid(nsub,nsub-k,  j+1,k  )];
      cell[1] = ip[pyramid(nsub,nsub-k-1,j+1,k+1)];
      cell[2] = ip[pyramid(nsub,nsub-k-1,j,  k+1)];
      sur->InsertNextCell(3,cell);
      surTyp->InsertNextType(VTK_TRIANGLE);
      surElt->InsertNextValue(number);
      surVis->InsertNextValue(1);
    }
  }}

  // face 4 (j=nsub-k)
  for (k=0; k<nsub; k++){
  for (i=0; i<nsub-k; i++){

    // lower triangle
    cell[0] = ip[pyramid(nsub,i,  nsub-k,  k  )];
    cell[1] = ip[pyramid(nsub,i,  nsub-k-1,k+1)];
    cell[2] = ip[pyramid(nsub,i+1,nsub-k,  k  )];
    sur->InsertNextCell(3,cell);
    surTyp->InsertNextType(VTK_TRIANGLE);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);

    // upper triangle
    if (i<nsub-k-1){
      cell[0] = ip[pyramid(nsub,i,  nsub-k-1,k+1)];
      cell[1] = ip[pyramid(nsub,i+1,nsub-k-1,k+1)];
      cell[2] = ip[pyramid(nsub,i+1,nsub-k,  k  )];
      sur->InsertNextCell(3,cell);
      surTyp->InsertNextType(VTK_TRIANGLE);
      surElt->InsertNextValue(number);
      surVis->InsertNextValue(1);
    }
  }}

  // face 5 (i=0)
  for (k=0; k<nsub; k++){
  for (j=0; j<nsub-k; j++){

    // lower triangle
    cell[0] = ip[pyramid(nsub,0,j,  k  )];
    cell[1] = ip[pyramid(nsub,0,j,  k+1)];
    cell[2] = ip[pyramid(nsub,0,j+1,k  )];
    sur->InsertNextCell(3,cell);
    surTyp->InsertNextType(VTK_TRIANGLE);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);

    // upper triangle
    if (j<nsub-k-1){
      cell[0] = ip[pyramid(nsub,0,j,  k+1)];
      cell[1] = ip[pyramid(nsub,0,j+1,k+1)];
      cell[2] = ip[pyramid(nsub,0,j+1,k  )];
      sur->InsertNextCell(3,cell);
      surTyp->InsertNextType(VTK_TRIANGLE);
      surElt->InsertNextValue(number);
      surVis->InsertNextValue(1);
    }
  }}

  //--------------------------------------------------------------------
  // Linear cells for element edges
  //--------------------------------------------------------------------

  // edge 1 (k=0, j=0)
  for (i=0; i<nsub; i++){
    cell[0] = ip[pyramid(nsub,i,0,0)];
    cell[1] = ip[pyramid(nsub,i+1,0,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 2 (k=0, i=nsub)
  for (j=0; j<nsub; j++){
    cell[0] = ip[pyramid(nsub,nsub,j,0)];
    cell[1] = ip[pyramid(nsub,nsub,j+1,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 3 (k=0, j=nsub)
  for (i=0; i<nsub; i++){
    cell[0] = ip[pyramid(nsub,i,nsub,0)];
    cell[1] = ip[pyramid(nsub,i+1,nsub,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 4 (k=0, i=0)
  for (j=0; j<nsub; j++){
    cell[0] = ip[pyramid(nsub,0,j,0)];
    cell[1] = ip[pyramid(nsub,0,j+1,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 5 (i=0, j=0)
  for (k=0; k<nsub; k++){
    cell[0] = ip[pyramid(nsub,0,0,k)];
    cell[1] = ip[pyramid(nsub,0,0,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 6 (i=nsub-k, j=0)
  for (k=0; k<nsub; k++){
    cell[0] = ip[pyramid(nsub,nsub-k,0,k)];
    cell[1] = ip[pyramid(nsub,nsub-k-1,0,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 7 (i=nsub-k, j=nsub-k)
  for (k=0; k<nsub; k++){
    cell[0] = ip[pyramid(nsub,nsub-k,nsub-k,k)];
    cell[1] = ip[pyramid(nsub,nsub-k-1,nsub-k-1,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 8 (i=0, j=nsub-k)
  for (k=0; k<nsub; k++){
    cell[0] = ip[pyramid(nsub,0,nsub-k,k)];
    cell[1] = ip[pyramid(nsub,0,nsub-k-1,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }
}
