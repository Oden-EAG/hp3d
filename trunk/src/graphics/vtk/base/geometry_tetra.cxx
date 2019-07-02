
#include "geometry.h"
#include "util.h"

#include "vtkCellType.h"

int geometry::ReadTetra(ifstream &ifs){
  if (IPRINT) cout << "geometry::ReadTetra"<<endl;

  // Read tetrahedron number and number of subdivisions
  int number,nsub;
  ifs >> number >> nsub;
  if (IPRINT){
    cout << "number = " << number << endl;
    cout << "nsub   = " << nsub   << endl;
  }

  // Read points
  int nrpt = nrpt_tetra(nsub);
  vtkIdType *ip = new vtkIdType[nrpt];
  ReadPoints(ip,nrpt,ifs);

  // Add cells
  TetraCells(ip,nsub,number);

  delete [] ip;

  return nrpt;
}

void geometry::TetraCells(vtkIdType *ip,int nsub,int number){
  if (IPRINT) cout << "geometry::TetraCells"<<endl;

  int i,j,k;
  vtkIdType cell[4],hex[8];

  //--------------------------------------------------------------------
  // Tetra cells for element volume
  //--------------------------------------------------------------------

  for (k=0; k<nsub;     k++){
  for (j=0; j<nsub-k;   j++){
  for (i=0; i<nsub-k-j; i++){

    hex[0] = ip[tetra(nsub,i,  j,  k)];
    hex[1] = ip[tetra(nsub,i+1,j,  k)];
    hex[2] = ip[tetra(nsub,i+1,j+1,k)];
    hex[3] = ip[tetra(nsub,i,  j+1,k)];
    hex[4] = ip[tetra(nsub,i,  j,  k+1)];
    hex[5] = ip[tetra(nsub,i+1,j,  k+1)];
    hex[6] = ip[tetra(nsub,i+1,j+1,k+1)];
    hex[7] = ip[tetra(nsub,i,  j+1,k+1)];

    cell[0] = hex[0];
    cell[1] = hex[1];
    cell[2] = hex[4];
    cell[3] = hex[3];
    vol->InsertNextCell(4,cell);
    volTyp->InsertNextType(VTK_TETRA);
    volElt->InsertNextValue(number);
    volVis->InsertNextValue(1);

    if (i<nsub-k-j-1){

      cell[0] = hex[5];
      cell[1] = hex[1];
      cell[2] = hex[3];
      cell[3] = hex[4];
      vol->InsertNextCell(4,cell);
      volTyp->InsertNextType(VTK_TETRA);
      volElt->InsertNextValue(number);
      volVis->InsertNextValue(1);

      cell[0] = hex[1];
      cell[1] = hex[2];
      cell[2] = hex[5];
      cell[3] = hex[3];
      vol->InsertNextCell(4,cell);
      volTyp->InsertNextType(VTK_TETRA);
      volElt->InsertNextValue(number);
      volVis->InsertNextValue(1);

      cell[0] = hex[7];
      cell[1] = hex[5];
      cell[2] = hex[3];
      cell[3] = hex[4];
      vol->InsertNextCell(4,cell);
      volTyp->InsertNextType(VTK_TETRA);
      volElt->InsertNextValue(number);
      volVis->InsertNextValue(1);

      cell[0] = hex[3];
      cell[1] = hex[7];
      cell[2] = hex[2];
      cell[3] = hex[5];
      vol->InsertNextCell(4,cell);
      volTyp->InsertNextType(VTK_TETRA);
      volElt->InsertNextValue(number);
      volVis->InsertNextValue(1);
    }

    if (i<nsub-k-j-2){
      cell[0] = hex[6];
      cell[1] = hex[5];
      cell[2] = hex[2];
      cell[3] = hex[7];
      vol->InsertNextCell(4,cell);
      volTyp->InsertNextType(VTK_TETRA);
      volElt->InsertNextValue(number);
      volVis->InsertNextValue(1);
    }
  }}}

  //--------------------------------------------------------------------
  // Triangle cells for element faces
  //--------------------------------------------------------------------

  // face 1 (k=0)
  for (j=0; j<nsub; j++){
  for (i=0; i<nsub-j; i++){

    cell[0] = ip[tetra(nsub,i,  j,  0)];
    cell[1] = ip[tetra(nsub,i+1,j,  0)];
    cell[2] = ip[tetra(nsub,i,  j+1,0)];
    sur->InsertNextCell(3,cell);
    surTyp->InsertNextType(VTK_TRIANGLE);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);

    if (i<nsub-j-1){
      cell[0] = ip[tetra(nsub,i+1,j+1,0)];
      cell[1] = ip[tetra(nsub,i,  j+1,0)];
      cell[2] = ip[tetra(nsub,i+1,j,  0)];
      sur->InsertNextCell(3,cell);
      surTyp->InsertNextType(VTK_TRIANGLE);
      surElt->InsertNextValue(number);
      surVis->InsertNextValue(1);
    }
  }}

  // face 2 (i=0)
  for (k=0; k<nsub; k++){
  for (j=0; j<nsub-k; j++){

    cell[0] = ip[tetra(nsub,0,j,  k  )];
    cell[1] = ip[tetra(nsub,0,j+1,k  )];
    cell[2] = ip[tetra(nsub,0,j,  k+1)];
    sur->InsertNextCell(3,cell);
    surTyp->InsertNextType(VTK_TRIANGLE);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);

    if (j<nsub-k-1){
      cell[0] = ip[tetra(nsub,0,j+1,k+1)];
      cell[1] = ip[tetra(nsub,0,j,  k+1)];
      cell[2] = ip[tetra(nsub,0,j+1,k  )];
      sur->InsertNextCell(3,cell);
      surTyp->InsertNextType(VTK_TRIANGLE);
      surElt->InsertNextValue(number);
      surVis->InsertNextValue(1);
    }
  }}

  // face 3 (j=0)
  for (k=0; k<nsub; k++){
  for (i=0; i<nsub-k; i++){

    cell[0] = ip[tetra(nsub,i,  0,k  )];
    cell[1] = ip[tetra(nsub,i+1,0,k  )];
    cell[2] = ip[tetra(nsub,i,  0,k+1)];
    sur->InsertNextCell(3,cell);
    surTyp->InsertNextType(VTK_TRIANGLE);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);

    if (i<nsub-k-1){
      cell[0] = ip[tetra(nsub,i+1,0,k+1)];
      cell[1] = ip[tetra(nsub,i,  0,k+1)];
      cell[2] = ip[tetra(nsub,i+1,0,k  )];
      sur->InsertNextCell(3,cell);
      surTyp->InsertNextType(VTK_TRIANGLE);
      surElt->InsertNextValue(number);
      surVis->InsertNextValue(1);
    }
  }}

  // face 4 (i=nsub-j-k)
  for (k=0; k<nsub; k++){
  for (j=0; j<nsub-k; j++){

    cell[0] = ip[tetra(nsub,nsub-j-k,  j,  k  )];
    cell[1] = ip[tetra(nsub,nsub-j-k-1,j+1,k  )];
    cell[2] = ip[tetra(nsub,nsub-j-k-1,j,  k+1)];
    sur->InsertNextCell(3,cell);
    surTyp->InsertNextType(VTK_TRIANGLE);
    surElt->InsertNextValue(number);
    surVis->InsertNextValue(1);

    if (j<nsub-k-1){
      cell[0] = ip[tetra(nsub,nsub-j-k-2,j+1,k+1)];
      cell[1] = ip[tetra(nsub,nsub-j-k-1,j,  k+1)];
      cell[2] = ip[tetra(nsub,nsub-j-k-1,j+1,k  )];
      sur->InsertNextCell(3,cell);
      surTyp->InsertNextType(VTK_TRIANGLE);
      surElt->InsertNextValue(number);
      surVis->InsertNextValue(1);
    }
  }}

  //--------------------------------------------------------------------
  // Linear cells for element edges
  //--------------------------------------------------------------------

  // edge 1 (j=0, k=0)
  for (i=0; i<nsub; i++){
    cell[0] = ip[tetra(nsub,i,0,0)];
    cell[1] = ip[tetra(nsub,i+1,0,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 2 (i=0, k=0)
  for (j=0; j<nsub; j++){
    cell[0] = ip[tetra(nsub,0,j,0)];
    cell[1] = ip[tetra(nsub,0,j+1,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 3 (i=0, j=0)
  for (k=0; k<nsub; k++){
    cell[0] = ip[tetra(nsub,0,0,k)];
    cell[1] = ip[tetra(nsub,0,0,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 4 (i=nsub-j, k=0)
  for (j=0; j<nsub; j++){
    cell[0] = ip[tetra(nsub,nsub-j,j,0)];
    cell[1] = ip[tetra(nsub,nsub-j-1,j+1,0)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 5 (j=nsub-k, i=0)
  for (k=0; k<nsub; k++){
    cell[0] = ip[tetra(nsub,0,nsub-k,k)];
    cell[1] = ip[tetra(nsub,0,nsub-k-1,k+1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }

  // edge 6 (k=nsub-i, j=0)
  for (i=0; i<nsub; i++){
    cell[0] = ip[tetra(nsub,i,0,nsub-i)];
    cell[1] = ip[tetra(nsub,i+1,0,nsub-i-1)];
    edg->InsertNextCell(2,cell);
    edgElt->InsertNextValue(number);
    edgVis->InsertNextValue(1);
  }
}
