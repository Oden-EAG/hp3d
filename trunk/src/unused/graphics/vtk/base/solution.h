#ifndef SOLUTION_H
#define SOLUTION_H

#include "geometry.h"

#include "vtkDoubleArray.h"

class solution: public geometry{

protected:

  int complex;                // 0 - real data, 1 - complex data
  int imaginary;              // 0 - display real part, 1 - display imaginary part
  int icomp;                  // real or complex valued component to be displayed
  vtkDoubleArray *pointData;
  double *rmin;
  double *rmax;
  double *rmag;

public:

  solution(int,int);
  ~solution();

  void Read(ifstream&);

  void Render();

  void FilterElements(vtkDataSetMapper *);

  void ToggleImaginary(void){
    imaginary++;
    if (imaginary==3) imaginary = 0;
  };

  void IncrementIcomp(){
    int nrcomp = pointData->GetNumberOfComponents();
    if (complex) nrcomp = nrcomp/2;
    icomp=icomp+1;
    if (icomp>nrcomp) icomp=1;
  };

  void DecrementIcomp(){
    int nrcomp = pointData->GetNumberOfComponents();
    if (complex) nrcomp = nrcomp/2;
    icomp=icomp-1;
    if (icomp<1) icomp=nrcomp;
  };

  void RescaleComponent();

  void Animate(vtkRenderWindowInteractor *);

  void ExplainKeys();
};

void LeftButtonPressS(vtkObject*,unsigned long,void*,void*);
void KeyPressS(vtkObject*,unsigned long,void*,void*);

#endif // SOLUTION_H
