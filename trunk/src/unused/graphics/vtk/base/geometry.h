#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkCellTypes.h"
#include "vtkIntArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkIdList.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkTextActor.h"
#include "vtkPlane.h"
#include "vtkCamera.h"
#include "vtkRenderWindow.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

class geometry{

protected:

  vtkPoints    *points;

  vtkCellArray         *vol;
  vtkCellTypes         *volTyp;
  vtkIntArray          *volElt;
  vtkUnsignedCharArray *volVis;

  vtkCellArray         *sur;
  vtkCellTypes         *surTyp;
  vtkIntArray          *surElt;
  vtkUnsignedCharArray *surVis;

  vtkCellArray         *edg;
  vtkIntArray          *edgElt;
  vtkUnsignedCharArray *edgVis;

  int display_edges;
  int clip_data;
  int interaction_mode;

  vtkTextActor *textActor;
  vtkPlane     *clipPlane;
  vtkCamera    *camera;

public:

  geometry();
  ~geometry();

  int GetInteractionMode(){return interaction_mode;};

  int ReadHexa(ifstream&);
  void HexaCells(vtkIdType*,int,int);

  int ReadTetra(ifstream&);
  void TetraCells(vtkIdType*,int,int);

  int ReadPrism(ifstream&);
  void PrismCells(vtkIdType*,int,int);

  int ReadPyramid(ifstream&);
  void PyramidCells(vtkIdType*,int,int);

  void Read(ifstream&);
  void ReadPoints(vtkIdType*,int,ifstream&);

  void MakeElementInvisible(int);

  void Render();

  void FilterElements(vtkDataSetMapper*);
  void FilterEdges(vtkDataSetMapper*);

  void ExplainKeys();

  void MakeAllElementsVisible(){
    vtkIdType i;
    for (i=0; i<volVis->GetNumberOfTuples(); i++) volVis->SetValue(i,1);
    for (i=0; i<surVis->GetNumberOfTuples(); i++) surVis->SetValue(i,1);
    for (i=0; i<edgVis->GetNumberOfTuples(); i++) edgVis->SetValue(i,1);
  };

  void ToggleClipping(){
    if (clip_data==0){
      clip_data = 1;
    }else{
      clip_data = 0;
    }
  };

  void ToggleInteractionMode(){
    if (interaction_mode==0){
      interaction_mode = 1;
      textActor->SetInput("Click elements to make them invisible ('v' to make all visible)");
    }else{
      interaction_mode = 0;
      textActor->SetInput("Press q to quit");
    }
  };

  void ToggleParallelProjection(){
    if (camera->GetParallelProjection()){
      camera->ParallelProjectionOff();
    }else{
      camera->ParallelProjectionOn();
    }
  };

  void Screenshot(vtkRenderWindow *rw){
    textActor->VisibilityOff();
    vtkWindowToImageFilter *imageFilter = vtkWindowToImageFilter::New();
    imageFilter->SetInput(rw);
    vtkPNGWriter *imageWriter = vtkPNGWriter::New();
    imageWriter->SetInput(imageFilter->GetOutput());
    imageWriter->SetFileName("screenshot.png");
    imageWriter->Write();
    imageFilter->Delete();
    imageWriter->Delete();
    textActor->VisibilityOn();
  };
};

void GetMappers(vtkRenderWindowInteractor*,vtkDataSetMapper*&,vtkDataSetMapper*&);
void LeftButtonPressG(vtkObject*,unsigned long,void*,void*);
void KeyPressG(vtkObject*,unsigned long,void*,void*);

#endif // GEOMETRY_H
