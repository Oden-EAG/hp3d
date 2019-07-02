
#include "solution.h"
#include "util.h"

#include "vtkCellData.h"
#include "vtkTextProperty.h"
#include "vtkRenderer.h"
#include "vtkProperty.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkCellPicker.h"
#include "vtkCallbackCommand.h"
#include "vtkClipDataSet.h"
#include "vtkCutter.h"
#include "vtkAppendFilter.h"
#include "vtkRendererCollection.h"
#include "vtkDataSetMapper.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCellType.h"
#include "vtkColorTransferFunction.h"
#include "vtkPointData.h"

#include <string>
using namespace std;

// Constructor

solution::solution(int complex,int nrcomp){
  if (IPRINT) cout << "solution::solution\n";

  int i;

  this->complex = complex;
  imaginary = 0;
  icomp = 1;

  pointData = vtkDoubleArray::New();
  if (complex){
    pointData->SetNumberOfComponents(2*nrcomp);
    rmin = new double[2*nrcomp];
    rmax = new double[2*nrcomp];
    rmag = new double[nrcomp];

    for (i=0; i<2*nrcomp; i++){ rmin[i] = 1.e10; rmax[i] = -1.e10; }
    for (i=0; i<nrcomp; i++){ rmag[i] = 0.0; }

  }else{
    pointData->SetNumberOfComponents(nrcomp);
    rmin = new double[nrcomp];
    rmax = new double[nrcomp];
    rmag = new double[nrcomp];

    for (i=0; i<nrcomp; i++){ rmin[i] = 1.e10; rmax[i] = -1.e10; rmag[i] = 0.0; }

  }
}

// Destructor

solution::~solution(){
  if (IPRINT) cout << "solution::~solution\n";

  pointData->Delete();
  delete [] rmin;
  delete [] rmax;
  delete [] rmag;
}

// Read from a file

void solution::Read(ifstream &ifs){
  if (IPRINT) cout << "solution::Read"<<endl;

  int nrfig,nfig;
  ifs >> nrfig;
  if (IPRINT) cout << "nrfig = "<<nrfig<<endl;

  int nrpt;
  int nrcomp = pointData->GetNumberOfComponents();
  int n,i,j;
  double *tuple = new double[nrcomp];
  for (int n=0; n<nrfig; n++){

    // Read points and establish cell connectivities
    ifs >> nfig;
    switch(nfig){
    case 1:
      nrpt = ReadHexa(ifs); break;
    case 2:
      nrpt = ReadTetra(ifs); break;
    case 3:
      nrpt = ReadPrism(ifs); break;
    case 4:
      nrpt = ReadPyramid(ifs); break;
    default:
      cout << "geometry::read: nfig = "<<nfig<<endl;
    }

    // Read field data, keeping track of min, max and magnitude of each component
    for (i=0; i<nrpt; i++){
      if (IPRINT) cout << "tuple =";
      for (j=0; j<nrcomp; j++){
        ifs >> tuple[j];
        if (IPRINT) cout << " "<<tuple[j];

        rmin[j] = min(rmin[j],tuple[j]);
        rmax[j] = max(rmax[j],tuple[j]);
        if (complex==0) rmag[j] = max(rmag[j],abs(tuple[j]));

      }
      if (IPRINT) cout << endl;

      if (complex==1){
        for (j=0; j<nrcomp/2; j++){
          rmag[j] = max(rmag[j],sqrt(tuple[2*j]*tuple[2*j]+tuple[2*j+1]*tuple[2*j+1]));
	}
      }

      pointData->InsertNextTuple(tuple);
    }
  }
  delete [] tuple;
}


void solution::RescaleComponent(){
  if (complex){
    cout << "REAL MIN = "<<rmin[2*(icomp-1)]<<", MAX = "<<rmax[2*(icomp-1)]<<endl;
    cout << "IMAG MIN = "<<rmin[2*icomp-1]<<", MAX = "<<rmax[2*icomp-1]<<endl;
    cout << "MAX MAGNITUDE = "<<rmag[icomp-1]<<endl;
    cout << "ENTER NEW REAL MIN:"<<endl;
    cin >> rmin[2*(icomp-1)];
    cout << "ENTER NEW REAL MAX:"<<endl;
    cin >> rmax[2*(icomp-1)];
    cout << "ENTER NEW IMAG MIN:"<<endl;
    cin >> rmin[2*icomp-1];
    cout << "ENTER NEW IMAG MAX:"<<endl;
    cin >> rmax[2*icomp-1];
    cout << "ENTER NEW MAX MAGNITUDE:"<<endl;
    cin >> rmag[icomp-1];
  }else{
    cout << "MIN = "<<rmin[icomp-1]<<", MAX = "<<rmax[icomp-1]<<endl;
    cout << "ENTER NEW MIN:"<<endl;
    cin >> rmin[icomp-1];
    cout << "ENTER NEW MAX:"<<endl;
    cin >> rmax[icomp-1];
  }
}

// Render

void solution::Render(){
  if (IPRINT) cout << "solution::Render"<<endl;

  //----------------------------//
  // CREATE RENDERER AND ACTORS //
  //----------------------------//

  // Create renderer
  vtkRenderer *ren = vtkRenderer::New();
  if (IPRINT) ren->DebugOn();
  ren->SetActiveCamera(camera);
  ren->SetBackground(1,1,1);
  ren->TwoSidedLightingOn();

  // Map from scalars to colors
  vtkColorTransferFunction *colorMap = vtkColorTransferFunction::New();
  colorMap->SetColorSpaceToRGB();
  colorMap->ClampingOn();

  // ELEMENT ACTOR

  // Create element mapper
  vtkDataSetMapper *elementMapper = vtkDataSetMapper::New();
  if (IPRINT) elementMapper->DebugOn();
  elementMapper->SetResolveCoincidentTopologyToPolygonOffset();
  elementMapper->ScalarVisibilityOn();
  elementMapper->SetScalarModeToUsePointData();
  elementMapper->SetLookupTable(colorMap);
  elementMapper->UseLookupTableScalarRangeOn();
  FilterElements(elementMapper);

  // Create element actor and add it to the renderer
  vtkActor *elementActor = vtkActor::New();
  if (IPRINT) elementActor->DebugOn();
  elementActor->SetMapper(elementMapper);
  //elementActor->GetProperty()->SetInterpolationToGouraud();  // Flat,Gouraud,Phong
  elementActor->GetProperty()->SetRepresentationToSurface(); // Points,Wireframe,Surface
  elementActor->GetProperty()->SetAmbient(1.0);
  elementActor->GetProperty()->SetDiffuse(0.0);
  elementActor->GetProperty()->SetSpecular(0.0);
  elementActor->GetProperty()->BackfaceCullingOff();
  elementActor->PickableOn();
  elementActor->VisibilityOn();
  ren->AddActor(elementActor);

  // EDGE ACTOR

  // Create edge mapper
  vtkDataSetMapper *edgeMapper = vtkDataSetMapper::New();
  if (IPRINT) edgeMapper->DebugOn();
  edgeMapper->ScalarVisibilityOff();
  FilterEdges(edgeMapper);

  // Create edge actor and add it to the renderer
  vtkActor *edgeActor = vtkActor::New();
  if (IPRINT) edgeActor->DebugOn();
  edgeActor->SetMapper(edgeMapper);
  edgeActor->GetProperty()->SetColor(0,0,0);
  edgeActor->GetProperty()->BackfaceCullingOn();
  edgeActor->PickableOff();
  edgeActor->SetVisibility(display_edges);
  ren->AddActor(edgeActor);

  ren->AddActor2D(textActor);

  ren->ResetCamera();

  //----------------------------------------//
  // CREATE RENDERING WINDOW AND INTERACTOR //
  //----------------------------------------//

  // Create rendering window and add renderer
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  if (IPRINT) renWin->DebugOn();
  renWin->AddRenderer(ren);
  renWin->SetSize(800,600);
  renWin->SwapBuffersOn();

  // Create an interactor and add the rendering window
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  if (IPRINT) iren->DebugOn();
  iren->SetRenderWindow(renWin);

  // Set the style of interaction
  vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
  if (IPRINT) style->DebugOn();
  iren->SetInteractorStyle(style);

  // Set the type of picker to be used
  vtkCellPicker *picker = vtkCellPicker::New();
  if (IPRINT) picker->DebugOn();
  picker->SetTolerance(0.001);
  iren->SetPicker(picker);

  // Set subroutine to be called on left button press event
  vtkCallbackCommand *leftButtonPressCmd = vtkCallbackCommand::New();
  leftButtonPressCmd->SetCallback(LeftButtonPressS);
  leftButtonPressCmd->SetClientData((void*)this);
  iren->AddObserver(vtkCommand::LeftButtonPressEvent,leftButtonPressCmd);

  // Set subroutine to be called on key press event
  vtkCallbackCommand *keyPressCmd = vtkCallbackCommand::New();
  keyPressCmd->SetCallback(KeyPressS);
  keyPressCmd->SetClientData((void*)this);
  iren->AddObserver(vtkCommand::KeyPressEvent,keyPressCmd);

  ExplainKeys();

  //------------------------//
  // START INTERACTIVE MODE //
  //------------------------//

  iren->Initialize();
  iren->Start();

  //----------//
  // CLEAN UP //
  //----------//

  // Render window and interactor
  renWin->Delete();
  iren->Delete();
  style->Delete();
  picker->Delete();
  leftButtonPressCmd->Delete();
  keyPressCmd->Delete();

  // Renderer and Actors
  ren->Delete();
  elementMapper->Delete(); elementActor->Delete();
  edgeMapper->Delete();    edgeActor->Delete();
  colorMap->Delete();
}


void solution::FilterElements(vtkDataSetMapper *elementMapper){
  if (IPRINT) cout << "solution::FilterElements"<<endl;

  vtkIdType i,*cell;
  unsigned char type;

  // Extract scalar point data
  vtkIdType nrpt = points->GetNumberOfPoints();
  vtkDoubleArray *pointScalars = vtkDoubleArray::New();
  pointScalars->SetNumberOfValues(nrpt);
  if (complex){
    if (imaginary==0 || imaginary==1){
      for (i=0; i<nrpt; i++)
        pointScalars->SetValue(i,pointData->GetComponent(i,2*(icomp-1)+imaginary));
    }else{
      for (i=0; i<nrpt; i++)
        pointScalars->SetValue(i,mag(pointData->GetComponent(i,2*(icomp-1)),
                                     pointData->GetComponent(i,2*icomp-1)));
    }
  }else{
    for (i=0; i<nrpt; i++)
      pointScalars->SetValue(i,pointData->GetComponent(i,icomp));
  }

  // Modify the color map
  if (complex){
    vtkColorTransferFunction *colorMap = (vtkColorTransferFunction*)elementMapper->GetLookupTable();
    colorMap->RemoveAllPoints();
    if (imaginary==0 || imaginary==1){
      colorMap->AddRGBPoint(rmin[2*(icomp-1)+imaginary],                                   0.0,0.0,1.0); // Blue  = min
      colorMap->AddRGBPoint((rmin[2*(icomp-1)+imaginary]+rmax[2*(icomp-1)+imaginary])/2.0, 0.0,1.0,0.0); // Green = mid
      colorMap->AddRGBPoint(rmax[2*(icomp-1)+imaginary],                                   1.0,0.0,0.0); // Red   = max
    }else{
      colorMap->AddRGBPoint(0.0,           0.0,1.0,0.0); // Green = zero
      colorMap->AddRGBPoint(rmag[icomp-1], 1.0,0.0,0.0); // Red   = max
    }
  }

  // ELEMENTS

  // Create unstructured grid of visible elements
  vtkUnstructuredGrid *elements = vtkUnstructuredGrid::New();
  vtkIntArray *elementNumbers   = vtkIntArray::New();
  elements->SetPoints(points);
  elements->GetPointData()->SetScalars(pointScalars);
  vol->InitTraversal();
  for (i=0; i<vol->GetNumberOfCells(); i++){
    vol->GetNextCell(nrpt,cell);
    type = volTyp->GetCellType(i);
    if (volVis->GetValue(i)){
      elements->InsertNextCell(type,nrpt,cell);
      elementNumbers->InsertNextValue(volElt->GetValue(i));
    }
  }
  elements->GetCellData()->SetScalars(elementNumbers);
  elementNumbers->Delete();
  pointScalars->Delete();

  // Pass visible elements (possibly after clipping) to element mapper
  if (clip_data==1){
    vtkClipDataSet *elementClipper = vtkClipDataSet::New();
    if (IPRINT) elementClipper->DebugOn();
    elementClipper->SetClipFunction(clipPlane);
    elementClipper->InsideOutOn();
    elementClipper->GenerateClipScalarsOff();
    elementClipper->SetInput(elements);
    elementMapper->SetInput(elementClipper->GetOutput());
    elementClipper->Delete();
  }else{
    elementMapper->SetInput(elements);
  }
  elements->Delete();
}

void solution::Animate(vtkRenderWindowInteractor *iren){
  if (IPRINT) cout << "solution::Animate"<<endl;

  system("mkdir frames");

  vtkIdType i,*cell;
  unsigned char type;

  // Get element and edge mappers
  vtkDataSetMapper *elementMapper,*edgeMapper;
  GetMappers(iren, elementMapper,edgeMapper);

  // Create scalar point data array
  vtkIdType nrpt = points->GetNumberOfPoints();
  vtkDoubleArray *pointScalars = vtkDoubleArray::New();
  pointScalars->SetNumberOfValues(nrpt);

  // Create unstructured grid of visible elements
  vtkUnstructuredGrid *elements = vtkUnstructuredGrid::New();
  elements->SetPoints(points);
  elements->GetPointData()->SetScalars(pointScalars);
  vol->InitTraversal();
  for (i=0; i<vol->GetNumberOfCells(); i++){
    vol->GetNextCell(nrpt,cell);
    type = volTyp->GetCellType(i);
    if (volVis->GetValue(i)) elements->InsertNextCell(type,nrpt,cell);
  }

  // Pass visible elements (possibly after clipping) to element mapper
  vtkClipDataSet *elementClipper;
  if (clip_data){
    elementClipper = vtkClipDataSet::New();
    if (IPRINT) elementClipper->DebugOn();
    elementClipper->SetClipFunction(clipPlane);
    elementClipper->InsideOutOn();
    elementClipper->GenerateClipScalarsOff();
    elementClipper->SetInput(elements);
    elementMapper->SetInput(elementClipper->GetOutput());
  }else{
    elementMapper->SetInput(elements);
  }

  // Modify the color map
  vtkColorTransferFunction *colorMap = (vtkColorTransferFunction*)elementMapper->GetLookupTable();
  colorMap->RemoveAllPoints();
  colorMap->AddRGBPoint(-rmag[icomp-1], 0.0,0.0,1.0); // Blue  = min
  colorMap->AddRGBPoint(0.0,            0.0,1.0,0.0); // Green = zero
  colorMap->AddRGBPoint(rmag[icomp-1],  1.0,0.0,0.0); // Red   = max

  // Retrieve render window and setup PNG writer
  vtkRenderWindow *renWin = iren->GetRenderWindow();
  vtkWindowToImageFilter *imageFilter = vtkWindowToImageFilter::New();
  imageFilter->SetInput(renWin);
  vtkPNGWriter *writer = vtkPNGWriter::New();
  writer->SetInput(imageFilter->GetOutput());

  // Loop through frames
  int iframe,nrframes = 36;
  double pi = 4.0*atan(1.0), phi,cosphi,sinphi;
  string str="frames/image00.png";
  string digits="0123456789";
  int ten,one;
  for (iframe=0; iframe<nrframes; iframe++){
    phi = (2.0*pi*iframe)/nrframes;
    cosphi = cos(phi);
    sinphi = sin(phi);

    // Update scalar point data
    nrpt = points->GetNumberOfPoints();
    for (i=0; i<nrpt; i++)
      pointScalars->SetValue(i,cosphi*pointData->GetComponent(i,2*icomp-2)
                             - sinphi*pointData->GetComponent(i,2*icomp-1));

    pointScalars->Modified();

    renWin->Render();

    imageFilter->Modified();

    ten = iframe/10;
    one = iframe - 10*ten;
    str.replace(12,1,digits,ten,1);
    str.replace(13,1,digits,one,1);
    cout << "Writing "<<str<<" ..."<<endl;
    writer->SetFileName(str.c_str());

    writer->Write();
  }

  writer->Delete();
  imageFilter->Delete();

  if (clip_data) elementClipper->Delete();
  elements->Delete();
  pointScalars->Delete();
}

void LeftButtonPressS(vtkObject *obj,
                      unsigned long eid,
                      void *clientdata,
                      void *calldata){
  vtkRenderWindowInteractor *iren = (vtkRenderWindowInteractor*)obj;
  solution *cd = (solution*)clientdata;

  int ix,iy,res;
  vtkIdType cell_id;
  vtkIntArray *ia;
  int iel;

  switch(cd->GetInteractionMode()){
  case 1:
    //
    // if element was clicked, make it invisible
    //
    // get mouse coordinates from render window interactor
    iren->GetMousePosition(&ix,&iy);
    if (IPRINT) cout << "LeftButtonPress: ix,iy = "<<ix<<" "<<iy<<endl;

    // get picker
    vtkCellPicker *picker = (vtkCellPicker*)iren->GetPicker();

    // use picker to pick mouse coordinates
    res = picker->Pick(ix,iy,0,iren->FindPokedRenderer(ix,iy));
    if (IPRINT) cout << "LeftButtonPress: res = "<<res<<endl;

    // return if nothing was picked
    if (res==0) return;

    // get the id of the picked cell
    cell_id = picker->GetCellId();
    if (IPRINT) cout << "LeftButtonPress: cell_id = "<<cell_id<<endl;

    // get associated integer data (element number)
    ia  = (vtkIntArray*)picker->GetDataSet()->GetCellData()->GetScalars();
    iel = ia->GetValue(cell_id);
    if (IPRINT) cout << "LeftButtonPress: iel = "<<iel<<endl;

    cd->MakeElementInvisible(iel);

    vtkDataSetMapper *elementMapper,*edgeMapper;
    GetMappers(iren, elementMapper,edgeMapper);
    cd->FilterElements(elementMapper);
    cd->FilterEdges(edgeMapper);

    iren->Render();
    break;
  }
}

void KeyPressS(vtkObject *obj,
               unsigned long eid,
               void* clientdata,
               void *calldata){
  vtkRenderWindowInteractor *iren = (vtkRenderWindowInteractor*)obj;
  solution *cd = (solution*)clientdata;

  vtkDataSetMapper       *elementMapper,*edgeMapper;
  vtkWindowToImageFilter *imageFilter;
  vtkPNGWriter           *imageWriter;

  char key = iren->GetKeyCode();
  switch(key){
  case '+':
    //
    // increment component to be displayed
    //
    cd->IncrementIcomp();
    GetMappers(iren, elementMapper,edgeMapper);
    cd->FilterElements(elementMapper);
    iren->Render();
    break;
  case '-':
    //
    // decrement component to be displayed
    //
    cd->DecrementIcomp();
    GetMappers(iren, elementMapper,edgeMapper);
    cd->FilterElements(elementMapper);
    iren->Render();
    break;
  case 'a':
    //
    // animate
    //
    cd->Animate(iren);
    GetMappers(iren, elementMapper,edgeMapper);
    cd->FilterElements(elementMapper);
    iren->Render();
    break;
  case 'c':
    //
    // toggle Clipping of data
    //
    cd->ToggleClipping();
    GetMappers(iren, elementMapper,edgeMapper);
    cd->FilterElements(elementMapper);
    cd->FilterEdges(edgeMapper);
    iren->Render();
    break;
  case 'd':
    //
    // toggle between real, imaginary and magnitude
    //
    cd->ToggleImaginary();
    GetMappers(iren, elementMapper,edgeMapper);
    cd->FilterElements(elementMapper);
    iren->Render();
    break;
  case 'i':
    //
    // interactively make elements Invisible
    //
    cd->ToggleInteractionMode();
    iren->Render();
    break;
  case 's':
    //
    // change scaling
    //
    cd->RescaleComponent();
    GetMappers(iren, elementMapper,edgeMapper);
    cd->FilterElements(elementMapper);
    iren->Render();
    break;
  case 't':
    //
    // Take a screenshot
    //
    cd->Screenshot(iren->GetRenderWindow());
    iren->Render();
    break;
  case 'v':
    //
    // make all elements Visible
    //
    cd->MakeAllElementsVisible();
    GetMappers(iren, elementMapper,edgeMapper);
    cd->FilterElements(elementMapper);
    cd->FilterEdges(edgeMapper);
    iren->Render();
    break;
  }
}

void solution::ExplainKeys(){
  cout << "Explanation of keys:"<<endl;
  cout << "+  Increments component being displayed"<<endl;
  cout << "-  Decrements component being displayed"<<endl;
  cout << "a  Animate current component"<<endl;
  cout << "c  Toggle clipping of data"<<endl;
  cout << "d  Toggle display between real, imaginary and magnitude"<<endl;
  cout << "e  Exit (same as Quit)"<<endl;
  cout << "i  Toggle mode to make elements invisible"<<endl;
  cout << "q  Quit (same as Exit)"<<endl;
  cout << "r  Recenter object in middle of window"<<endl;
  cout << "s  Change default scalar range"<<endl;
  cout << "t  Take a screenshot (screenshot.png)"<<endl;
  cout << "v  Make all elements visible"<<endl;
}
