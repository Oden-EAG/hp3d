
#include "geometry.h"
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

#include <time.h>

// Constructor

geometry::geometry(){

  if (IPRINT) cout << "geometry::geometry\n";

  points = vtkPoints::New();

  vol    = vtkCellArray::New();
  volTyp = vtkCellTypes::New();
  volElt = vtkIntArray::New();
  volVis = vtkUnsignedCharArray::New();

  sur    = vtkCellArray::New();
  surTyp = vtkCellTypes::New();
  surElt = vtkIntArray::New();
  surVis = vtkUnsignedCharArray::New();

  edg    = vtkCellArray::New();
  edgElt = vtkIntArray::New();
  edgVis = vtkUnsignedCharArray::New();

  display_edges = 1;
  clip_data = 1;
  interaction_mode = 0;

  textActor = vtkTextActor::New();
  textActor->SetInput("Press q to quit");
  textActor->SetTextScaleModeToNone();
  textActor->SetPosition(5,5);
  textActor->PickableOff();
  textActor->GetTextProperty()->SetColor(0,0,0);

  clipPlane = vtkPlane::New();
  clipPlane->SetOrigin(0.0,0.0,0.0);
  clipPlane->SetNormal(0.0,1.0,0.0);

  camera = vtkCamera::New();
  camera->SetPosition(8.0,8.0,8.0);
  camera->SetFocalPoint(0.0,0.0,0.0);
  camera->SetViewUp(0.0,0.0,1.0);

}

// Destructor

geometry::~geometry(){

  if (IPRINT) cout << "geometry::~geometry\n";

  points->Delete();

  vol->Delete();
  volTyp->Delete();
  volElt->Delete();
  volVis->Delete();

  sur->Delete();
  surTyp->Delete();
  surElt->Delete();
  surVis->Delete();

  edg->Delete();
  edgElt->Delete();
  edgVis->Delete();

  textActor->Delete();

  clipPlane->Delete();

  camera->Delete();

}

// Read from file

void geometry::Read(ifstream &ifs){
  clock_t c1 = clock();
  if (IPRINT) cout << "geometry::Read"<<endl;

  int nrfig,ifig;
  ifs >> nrfig;
  if (IPRINT) cout << "nrfig = "<<nrfig<<endl;

  for (int i=0; i<nrfig; i++){
    ifs >> ifig;
    if (IPRINT) cout << "ifig = "<<ifig<<endl;

    switch(ifig){
    case 1:
      ReadHexa(ifs); break;
    case 2:
      ReadTetra(ifs); break;
    case 3:
      ReadPrism(ifs); break;
    case 4:
      ReadPyramid(ifs); break;
    default:
      cout << "geometry::read: ifig = "<<ifig<<endl;
    }
  }
  cout << "Read: points->GetNumberOfPoints()="<<points->GetNumberOfPoints()<<endl;
  cout << "Read: vol->GetNumberOfCells()    ="<<vol->GetNumberOfCells()<<endl;
  cout << "Read: sur->GetNumberOfCells()    ="<<sur->GetNumberOfCells()<<endl;
  cout << "Read: edg->GetNumberOfCells()    ="<<edg->GetNumberOfCells()<<endl;
  clock_t c2 = clock();
  cout << "Read: time = "<<double(c2-c1)/double(CLOCKS_PER_SEC)<<" s"<<endl;
}
//Read: points->GetNumberOfPoints()=135388
//Read: vol->GetNumberOfCells()    =101044
//Read: sur->GetNumberOfCells()    =206592
//Read: edg->GetNumberOfCells()    =158372
//Read: time = 0.74 s

void geometry::ReadPoints(vtkIdType *ip,int nrpt,ifstream &ifs){
  if (IPRINT) cout << "geometry::ReadPoints"<<endl;

  // Read points and record their global locations
  double x[3];
  for (int i=0; i<nrpt; i++){
    ifs >> x[0] >> x[1] >> x[2];
    if (IPRINT) cout << "x["<<i<<"]=("<<x[0]<<","<<x[1]<<","<<x[2]<<")"<<endl;
    ip[i] = points->InsertNextPoint(x);
  }
}

void geometry::MakeElementInvisible(int iel){
  if (IPRINT) cout << "geometry::MakeElementInvisible("<<iel<<")"<<endl;
  vtkIdType i;

  for (i=0; i<volElt->GetNumberOfTuples(); i++)
    if (volElt->GetValue(i)==iel) volVis->SetValue(i,0);

  for (i=0; i<surElt->GetNumberOfTuples(); i++)
    if (surElt->GetValue(i)==iel) surVis->SetValue(i,0);

  for (i=0; i<edgElt->GetNumberOfTuples(); i++)
    if (edgElt->GetValue(i)==iel) edgVis->SetValue(i,0);
}

////////////////////////////////////////////////////////////////////////
///////////////////////// RENDER GEOMETRY //////////////////////////////
////////////////////////////////////////////////////////////////////////

void geometry::Render(){
  if (IPRINT) cout << "geometry::Render"<<endl;

  ExplainKeys();

  //----------------------------//
  // CREATE RENDERER AND ACTORS //
  //----------------------------//

  // Create renderer
  vtkRenderer *ren = vtkRenderer::New();
  if (IPRINT) ren->DebugOn();
  ren->SetActiveCamera(camera);
  ren->SetBackground(1,1,1);
  ren->TwoSidedLightingOn();

  // ELEMENT ACTOR

  // Create element mapper
  vtkDataSetMapper *elementMapper = vtkDataSetMapper::New();
  if (IPRINT) elementMapper->DebugOn();
  elementMapper->SetResolveCoincidentTopologyToPolygonOffset();
  elementMapper->ScalarVisibilityOff();
  FilterElements(elementMapper);

  // Create element actor and add it to the renderer
  vtkActor *elementActor = vtkActor::New();
  if (IPRINT) elementActor->DebugOn();
  elementActor->SetMapper(elementMapper);
  elementActor->GetProperty()->SetInterpolationToFlat();     // Flat,Gouraud,Phong
  elementActor->GetProperty()->SetRepresentationToSurface(); // Points,Wireframe,Surface
  elementActor->GetProperty()->SetColor(0,0,1);
  elementActor->GetProperty()->SetAmbient(0.5);
  elementActor->GetProperty()->SetDiffuse(0.6);
  elementActor->GetProperty()->SetSpecular(1.0);
  elementActor->GetProperty()->SetSpecularPower(10.0);
  elementActor->GetProperty()->SetOpacity(1.0);
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
  leftButtonPressCmd->SetCallback(LeftButtonPressG);
  leftButtonPressCmd->SetClientData((void*)this);
  iren->AddObserver(vtkCommand::LeftButtonPressEvent,leftButtonPressCmd);

  // Set subroutine to be called on key press event
  vtkCallbackCommand *keyPressCmd = vtkCallbackCommand::New();
  keyPressCmd->SetCallback(KeyPressG);
  keyPressCmd->SetClientData((void*)this);
  iren->AddObserver(vtkCommand::KeyPressEvent,keyPressCmd);

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
}

void geometry::FilterElements(vtkDataSetMapper *elementMapper){
  if (IPRINT) cout << "geometry::FilterElements"<<endl;

  vtkIdType i,nrpt,*cell;
  unsigned char type;

  // ELEMENTS

  // Create unstructured grid of visible elements
  vtkUnstructuredGrid *elements = vtkUnstructuredGrid::New();
  vtkIntArray *elementNumbers   = vtkIntArray::New();
  elements->SetPoints(points);
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

void geometry::FilterEdges(vtkDataSetMapper *edgeMapper){
  if (IPRINT) cout << "geometry::FilterEdges"<<endl;

  vtkIdType i,nrpt,*cell;
  unsigned char type;

  // EDGES

  // Create polydata of visible edges
  vtkUnstructuredGrid *edges = vtkUnstructuredGrid::New();
  if (IPRINT) edges->DebugOn();
  edges->SetPoints(points);
  edg->InitTraversal();
  for (i=0; i<edg->GetNumberOfCells(); i++){
    edg->GetNextCell(nrpt,cell);
    if (edgVis->GetValue(i)) edges->InsertNextCell(VTK_LINE,2,cell);
  }

  // Pass visible edges (possibly after clipping) to edge mapper
  if (clip_data==1){
    vtkClipDataSet *edgeClipper = vtkClipDataSet::New();
    if (IPRINT) edgeClipper->DebugOn();
    edgeClipper->SetClipFunction(clipPlane);
    edgeClipper->InsideOutOn();
    edgeClipper->GenerateClipScalarsOff();
    edgeClipper->SetInput(edges);

    // Create polydata of visible faces
    vtkUnstructuredGrid *faces = vtkUnstructuredGrid::New();
    faces->SetPoints(points);
    sur->InitTraversal();
    for (i=0; i<sur->GetNumberOfCells(); i++){
      sur->GetNextCell(nrpt,cell);
      type = surTyp->GetCellType(i);
      if (surVis->GetValue(i)) faces->InsertNextCell(type,nrpt,cell);
    }

    // Cut through visible faces to obtain cut-edges and pass them to the cut-edge mapper
    vtkCutter *faceCutter = vtkCutter::New();
    if (IPRINT) faceCutter->DebugOn();
    faceCutter->SetCutFunction(clipPlane);
    faceCutter->GenerateCutScalarsOff();
    faceCutter->SetInput(faces);

    vtkAppendFilter *append = vtkAppendFilter::New();
    append->SetInput(edgeClipper->GetOutput());
    append->AddInput(faceCutter->GetOutput());
    edgeMapper->SetInput(append->GetOutput());

    append->Delete();
    faceCutter->Delete();
    faces->Delete();
    edgeClipper->Delete();
  }else{
    edgeMapper->SetInput(edges);
  }

  edges->Delete();
}

void geometry::ExplainKeys(){
  cout << "Explanation of keys:"<<endl;
  cout << "c - toggle Clipping of data"<<endl;
  cout << "e - Exit (same as Quit)"<<endl;
  cout << "i - make elements Invisible"<<endl;
  cout << "l - toggle Look between projection and perspective"<<endl;
  cout << "q - Quit (same as Exit)"<<endl;
  cout << "r - Recenter object in middle of window"<<endl;
  cout << "t - Take a screenshot (screenshot.png)"<<endl;
  cout << "v - make all elements Visible"<<endl;
}

void GetMappers(vtkRenderWindowInteractor *iren,
                vtkDataSetMapper *&elementMapper,
                vtkDataSetMapper *&edgeMapper){
  vtkRenderWindow       *renWin = iren->GetRenderWindow();
  vtkRendererCollection *renCol = renWin->GetRenderers();
  vtkRenderer           *ren    = renCol->GetFirstRenderer();
  vtkActorCollection    *actCol = ren->GetActors();
  actCol->InitTraversal();
  elementMapper = (vtkDataSetMapper*)actCol->GetNextActor()->GetMapper();
  edgeMapper    = (vtkDataSetMapper*)actCol->GetNextActor()->GetMapper();
}

void LeftButtonPressG(vtkObject *obj,
                      unsigned long eid,
                      void *clientdata,
                      void *calldata){
  vtkRenderWindowInteractor *iren = (vtkRenderWindowInteractor*)obj;
  geometry *cd = (geometry*)clientdata;

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

void KeyPressG(vtkObject *obj,
               unsigned long eid,
               void* clientdata,
               void *calldata){
  vtkRenderWindowInteractor *iren = (vtkRenderWindowInteractor*)obj;
  geometry *cd = (geometry*)clientdata;

  vtkDataSetMapper       *elementMapper,*edgeMapper;
  vtkWindowToImageFilter *imageFilter;
  vtkPNGWriter           *imageWriter;

  char key = iren->GetKeyCode();
  switch(key){
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

  case 'i':
    //
    // interactively make elements Invisible
    //
    cd->ToggleInteractionMode();
    iren->Render();
    break;

  case 'l':
    //
    // toggle camera projection
    //
    cd->ToggleParallelProjection();
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
