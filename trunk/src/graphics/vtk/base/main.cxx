
#include "geometry.h"
#include "solution.h"

int main(){

  ifstream ifs;
  geometry *geom;
  solution *soln;
  int idec,icomplex,nrcomp;

  idec = 1;
  while (idec>0){
    cout << "SELECT:"<<endl;
    cout << "QUIT................ 0"<<endl;
    cout << "GEOMETRY GRAPHICS... 1"<<endl;
    cout << "SOLUTION GRAPHICS... 2"<<endl;
    cin >> idec;

    switch(idec){
    case 1:

      ifs.open("geom",ifstream::in);
      geom = new geometry();
      geom->Read(ifs);
      ifs.close();

      geom->Render();

      delete geom;

      break;

    case 2:

      ifs.open("soln",ifstream::in);
      ifs >> icomplex>>nrcomp;
      //cout << "icomplex = "<<icomplex<<endl;
      //cout << "nrcomp   = "<<nrcomp<<endl;
      soln = new solution(icomplex,nrcomp);
      soln->Read(ifs);
      ifs.close();

      soln->Render();

      delete soln;

      break;

    }
  }

  return 0;
}
