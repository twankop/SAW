
/* Title: Self avoiding random walks
 * Course: Computational Physics 2015
 * By: Twan Koperberg
 * Student ID: 0713309
 * Leiden University
 */

/* Notes:
 * 
 * SIZE gives the number of monomers, meaning that the polymer has length, SIZE-1.
 * 
 * Monomer identifiers range from 1 to SIZE,
 * and the coordinates for monomer i are coordinates[i-1].
 * The coordinate values have no bounds. (They could possibly overflow!)
 * 
 * Lattice coordinates range from 0 to SIZE-1,
 * and lattice[i][j] contains the site with coordinates x=i, y=j.
 *
 * ALGORITHM_CHOICE==1 is the local moves algorithm.
 * ALGORITHM_CHOICE==2 is the slithering snake algorithm.
 * ALGORITHM_CHOICE==3 is the pivot algorithm. 
 */

#include <iostream> //cout
#include <fstream>  //instream, ofstream
#include <cstdlib>  //rand, srand
#include <ctime>    //time_t, time
#include <cmath>    //abs, exp, pow
//#include <cstdio>
//#include "SFMT/SFMT.o"

using namespace std; 

const int SIZE=11;
const int SAMPLE_SIZE=100000;
const int ALGORITHM_CHOICE=3;
const int WITH_THERMALIZATION_OUTPUT=0;
const int MEASURE_FREQUENCY=pow(SIZE,0);
const string thermalizationFilename="thermalization_data_part?_S=?_A=?.txt";
const string dataFilename="data_S=11_A=3.txt";

class Coordinates{
  friend class RandomWalk;
  public:
    Coordinates();
  //private:
    int x;
    int y;
};//Coordinates

class RandomWalk{
  public:
    RandomWalk();
    void printLattice();
    void printCoordinates();
    bool isSelfAvoiding();
    double computeEndToEndLength();
    double computeGyrationRadius();
    int computeVerticalityDegree();
    void moveLocally(int movingPoint);
    void slither();
    void pivot(int pivotPoint, int pivotDirection);
    void doTimestep();
    void thermalize();
  //private:
    int lattice[SIZE][SIZE];
    Coordinates coordinates[SIZE];
    bool isOneHead;
};//RandomWalk

int doMod(int x, int m) {
    return ((x%=m)<0)?x+m:x;
}//doMod

Coordinates::Coordinates(){
  x=0;
  y=0;
}//Coordinates::constructor

RandomWalk::RandomWalk(){
  for(int i=0;i<SIZE;i++){
    lattice[i][0]=i+1;
    coordinates[i].x=i;
    coordinates[i].y=0;
  }
  for(int j=1;j<SIZE;j++){
    for(int i=0;i<SIZE;i++){
      lattice[i][j]=0;
    }
  }
  isOneHead=true;
}//RandomWalk::constructor

void RandomWalk::printLattice(){
  ofstream output("random_walk_lattice.txt");
  for(int j=SIZE-1;j>=0;j--){
    for(int i=0;i<SIZE;i++){
      output <<lattice[i][j]<<" ";
    }
    output <<"\n";
  }
  output.close();
}//RandomWalk::printLattice

void RandomWalk::printCoordinates(){
  ofstream output("random_walk_coordinates.txt");
  for(int i=0;i<SIZE;i++){
    output <<"("<<coordinates[i].x<<","<<coordinates[i].y<<")\n";
  }
  output.close();
}//RandomWalk::printLattice

bool RandomWalk::isSelfAvoiding(){
  bool isSelfAvoiding=true;
  for(int i=0;isSelfAvoiding and i<SIZE-1;i++){
    if(lattice[doMod(coordinates[i].x,SIZE)]
              [doMod(coordinates[i].y,SIZE)]!=i+1){
      isSelfAvoiding=false;
    }
  }
  return isSelfAvoiding;
}//RandomWalk::isSelfAvoiding

double RandomWalk::computeEndToEndLength(){
  return pow(coordinates[0].x-coordinates[SIZE-1].x,2)+
         pow(coordinates[0].y-coordinates[SIZE-1].y,2);
}//RandomWalk::computeEndToEndLength

double RandomWalk::computeGyrationRadius(){
  Coordinates centreOfMass;
  double radiusOfGyration=0;
  
  for(int i=0;i<SIZE;i++){
    centreOfMass.x+=coordinates[i].x;
    centreOfMass.y+=coordinates[i].y;
  }
  centreOfMass.x/=SIZE;
  centreOfMass.y/=SIZE;
  
  for(int i=0;i<SIZE;i++){
    radiusOfGyration+=pow(coordinates[i].x-centreOfMass.x,2)+
                      pow(coordinates[i].y-centreOfMass.y,2);
  }

  return radiusOfGyration/SIZE;
}//RandomWalk::computeGyrationRadius

int RandomWalk::computeVerticalityDegree(){
  int verticalityDegree=0;
  for(int i=0;i<SIZE-1;i++){
    if(coordinates[i].x==coordinates[i+1].x){
      verticalityDegree++;
    }
  }
  return verticalityDegree;
}//RandomWalk::computeVerticalityDegree

void RandomWalk::moveLocally(int movingPoint){
  int nPossibilities=0, choice=0;
  Coordinates possibleCoordinates[3];
  if(movingPoint==1){
    if(lattice[doMod(coordinates[1].x+1,SIZE)]
              [doMod(coordinates[1].y,SIZE)]==0){
      possibleCoordinates[nPossibilities].x=coordinates[1].x+1;
      possibleCoordinates[nPossibilities].y=coordinates[1].y;
      nPossibilities++;
    }
    if(lattice[doMod(coordinates[1].x-1,SIZE)]
              [doMod(coordinates[1].y,SIZE)]==0){
      possibleCoordinates[nPossibilities].x=coordinates[1].x-1;
      possibleCoordinates[nPossibilities].y=coordinates[1].y;
      nPossibilities++;
    }
    if(lattice[doMod(coordinates[1].x,SIZE)]
              [doMod(coordinates[1].y+1,SIZE)]==0){
      possibleCoordinates[nPossibilities].x=coordinates[1].x;
      possibleCoordinates[nPossibilities].y=coordinates[1].y+1;
      nPossibilities++;
    }
    if(lattice[doMod(coordinates[1].x,SIZE)]
              [doMod(coordinates[1].y-1,SIZE)]==0){
      possibleCoordinates[nPossibilities].x=coordinates[1].x;
      possibleCoordinates[nPossibilities].y=coordinates[1].y-1;
      nPossibilities++;
    }
  }
  else if(movingPoint==SIZE){
    if(lattice[doMod(coordinates[SIZE-2].x+1,SIZE)]
              [doMod(coordinates[SIZE-2].y,SIZE)]==0){
      possibleCoordinates[nPossibilities].x=coordinates[SIZE-2].x+1;
      possibleCoordinates[nPossibilities].y=coordinates[SIZE-2].y;
      nPossibilities++;
    }
    if(lattice[doMod(coordinates[SIZE-2].x-1,SIZE)]
              [doMod(coordinates[SIZE-2].y,SIZE)]==0){
      possibleCoordinates[nPossibilities].x=coordinates[SIZE-2].x-1;
      possibleCoordinates[nPossibilities].y=coordinates[SIZE-2].y;
      nPossibilities++;
    }
    if(lattice[doMod(coordinates[SIZE-2].x,SIZE)]
              [doMod(coordinates[SIZE-2].y+1,SIZE)]==0){
      possibleCoordinates[nPossibilities].x=coordinates[SIZE-2].x;
      possibleCoordinates[nPossibilities].y=coordinates[SIZE-2].y+1;
      nPossibilities++;
    }
    if(lattice[doMod(coordinates[SIZE-2].x,SIZE)]
              [doMod(coordinates[SIZE-2].y-1,SIZE)]==0){
      possibleCoordinates[nPossibilities].x=coordinates[SIZE-2].x;
      possibleCoordinates[nPossibilities].y=coordinates[SIZE-2].y-1;
      nPossibilities++;
    }
  }
  else if((coordinates[movingPoint].x!=coordinates[movingPoint-2].x) and 
         (coordinates[movingPoint].y!=coordinates[movingPoint-2].y)){
    if((coordinates[movingPoint].y==coordinates[movingPoint-1].y) and
      (lattice[doMod(coordinates[movingPoint].x,SIZE)]
              [doMod(coordinates[movingPoint-2].y,SIZE)]==0)){
        possibleCoordinates[nPossibilities].x=coordinates[movingPoint].x;
        possibleCoordinates[nPossibilities].y=coordinates[movingPoint-2].y;
        nPossibilities++;
    }
    else if((coordinates[movingPoint].x==coordinates[movingPoint-1].x) and
           (lattice[doMod(coordinates[movingPoint-2].x,SIZE)]
                   [doMod(coordinates[movingPoint].y,SIZE)]==0)){
        possibleCoordinates[nPossibilities].x=coordinates[movingPoint-2].x;
        possibleCoordinates[nPossibilities].y=coordinates[movingPoint].y;
        nPossibilities++;
    }
  }
  
  if(nPossibilities>0){
    if(nPossibilities>1){
      choice=rand()%nPossibilities;
    }
    lattice[doMod(coordinates[movingPoint-1].x,SIZE)]
           [doMod(coordinates[movingPoint-1].y,SIZE)]=0;
    lattice[doMod(possibleCoordinates[choice].x,SIZE)]
           [doMod(possibleCoordinates[choice].y,SIZE)]=movingPoint;
    coordinates[movingPoint-1]=possibleCoordinates[choice];
  }

}//RandomWalk::moveLocally

void RandomWalk::slither(){
  int direction, cosAngle=0, sinAngle=0;
  direction=(rand()%3)+1;
  Coordinates possibleCoordinates;
  int translatedXCoordinate=0;
  int translatedYCoordinate=0;
  
  switch(direction){
    case 1:
      cosAngle=0;
      sinAngle=1;
      break;
    case 2:
      cosAngle=-1;
      sinAngle=0;
      break;
    case 3:
      cosAngle=0;
      sinAngle=-1;
  }
  
  if(isOneHead){
    translatedXCoordinate=coordinates[1].x-coordinates[0].x;
    translatedYCoordinate=coordinates[1].y-coordinates[0].y;
    possibleCoordinates.x=translatedXCoordinate*cosAngle-
                          translatedYCoordinate*sinAngle+coordinates[0].x;
    possibleCoordinates.y=translatedXCoordinate*sinAngle+
                          translatedYCoordinate*cosAngle+coordinates[0].y;
    if(lattice[doMod(possibleCoordinates.x,SIZE)]
              [doMod(possibleCoordinates.y,SIZE)]==0 or
      lattice[doMod(possibleCoordinates.x,SIZE)]
             [doMod(possibleCoordinates.y,SIZE)]==SIZE){
      lattice[doMod(coordinates[SIZE-1].x,SIZE)]
             [doMod(coordinates[SIZE-1].y,SIZE)]=0;
      for(int i=SIZE;i>1;i--){
	lattice[doMod(coordinates[i-2].x,SIZE)]
	       [doMod(coordinates[i-2].y,SIZE)]=i;
	coordinates[i-1]=coordinates[i-2];
      }
      lattice[doMod(possibleCoordinates.x,SIZE)]
             [doMod(possibleCoordinates.y,SIZE)]=1;
      coordinates[0]=possibleCoordinates;
    }
    else{
      isOneHead=false;
    }
  }
  else{
    translatedXCoordinate=coordinates[SIZE-2].x-coordinates[SIZE-1].x;
    translatedYCoordinate=coordinates[SIZE-2].y-coordinates[SIZE-1].y;
    possibleCoordinates.x=translatedXCoordinate*cosAngle-
                          translatedYCoordinate*sinAngle+coordinates[SIZE-1].x;
    possibleCoordinates.y=translatedXCoordinate*sinAngle+
                          translatedYCoordinate*cosAngle+coordinates[SIZE-1].y;
    if(lattice[doMod(possibleCoordinates.x,SIZE)]
              [doMod(possibleCoordinates.y,SIZE)]==0 or
      lattice[doMod(possibleCoordinates.x,SIZE)]
             [doMod(possibleCoordinates.y,SIZE)]==1){
      lattice[doMod(coordinates[0].x,SIZE)]
             [doMod(coordinates[0].y,SIZE)]=0;
      for(int i=1;i<SIZE;i++){
	lattice[doMod(coordinates[i].x,SIZE)][doMod(coordinates[i].y,SIZE)]=i;
	coordinates[i-1]=coordinates[i];
      }
      lattice[doMod(possibleCoordinates.x,SIZE)]
             [doMod(possibleCoordinates.y,SIZE)]=SIZE;
      coordinates[SIZE-1]=possibleCoordinates;
    }
    else{
      isOneHead=true;
    }
  }  
}//RandomWalk::slither

void RandomWalk::pivot(int pivotPoint, int pivotAngle){
  int cosAngle=0;
  int sinAngle=0;
  Coordinates translatedCoordinates;
  Coordinates possibleCoordinates;
  bool isSelfAvoiding=true;
  int i=pivotPoint;
  
  switch(pivotAngle){
    case 1:
      cosAngle=0;
      sinAngle=1;
      break;
    case 2:
      cosAngle=-1;
      sinAngle=0;
      break;
    case 3:
      cosAngle=0;
      sinAngle=-1;
  }
  while(i<SIZE and isSelfAvoiding){
    translatedCoordinates.x=coordinates[i].x-coordinates[pivotPoint-1].x;
    translatedCoordinates.y=coordinates[i].y-coordinates[pivotPoint-1].y;
    possibleCoordinates.x=translatedCoordinates.x*cosAngle-
                          translatedCoordinates.y*sinAngle+
                          coordinates[pivotPoint-1].x;
    possibleCoordinates.y=translatedCoordinates.x*sinAngle+
                          translatedCoordinates.y*cosAngle+
                          coordinates[pivotPoint-1].y;
    if(lattice[doMod(possibleCoordinates.x,SIZE)]
              [doMod(possibleCoordinates.y,SIZE)]!=0 and
      lattice[doMod(possibleCoordinates.x,SIZE)]
             [doMod(possibleCoordinates.y,SIZE)]<=i+1){
      isSelfAvoiding=false;
    }
    else{
      if(lattice[doMod(coordinates[i].x,SIZE)]
	        [doMod(coordinates[i].y,SIZE)]==i+1){
	lattice[doMod(coordinates[i].x,SIZE)][doMod(coordinates[i].y,SIZE)]=0;
      }
      lattice[doMod(possibleCoordinates.x,SIZE)]
             [doMod(possibleCoordinates.y,SIZE)]=i+1;
      coordinates[i].x=possibleCoordinates.x;
      coordinates[i].y=possibleCoordinates.y;
      i++;
    }
  }
  if(!isSelfAvoiding){
    switch(pivotAngle){
      case 1:
        cosAngle=0;
        sinAngle=-1;
        break;
      case 2:
        cosAngle=-1;
        sinAngle=0;
        break;
      case 3:
        cosAngle=0;
        sinAngle=1;
    }
    for(int j=pivotPoint;j<i;j++){
      lattice[doMod(coordinates[j].x,SIZE)][doMod(coordinates[j].y,SIZE)]=0;
      translatedCoordinates.x=coordinates[j].x-coordinates[pivotPoint-1].x;
      translatedCoordinates.y=coordinates[j].y-coordinates[pivotPoint-1].y;
      coordinates[j].x=translatedCoordinates.x*cosAngle-
                       translatedCoordinates.y*sinAngle+
                       coordinates[pivotPoint-1].x;
      coordinates[j].y=translatedCoordinates.x*sinAngle+
                       translatedCoordinates.y*cosAngle+
                       coordinates[pivotPoint-1].y;
    }
    for(int j=pivotPoint;j<SIZE;j++){
      lattice[doMod(coordinates[j].x,SIZE)][doMod(coordinates[j].y,SIZE)]=j+1;
    }
  }
}//RandomWalk::pivot

void RandomWalk::doTimestep(){
  switch(ALGORITHM_CHOICE){
    case 1:
      int movingPoint;
      movingPoint=(rand()%SIZE)+1;
      moveLocally(movingPoint);
      break;
    case 2:
      slither();
      break;
    case 3:
      int pivotPoint, pivotAngle;
      pivotPoint=(rand()%(SIZE-1))+1;
      pivotAngle=(rand()%3)+1;
      pivot(pivotPoint, pivotAngle);
  }
}//randomWalk::doTimestep

void RandomWalk::thermalize(){
  int deltaX=coordinates[SIZE-1].x-coordinates[0].x;
  int verticalityDegree=0;
  bool isDeltaXThermalized=false;
  bool isVerticalityDegreeThermalized=false;
  
  switch(WITH_THERMALIZATION_OUTPUT){
    case 0:{
      int thermalizationTime=0;
      while(!isDeltaXThermalized or !isVerticalityDegreeThermalized){
	thermalizationTime++;
        for(int i=0;i<MEASURE_FREQUENCY;i++){
          doTimestep();
        }
        deltaX=coordinates[SIZE-1].x-coordinates[0].x;
	verticalityDegree=computeVerticalityDegree();
	if(deltaX<0){
	  isDeltaXThermalized=true;
	}
	if(verticalityDegree==(SIZE-1)/2){
	  isVerticalityDegreeThermalized=true;
	}
      }
      for(int i=0;i<thermalizationTime;i++){
        for(int j=0;j<MEASURE_FREQUENCY;j++){
          doTimestep();
        }
      }
      break;
    }
    
    case 1:{
      ofstream output(thermalizationFilename.c_str());
      output <<"Iteration, End_to_end_length, Delta_X\n";
      for(int j=1;j<=SAMPLE_SIZE;j++){
        while(deltaX>=0){
          output <<j<<", "<<computeEndToEndLength()<<", "<<deltaX<<"\n";
          for(int i=0;i<MEASURE_FREQUENCY;i++){
            doTimestep();
          }
          deltaX=coordinates[SIZE-1].x-coordinates[0].x;
        }
        *this=RandomWalk();
	deltaX=coordinates[SIZE-1].x-coordinates[0].x;
      }
      output.close();
      break;
    }
    
    case 2:{
      ofstream output(thermalizationFilename.c_str());
      output <<"Iteration, End_to_end_length, Verticality_Degree\n";
      for(int j=1;j<=SAMPLE_SIZE;j++){
        while(verticalityDegree!=(SIZE-1)/2){
          output <<j<<", "<<computeEndToEndLength()<<", "
	         <<verticalityDegree<<"\n";
          for(int i=0;i<MEASURE_FREQUENCY;i++){
            doTimestep();
          }
          verticalityDegree=computeVerticalityDegree();
        }
        *this=RandomWalk();
        verticalityDegree=computeVerticalityDegree();
      }
    output.close();
    }
  }  
}//RandomWalk::thermalize

int main(){
  srand(time(NULL));
  RandomWalk randomWalk;
  
  randomWalk.thermalize();
  
  if(!WITH_THERMALIZATION_OUTPUT){
    ofstream output(dataFilename.c_str());
    output <<"Size, Algorithm, End_to_end_length, Gyration_radius\n";
    for(int i=0;i<SAMPLE_SIZE;i++){
      output <<SIZE<<", "<<ALGORITHM_CHOICE<<", "
             <<randomWalk.computeEndToEndLength()<<", "
	     <<randomWalk.computeGyrationRadius()<<"\n";
      for(int j=0;j<SIZE;j++){
        randomWalk.doTimestep();
      }
    }
    output.close();
  }
    
  return EXIT_SUCCESS;
}//main