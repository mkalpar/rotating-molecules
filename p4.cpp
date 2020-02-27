// MERVE KALPAR
// 20200106
// this project creates cell, atom, and molecule classes, 
// along with many member functions,
// they read molecules from a xyz file, assign them to arrays, 
// finds their total mass, density, volume
// prints to the screen
// rotates them around axises

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <numeric>
#include <math.h>
using namespace std;

//function for cross product, used in calculating the volume
template <class atype>
vector<atype> crossProduct(vector<atype> vect1, vector<atype> vect2, vector<atype> crossP){
	crossP.push_back(vect1[1]*vect2[2]-vect1[2]*vect2[1]);
	crossP.push_back(vect1[0]*vect2[2]-vect1[2]*vect2[0]);
	crossP.push_back(vect1[0]*vect2[1]-vect1[1]*vect2[0]);
	return crossP;
}

//function for dot product, used in calculating the volume
template <class atype>
double dotProduct(vector<atype> vect1, vector<atype> vect2){
	
	double sum = 0;	
	for (int i = 0; i < 3; ++i)
	{
		sum += vect1[i]*vect2[i];
	}
	return sum;
}

class Cell {
public:
	double* a1;
	double* a2;
	double* a3;
private:
	double volume;
public:
	Cell() : a1(0),a2(0),a3(0) {};
	Cell(string FileName) {
		ifstream inputFile(FileName);
		a1 = new double[3];
		a2 = new double[3];
		a3 = new double[3];
		inputFile >> a1[0] >> a1[1] >> a1[2] >> 
					a2[0]>> a2[1]>> a2[2] >> a3[0]>> 
					a3[1]>> a3[2];
	};

	~Cell(){};

	void printCell(){
		cout <<"a1:" <<*(a1)<< " " << *(a1+1) << " " << *(a1+2)<<endl;
		cout <<"a2:" <<*(a2)<< " " << *(a2+1) << " " << *(a2+2)<<endl;
		cout <<"a3:" <<*(a3)<< " " << *(a3+1) << " " << *(a3+2)<<endl;
	};

	double calculateVolume(){
		vector<double> v1;
		vector<double> v2;
		vector<double> v3;
		for (int i = 0; i < 3; ++i) 
		{	
			v1.push_back(a1[i]);
			v2.push_back(a2[i]);
			v3.push_back(a3[i]);
		}

		vector<double> cpV;
		vector<double> crossPVector = crossProduct(v2,v3,cpV);
		volume = dotProduct(v1,crossPVector);
		return volume;
	};

	double getVolume(){
		return volume;
	};
};

class Atom {
public:
	string e;
	double x,y,z;

	Atom() : x(0),y(0),z(0),e("_") {};
	Atom(string e, double x, double y, double z) : e(e),x(x),y(y),z(z) {};
	~Atom() {};
	void printAtom(){
		cout <<x<<y<<z<<e<<endl;
	};
	void rotX(double radVal){
		x = x;
		double tmpY = y;
		y = tmpY*cos(radVal)-z*sin(radVal);
		z = tmpY*sin(radVal)+z*cos(radVal);
	};
	void rotY(double radVal){
		double tmpX = x;
		x = tmpX*cos(radVal) + z*sin(radVal);
		y = y;
		z = -tmpX*sin(radVal) + z*cos(radVal);
	};
	void rotZ(double radVal){
		double tmpX = x;
		x = x*cos(radVal) - y*sin(radVal);
		y = tmpX*sin(radVal) + y*cos(radVal);
		z = z;
	};
	
};

class Molecule {
public:
	int nAtom;
	Atom* atoms;
	string comment;
	double center;
	Cell* cell;
	double volume, totalMass, density;
	double CenterX,CenterY,CenterZ;
	double sumOfX,sumOfY,sumOfZ;
	Molecule(string FileName){
		ifstream inputFile(FileName);
		inputFile >> nAtom >> comment;
		cout << nAtom << " "<<comment;
		atoms = new Atom[nAtom];
		int i = 0;	
		while(!inputFile.eof()){
			inputFile >> atoms[i].e >> atoms[i].x >> atoms[i].y >> atoms[i].z; 
			//cout <<"e: " <<atoms[i].e <<" x: "<< atoms[i].x << " y: "<<atoms[i].y <<" z: "<< atoms[i].z <<endl;
			i++;
		}
	};
	~Molecule(){delete [] atoms;}; 

	void assignCell(string FileName){
		cell = new Cell(FileName);
	};
	void printMol(){
		cout <<"Printing molecule!";
		cout << endl << nAtom <<endl;
		cout <<comment<<endl;
		for (int i = 0; i < nAtom; ++i)
		{	
			cout << (*(atoms+i)).e << " "<<(*(atoms+i)).x 
			<< " " << (*(atoms+i)).y << " " <<
			(*(atoms+i)).z <<endl;	
		}
	};
	void writeMol(string FileName){
		ofstream outputFile(FileName,ios::out | ios::app);

		if(!outputFile) {throw "Output file not found!"; abort();} 

		outputFile << nAtom <<endl<< comment<<endl;
		int x = 0;	
		while(x<nAtom){
			outputFile<< atoms[x].e << " "<<atoms[x].x << " "<< atoms[x].y <<" " << atoms[x].z<<endl; 
			x++;
		}
	};
	void findCenter(){
		//sums the X,Y,Z coordinates each, and divides the sums by the number of atoms
		cout << "find center: "<<endl; 
		sumOfX= 0.0;sumOfY = 0.0;sumOfZ= 0.0; 

		for (int i = 0; i < nAtom; ++i)
		{	
			sumOfX+= (*(atoms+i)).x; 
			sumOfY+= (*(atoms+i)).y; 
			sumOfZ+= (*(atoms+i)).z;
		}
		cout << "sum of X coords:"<<sumOfX<<endl;
		cout << "sum of Y coords:"<<sumOfY<<endl;
		cout << "sum of Z coords:"<<sumOfZ<<endl;
		CenterX = sumOfX/nAtom;
		CenterY = sumOfY/nAtom;
		CenterZ = sumOfZ/nAtom;
		cout << "Center X: "<<CenterX<<endl;
		cout << "Center Y: "<<CenterY<<endl;
		cout << "Center Z: "<<CenterZ<<endl;
	};
	void center2origin(){
		// moves the center to origin
		// this method doesnt work correctly yet
		// it is supposed to be center of gravity?
		double moveX = CenterX/nAtom;
		double moveY = CenterY/nAtom;
		double moveZ = CenterZ/nAtom;
		cout << "move x's by : " << moveX<<endl;
		cout << "move y's by : " << moveY<<endl;
		cout << "move z's by : " << moveZ<<endl;

		for (int i = 0; i < nAtom; ++i)
		{	
			if((*(atoms+i)).x<0) {(*(atoms+i)).x= (*(atoms+i)).x+moveX;}
			else if ((*(atoms+i)).x>0) {(*(atoms+i)).x= (*(atoms+i)).x-moveX;}
			if((*(atoms+i)).y<0) {(*(atoms+i)).y= (*(atoms+i)).y+moveY; }
			else if ((*(atoms+i)).y>0) {(*(atoms+i)).y= (*(atoms+i)).y-moveY;}
			if((*(atoms+i)).z<0) {(*(atoms+i)).z= (*(atoms+i)).z+moveZ; }
			else if ((*(atoms+i)).z>0) {(*(atoms+i)).z= (*(atoms+i)).z-moveZ;}
			//(*(atoms+i)).y= (*(atoms+i)).y+moveY; 
			//(*(atoms+i)).z= (*(atoms+i)).z+moveZ;
		}
	};
	void rotX(double radVal){
		for (int i = 0; i < nAtom; ++i) atoms[i].rotX(radVal); 
	};
	void rotY(double radVal){
		for (int i = 0; i <  nAtom; ++i) atoms[i].rotY(radVal); 
	};
	void rotZ(double radVal){
		for (int i = 0; i < nAtom; ++i) atoms[i].rotZ(radVal); 
	};
	double calcTotMass(){
		//there are four distinct atoms in the file, 
		//calculates the total mass by adding up the mass for all these repeating atoms
		totalMass = 0.0;
		for (int i = 0; i < nAtom; ++i)
		{	
			//cout << (*(atoms+i)).e<< " ";
			if((*(atoms+i)).e == "C") {totalMass += 12.011;}
			else if((*(atoms+i)).e == "H") {totalMass += 1.008;}
			else if((*(atoms+i)).e == "N") {totalMass += 14.007;}
			else if((*(atoms+i)).e == "O") {totalMass += 15.999;}
			
		}
		return totalMass;
	};
	void printCell() {
		//redirects to the member function of the Cell class
		cell->printCell();
	};
	double getVolume(){
		//redirects to the member functions of the Cell class
		volume = 0.0;
		cell->calculateVolume();
		volume = cell->getVolume();
		return volume;
	};
	void calcDensity(){
		//calls other member functions of the same class
		calcTotMass();
		getVolume();
		density = totalMass/volume;

	};
	double getDensity(){
		return density;
	};

};

int main() {

	Cell c1;
	Cell c2("cell_file");

	cout << "volume cell:"<<c2.calculateVolume()<<endl;
	cout << "volume cell:"<<c2.getVolume()<<endl;
	
	Atom a1("test",1,2,3);
	cout << "atom: "<< a1.x <<" "<< a1.y<< " " << a1.z <<endl;
	Atom a1Rot;
	a1.rotX(1.00);
	cout << "atom rotated around X: "<< a1.x <<" "<< a1.y<< " " << a1.z <<endl;
	a1.rotY(0.778);
	cout << "atom rotated around Y: "<< a1.x <<" "<< a1.y<< " " << a1.z <<endl;


	Molecule m1("16fuCyt.xyz");
	m1.assignCell("cell_file");
	cout << endl <<"Volume of cells: " << m1.getVolume() << endl;
	cout << "Total mass of atoms: " << m1.calcTotMass() << endl;
	m1.calcDensity();
	cout << "Density: " << m1.getDensity() << endl;
	m1.printCell();
	//m1.printMol();
	m1.findCenter();
	//m1.center2origin();// this doesnt properly work 
	//m1.findCenter();
	 
	 for (int i = 0; i < 628; ++i)
	 {
	 	m1.rotZ(0.01);
	 	m1.writeMol("rotated.xyz");
	 }
	

	return 0;
}


