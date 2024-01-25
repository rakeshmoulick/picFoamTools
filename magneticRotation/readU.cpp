/*
Objective: This code extracts the values of U from various run files of picFoam-magneticRotation Tutorial. 
The extracted data is written on a file "results.txt"
*/
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
using namespace std;

// Define a file-ID
FILE *file;

// Function readFile extracts and writes data to a new file "results.txt"
void readFile(istream &strm, double time, FILE *file)
{
    string line;
    int N;
    // The while statement points to the line and gets the line number 
    // where starting word "object" is encountered 
    while ( getline( strm, line ) && line.find( "object" ) == string::npos ) ;
    
    // each getline command points to a new line of the file. 
    getline(strm, line);
    getline(strm, line);
    getline(strm, line);
    getline(strm, line); 
    //cout << "getline: " << line << endl;

    // strinstream streams values to N
    stringstream (line) >> N;
    cout << "N = " << N << endl;
    
    // You need to see the architecture of the line. In this case it starts as 1((. 
    // These three characters are assigned to c1, c2 and c3 respectively. 
    char c1, c2, c3;
    
    // The next integer values are assigned to x, y and z respectively.
    double x, y, z;

    for (int i=0; i<N; i++)
    {        

        stringstream(line) >> c1 >> c2 >> c3 >> x >> y >> z;
        //cout << "x=" << x << " " << "y=" << y << " " << "z=" << z << " " << endl; 

        fprintf(file, "%g \t %g \t %g \t %g\n",time, x, y, z);
        
        // For writing data on screen 
        //printf("%g \t %g \t %g \t %g\n",time, x, y, z);

        fflush(file);
    }

}

int main()
{
   char NAME[50];
   sprintf(NAME,"results.txt");
   // open a file and clean it. 
   file = fopen(NAME,"w");
   // Re-open the file to just append (write at end of the last file line) 
   file = fopen(NAME,"a");

   double k=0; 
   // you need to know the number of folders created based on writing interval
   for (int i=0; i<=80; i++)
   {
        /*-----------------------------------*/
        char NAME[50];

        sprintf(NAME,"%g/lagrangian/pic/U",k);
        //cout << "Name:"  << NAME << endl;
        string filename = NAME;
        
        //string filename = "2e-09/lagrangian/pic/U";             
        /*-----------------------------------*/
        ifstream strm( filename );                        
        /*-----------------------------------*/
        double time; 
        stringstream (filename) >> time;
        //cout << "time :" << time << endl; 
        /*-----------------------------------*/
        readFile( strm, time, file);
        k += 0.5e-09;  
   }
   
}