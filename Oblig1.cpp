# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
#include <string>
#include "time.h"
#include <string.h>
using namespace std;
ofstream ofile;


float f(float x) {
    return 100*exp(-10*x);
}


float exact(float x) {
    return 1-(1-exp(-10))*x-exp(-10*x);
}


void vectorToFile(char* name, float in[], int n) {
    // Open file and write results to file:
    ofile.open(name);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i=0;i<n;i++) {
       ofile<< setw(15) << setprecision(8) << in[i];
       ofile<<"\n";
    }
    ofile.close();
}


void ErrorToFile(char* name, float in[], int n, float h) {
    // Open file and write results to file:
    ofile.open(name);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i=0;i<n;i++) {
       ofile<< setw(15) << setprecision(8) << log10(abs((in[i]-exact(i*h))/exact(i*h)));
       ofile<<"\n";
    }
    ofile.close();
}



void exactToFile(char* name, float h, int rows) {
    // Open file and write results to file:
    ofile.open(name);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i=0;i<rows;i++) {
       ofile<< setw(15) << setprecision(8) << exact(i*h);
       ofile<<"\n";
    }
    ofile.close();
}


//compare the entries of b with the expected entries from exact
float err(int rows, float h, float b[]) {
    float error=0;
    float temp;
    for(int i=1; i<rows; i++) {
        error+=log10(abs((b[i]-exact(i*h))/exact(i*h)));
        //if (temp>error) {error=temp;}
    }
    return error;
}


float worst(int rows, float h, float b[]) {
    float error=log10(abs(b[1]-exact(h)/exact(h)));
    float temp;
    for(int i=1; i<rows; i++) {
        temp=log10(abs(b[i]-exact(i*h)/exact(i*h)));
        if (temp>error) {error=temp;}
    }
    return error;
}


void makeThomasSystem(int rows, float h, float lower[], float middle[], float upper[], float b[]) {
    for(int i=0; i<rows-1; i++) {
        lower[i]=-1;
        middle[i]=2;
        upper[i]=-1;
        b[i]=f(i*h)*h*h;
    }
    middle[rows-1]=2;
    b[rows-1]=f((rows-1)*h)*h*h;
}


void makeThomasSystem(int rows, float h, float b[]) {
    for(int i=0; i<rows-1; i++) {
        b[i]=f(i*h)*h*h;
    }
    b[rows-1]=f((rows-1)*h)*h*h;
}


/*void printTrimat(float lower[], float middle[], float upper[], int rows) {
    for(int i=0; i<rows; i++) {
        for(int j=0; j<rows; j++) {
            if (j==i-1) { cout<<lower[i-1]<<"  "; }
            else if (j==i) { cout<<middle[i]<<"  "; }
            else if (j==i+1) { cout<<upper[i]<<"  "; }
            else { cout<<0; }
        }
     cout<<"\n";
    }
}


void printVec(float b[], int rows) {
    for(int i=0; i<rows; i++) {
        cout<<b[i]<<"\n";
    }
}*/


void triSolve(float lower[], float middle[], float upper[], int rows, float b[]) {
    //solve the system given by the tridiagonal matrix mat, and vector b
    //forward substitution
    for (int i=0; i<rows-1; i++) {
            //add an appropriate amount (coeff) of row i to row i+1
            float coeff=lower[i]/middle[i];
            //update the middle entries, the lower entries are 0, and the upper unchanged.
            //then update b
            middle[i+1]-=upper[i]*coeff;
            b[i+1]-=b[i]*coeff;
    }
    //back substitution
    b[rows-1]/=middle[rows-1];
    for(int i=rows-2; i>=0; i--) {
        b[i]=(b[i]-upper[i]*b[i+1])/middle[i];
    }
}


void thomasSolve(float b[], int rows) {
    //solve the system given by the tridiagonal matrix mat, and vector b
    //the outline for this procedure is in the pdf
    //Forward substitution
    //update the b[i], each i is incremented by one (in relation to the formulas in the pdf)
    //to accomodate for C++s way of handling the 1st element in an array)
    for(int i=1; i<=rows-1; i++) {
        b[i]+=i/((float) i+1)*b[i-1];
    }

    //Back substitution
    b[rows-1]*=(rows)/((float) rows+1); //(rows-1)/rows might be the correct one,
    //but this way this method agrees with the general one on v[rows-1]
    //solve iteratively
    for(int i=rows-2; i>=0; i--) {
        b[i]=i/((float) i+1)*(b[i]+b[i+1]);
    }
}


bool exercise1b() {
    cout<<"EXERCISE 1B\n";
    int rows=1;
    float h;
    float * lower;    float * middle;    float * upper;    float * b;
    clock_t start, finish;
    char* name;
    char ints[10];
    for (int n=1; n<=6; n++) {
        //initialize the new system to solve
        rows=rows*10;
        h=1/((float)rows);
        lower=new float[rows-1]; middle=new float[rows]; upper=new float[rows-1]; b=new float[rows];
        makeThomasSystem(rows,h,lower,middle,upper,b);
        //solve and time
        start=clock();
        triSolve(lower,middle,upper,rows,b);
        finish=clock();
        //print some results
        cout<<"n="<<rows<<"\n";
        cout<<"time:"<<(finish-start)/CLOCKS_PER_SEC<<"\n";
        cout<<"err:"<<err(rows,h,b)<<"\n";
        cout<<"worst:"<<worst(rows,h,b)<<"\n";
        //print results a file
        name="=n-1b";
        sprintf(ints, "%i",rows);
        strcat(ints,name);
        vectorToFile(ints,b,rows);
    }
    delete [] lower; delete [] middle; delete [] upper; delete [] b;
    return true;
}


void exercise1c() {
    cout<<"EXERCISE 1C\n";
    float h;
    float * b;
    clock_t start, finish;
    char* name;
    char ints[10];
    int rows=1000000;
    h=1/((float)rows);
    b= new float[rows];
    makeThomasSystem(rows,h,b);
    //print exact to a file
    name="=n-exact";
    sprintf(ints, "%i",rows);
    strcat(ints,name);
    exactToFile(ints,h, rows);

    //solve and time
    start=clock();
    thomasSolve(b,rows);
    finish=clock();

    //print some results
    cout<<"n="<<rows<<"\n";
    cout<<"time:"<<setprecision(10)<<(finish-start)/((float)CLOCKS_PER_SEC)<<"\n";
    cout<<"Relative error:"<<err(rows,h,b)<<"\n";
    cout<<"worst:"<<worst(rows,h,b)<<"\n";
    delete [] b;
}


void exercise1d() {
    //run the thomassolve algorithm for n=10^7, and print the relative error to a table(file)
    cout<<"EXERCISE 1D\n";
    float h;
    float * b;
    int rows=10000000;
    h=1/((float)rows);
    b= new float[rows];
    makeThomasSystem(rows,h,b);
    //solve
    thomasSolve(b,rows);
    //print results a file
    ErrorToFile("1D-table of errors.txt",b,rows,h);
    delete[] b;
}

int main() {
    exercise1b();
    exercise1c();
    exercise1d();
    return 0;
}
