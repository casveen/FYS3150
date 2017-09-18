# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
#include <string>
#include "time.h"
#include <string.h>
#include <armadillo>
using namespace std;
using namespace arma;
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


float worst(int rows, float h, float b[]) {
    float error=log10(abs(b[1]-exact(h)/exact(h)));
    float temp;
    for(int i=1; i<rows; i++) {
        temp=log10(abs(b[i]-exact(i*h)/exact(i*h)));
        if (temp>error) {error=temp;}
    }
    return error;
}


//overloaded makeThomassystem functions to make thomas systems for specific inputs
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

void makeThomasSystem(mat A, mat b,float h) {
    int rows=size(A)[1];
    for(int i=0; i<rows-1; i++) {
        A(i,i+1)=-1;
        A(i,i)=2;
        A(i+1,i)=-1;
        b(i)=f(i*h)*h*h;
    }
    A(rows-1, rows-1)=2;
    b(rows-1)=f((rows-1)*h)*h*h;
}


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
    int rows=1;
    float h;
    float * lower;    float * middle;    float * upper;    float * b;
    clock_t start, finish;
    char* name;
    char ints[10];
    //solve for n=10, 100 and 1000 and make tables
    for (int n=1; n<=3; n++) {
        //initialize the new system to solve
        rows=rows*10;
        h=1/((float)rows);
        lower=new float[rows-1]; middle=new float[rows]; upper=new float[rows-1]; b=new float[rows];
        makeThomasSystem(rows,h,lower,middle,upper,b);
        //solve and time
        start=clock();
        triSolve(lower,middle,upper,rows,b);
        finish=clock();
        //print results a file
        name="=n-1b";
        sprintf(ints, "%i",rows);
        strcat(ints,name);
        vectorToFile(ints,b,rows);
    }
    //free memory
    delete [] lower; delete [] middle; delete [] upper; delete [] b;
    return true;
}


void exercise1c(int rows) {
    cout<<"EXERCISE 1C("<<"n="<<rows<<")\n";
    float h;
    float * lower;    float * middle;    float * upper;    float * b;
    clock_t start, finish;
    h=1/((float)rows);
    b= new float[rows];
    //make and solve the system
    makeThomasSystem(rows,h,b);
    //print exact to a file
    /*name="=n-exact";
    sprintf(ints, "%i",rows);
    strcat(ints,name);
    exactToFile(ints,h, rows);
    */
    //solve and time
    start=clock();
    thomasSolve(b,rows);
    finish=clock();

    //print some results
    cout<<"time:"<<setprecision(10)<<(finish-start)/((float)CLOCKS_PER_SEC)<<"\n";
    //make the corresponding system and solve with the general method
    lower=new float[rows-1]; middle=new float[rows]; upper=new float[rows-1]; b=new float[rows];
    makeThomasSystem(rows,h,lower,middle,upper,b);
    start=clock();
    triSolve(lower,middle,upper,rows,b);
    finish=clock();
    cout<<"General algorithm:\n";
    cout<<"time:"<<setprecision(10)<<(finish-start)/((float)CLOCKS_PER_SEC)<<"\n";
    delete [] lower; delete [] middle; delete [] upper; delete [] b;
}


void exercise1d() {
    //run the thomassolve algorithm for n, and print the relative error to a table(file)
    cout<<"\nEXERCISE 1D\n";
    ofile.open("1D-table of errors.txt");
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    float h;
    float * b;
    int rows=10;
    ofile<<"n:          worst relative error:\n";
    for(int i=0; i<7; i++) {
        h=1/((float)rows);
        b= new float[rows];
        makeThomasSystem(rows,h,b);

        thomasSolve(b,rows);
        ofile<<setw(12)<<left<<rows;
        ofile<<setprecision(8)<<setw(12)<<worst(rows,h,b);
        ofile<<"\n";
        rows*=10;
    }
    ofile.close();
    delete[] b;
}


void exercise1e() {
    //compare to armadillos LUsolve, unfinished
    mat * A;
    mat * b;
    mat * L;
    mat * U;
    mat * P;
    int rows=10;
    float h;
    //solve for rows =10,100,1000
    for(int i=0; i<3; i++) {
        cout<<"\n";
        //make the system
        h=1/((float)rows);
        A=new mat(rows,rows);
        b=new mat(rows,1);
        L=new mat(rows,rows); U=new mat(rows,rows); P= new mat(rows,rows);
        makeThomasSystem(*A,*b,h);
        //Using lu or solve results in main giving absolutely no output...
        rows*=10;
    }
}


int main() {
    exercise1b();
    exercise1c(1000000);
    exercise1d();
    exercise1e();
    return 0;
}
