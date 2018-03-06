 /*
 FORSYTHE (C++ port) 
 Original by
 Grigoriev D. I.
 gr 2081/4
 350-22-94

 Adaptation by
 Rump Sergey aka SR3u (I just corrected some compilation errors and added somthing useful)
 gr 2084-1
 
 contacts:
 forsythe.cpp@gmail.com
 */
/*
 */
#pragma once
#ifndef _FORSYTHE_H_
#define _FORSYTHE_H_ 


#define FORSYTHE_STABLE 0 /*Unstable*/
#define FORSYTHE_VERSION "1.0.0.2.7a2"
#define FORSYTHE_SPLINE_STABLE 1
#define FORSYTHE_DECOMP_STABLE 1
#define FORSYTHE_LAGRANGE_STABLE 1
#define FORSYTHE_NEWTON_STABLE 1
#define FORSYTHE_QUANC8_RESULT_STABLE 1
#define FORSYTHE_MATRIX_STABLE 0  /*Unstable*/
#define FORSYTHE_VECTOR_STABLE 0  /*Unstable*/
#define FORSYTHE_LINEARSYSTEM_STABLE 1
#define FORSYTHE_TABLEFILE_STABLE 1

#define FORSYTHE_MAXITER 1000000000//Maximum iterations for some methods
#define FORSYTHE_MAXITER_2 100//Maximum iterations for some methods
#define FORSYTHE_USE_MAXITER 1

#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <iostream>
#include <string>
//#define REAL double//replaced by typedef
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define EXP(x) exp(x)
#define SQRT(x) sqrt(x)
#define SIN(x) sin(x)
#define COS(x) cos(x)
#define TAN(x) tan(x)
#define ATAN(x) atan(x)

#define _FORSYTHE_DEBUG_ 1

typedef double REAL;
typedef int INTEGER;

using namespace std;

//enum bool{false, true}; //remove if using MS VC or if your compiler has already type bool

inline REAL absval(REAL x)
{
    return fabs(x); 
}
inline INTEGER absval(INTEGER x)
{
	return abs(x);
}
#define ABS(x) absval(x)
#define SIGN(a,b) (b<0) ? (-absval(a)) : absval(a)

REAL ZEROIN(REAL AX,REAL BX,REAL (F)(REAL X),REAL TOL);
void INIT_URAND(INTEGER I);
REAL URAND();
void SVD(INTEGER M,INTEGER N,REAL **A,REAL *W,bool MATU,REAL **U,bool MATV,REAL **V,INTEGER &IERR,REAL *RV1);
void QUANC8(REAL (FUN)(REAL X),REAL A,REAL B,REAL ABSERR,REAL RELERR,REAL &RESULT,REAL &ERREST,INTEGER &NOFUN,REAL &FLAG);
REAL FMIN(REAL AX,REAL BX,REAL (F)(REAL X),REAL TOL);
void RKF45(void(F)(REAL T,REAL*Y,REAL*YP),INTEGER NEQN,REAL *Y,REAL &T,REAL TOUT,REAL &RELERR,REAL &ABSERR,REAL *WORK,INTEGER &IFLAG);


struct QUANC8_Result // Results of QUANC8 subroutine
{
    REAL value;
    REAL errest;
    INTEGER no_fun;
    REAL flag;
    static bool IsStable(void){return FORSYTHE_QUANC8_RESULT_STABLE;}
};

class SPLINE
{
    REAL *X,
	     *Y,
	     *B,
	     *C,
	     *D;
    INTEGER N;
public:
	SPLINE(INTEGER _N,REAL *_X,REAL *_Y);
    ~SPLINE();
    REAL Eval(REAL U);
    static bool IsStable(void){return FORSYTHE_SPLINE_STABLE;}
};
class Vector// vector
{
private:
    friend class Matrix;
    friend class LinearSystem;
    friend class DECOMP;
    REAL *Data;
    INTEGER size;
public:
    Vector(void);
    Vector(INTEGER Size);
    Vector(const Vector &V);
    Vector(INTEGER Size,REAL *D);
    ~Vector(void);
    INTEGER Size(void){return size;}
    void Set(INTEGER Size,REAL *M);//Fills the vector with array M
    REAL Norm(void);// norm of the vector
    REAL& operator[](INTEGER i);//acces to element #i
    void Print(void);// Output on screen
    void Print(const char *FileName, bool overwrite);// Output in FileName file, overwrites the file if overwrite is true
    void Print(FILE *OUT);// Output in file OUT
    void Print(const char *FileName,bool overwrite,const char *accuracy);//same as above, but you can adjust the number of digits after the poINTEGER (parameter accuracy)
    void Print(FILE *OUT,const char *accuracy);// number of digits should be in quotes for example, if you want 20 digits after the poINTEGER accuracy shall be "20"
    void Print(const char *accuracy);//
    bool operator==(const Vector &A);// compare operators
    bool operator!=(const Vector &A);
    bool operator==(REAL *A);// compare operators !insecure! (size of array M is assumed equal to the size of vector)
    bool operator!=(REAL *A);//
    Vector operator=(const Vector &M);// assignment
    Vector operator=(REAL *M);// assignment !insecure! (size of array M is assumed equal to the size of vector)
    Vector operator+(const Vector &M);// addition
    Vector operator-(const Vector &M);// substraction
    Vector operator*(REAL num);// multiplication
    Vector operator/(REAL num);// division
    REAL operator*(const Vector &vec);//scalar multiplication
    static bool IsStable(void)//returns true if Vector class is stable
    {return FORSYTHE_VECTOR_STABLE;}
};
class Matrix// square matrix
{
    friend class DECOMP;
    friend class Vector;
    friend class LinearSystem;
    friend void tred2(Matrix &M,Vector &D,Vector &E);
private:
    string name;
    REAL **Data;// data array
    INTEGER size;// matrix size
    static INTEGER _CountSpaces(REAL num,REAL maxnum);
    REAL& _GetEl(INTEGER i,INTEGER j);// i,j element access
public:
    Matrix(void);
    Matrix(INTEGER Size);
    Matrix(const Matrix &M);
    Matrix(INTEGER Size, REAL **D);//creates matrix filled with elements of the array D
    Matrix(INTEGER Size, REAL num);//creates matrix filled with number num
    //constructors below are almost the same as above, but the matrix will be created with name
    //matrixName - name of the matrix
    Matrix(const char* matrixName);
    Matrix(INTEGER Size,const char *matrixName);
    Matrix(const Matrix &M,const char* matrixName);
    Matrix(INTEGER Size, REAL **D,const char *matrixName);
    Matrix(INTEGER Size, REAL num,const char *matrixName);
    Matrix(string matrixName);
    Matrix(INTEGER Size,string matrixName);
    Matrix(const Matrix &M,string matrixName);
    Matrix(INTEGER Size, REAL **D,string matrixName);
    Matrix(INTEGER Size, REAL num,string matrixName);
    ~Matrix(void);
    void SetDiag(INTEGER num);//fills the matrix as the diagonal with num
    void SetIdentity(void);//fills the matrix as the Identity matrix
    void Set(INTEGER Size,REAL **M);//Fills the matrix with array M
    void Set(INTEGER Size,REAL **M,const char *Name);//Fills the matrix with array M and set its name to Name
    void Set(INTEGER Size,REAL **M,string Name);
    string GetName(void);//returns matrix name in stl string format
    string GetNameString(void);//returns matrix name in stl string format
    const char* GetName_c(void);//returns matrix name in c string format
    void SetName(const char *newName);//changes the matrix name to newName
    void SetName(string newName);//changes the matrix name to newName
    void Print(void);// Output on screen
    void Print(const char *FileName, bool overwrite);// Output in FileName file, overwrites the file if overwrite is true
    void Print(FILE *OUT);// Output in file OUT
    void Print(const char *FileName,bool overwrite,const char *accuracy);//Same as above, but you can adjust the number of digits after the poINTEGER (parameter accuracy)
    void Print(FILE *OUT,const char *accuracy);// number of digits should be in quotes for example, if you want 20 digits after the poINTEGER accuracy shall be "20"
    void Print(const char *accuracy);// 
    void PrintReadable(const char *FileName);// Output in FileName file in format for Read methods
    void PrintReadable(FILE *OUT);// Output in file OUT in format for Read methods
    void Read(const char *FileName);// read matrix from file IN
    void Read(FILE *IN);// read matrix from file IN
    /* file for methods Read must be like this:
    "<MatrixName>"
    <MatrixSize>
    <Elem> <Elem> <Elem> ... <Elem>
    <Elem> <Elem> <Elem> ... <Elem>
    ...............................
    <Elem> <Elem> <Elem> ... <Elem>
    
    (MatrixName is optional)
    Examples
     
    "SampleMatrixName"
    3
    100 9 30
    5 8 0
    33 1024 1 
    
    or like this
    
    3
    100 9 30
    5 8 0
    33 1024 1
    */
    static void Read(const char *FileName,Matrix **M,INTEGER MatrixNum);// read MatrixNum matrix from FileName file
    static void Read(FILE *IN,Matrix *M[],INTEGER MatrixNum);// read MatrixNum matrix from file IN
    /* File IN must be
     <MatrixNum>
     "<MatrixName>"
     <MatrixSize>
     <Elem> <Elem> <Elem> ... <Elem>
     <Elem> <Elem> <Elem> ... <Elem>
     ...............................
     <Elem> <Elem> <Elem> ... <Elem>
     "<MatrixName>"
     <MatrixSize>
     <Elem> <Elem> <Elem> ... <Elem>
     <Elem> <Elem> <Elem> ... <Elem>
     ...............................
     <Elem> <Elem> <Elem> ... <Elem>
     ...............................
     MatrixNum is optional
     */
    REAL* GetRow(INTEGER i);//Returns the row number i
    Vector GetRowVec(INTEGER i);//Same as above, but returns a vector
    REAL* GetColumn(INTEGER j);//Returns the column number j
    Vector GetColumnVec(INTEGER j);//Same as above, but returns a vector
    REAL MaxEl(void);// Max element of the matrix
    REAL MaxRow(INTEGER i);// Max element in the row number i
    REAL MaxCol(INTEGER j);// Max element in the column number j
    REAL AbsMaxRow(INTEGER i);// Absolute value of max element in the row number i
    REAL AbsMaxCol(INTEGER j);// Absolute value of max element in the column number j
    INTEGER Size(void){return size;}// Size of the Matrix
    REAL& operator()(INTEGER i,INTEGER j);// i,j element access
    bool operator==(const Matrix &A);// compare operators
    bool operator!=(const Matrix &A);
    bool operator==(REAL **A);// compare operators !insecure! (size of array M is assumed equal to the size of matrix)
    bool operator!=(REAL **A);//
    Vector operator*(Vector M);//multiply by vector
    Matrix operator=(const Matrix &M);// assignment
    Matrix operator=(REAL **M);//assignment !insecure! (size of array M is assumed equal to the size of matrix)
    Matrix operator+(const Matrix &M);// additionw
    Matrix operator-(const Matrix &M);// substraction
    Matrix operator*(const Matrix &M);// multiplication
    Matrix operator/(const Matrix &M);// division
    Matrix operator!(void);// reverse matrix
    Matrix Transpose(void);// transpose matrix
    Matrix T(void);// transpose matrix
    REAL D(void);// determinant
    REAL Det(void);// determinant
    REAL Determinant(void);// determinant
    REAL Norm(void);// norm of the matrix (see void Norm_III(void))
    REAL Norm_III(void);// norm of the matrix (canonical number III sqrt(sum[i=0;i<size;i++](sum[j=0;j<size;j++]((abs(a[i][j]))^2))
    REAL Cond(void);// conditional
    void AddColumn(INTEGER ColNum,REAL *C);
    void AddRow(INTEGER RowNum,REAL *C);
    void AddColumn(INTEGER ColNum,Vector &V);
    void AddRow(INTEGER RowNum,Vector &V);
    void EigenvaluesAndEigenvectors(Matrix &EigenVectorsInCol,Vector &EigenValues);//Finds EigenValues and EigenVectors of the matrix
                                                                                   //(Only in REAL numbers!)
                                                                                   //EigenVectors are stored in rows of matrix EigenVectorsInCol
                                                                                   //Eigenvalues are in vector EigenValues
    static bool IsStable(void);// returns true if this class is stable
};

class DECOMP
{
    INTEGER     N;
    REAL    **A;
    INTEGER     *IPVT;
    REAL    *WORK;
    REAL    COND;
    public:
	    DECOMP(INTEGER _N,REAL **_A);
        DECOMP(Matrix &M);// see Matrix class
	    ~DECOMP();
    void decomp(INTEGER _N,REAL **_A);
    REAL    Cond(){return COND;};
    REAL    Det();
    void    Solve(REAL *B);
    void    Solve(Vector *B){Solve(B->Data);}
    static bool IsStable(void){return FORSYTHE_DECOMP_STABLE;}
};

class LinearSystem// linear system like A*X=B where A is a Matrix, X and B are vectors 
{//A*X=B
private:
    friend class DECOMP;
    friend class Vector;
    friend class Matrix;
    Matrix _A;
    REAL  *_B;
public:
    LinearSystem(void);
    LinearSystem(const Matrix &A,REAL *B);
    LinearSystem(REAL **A,REAL *B,INTEGER Size);
    LinearSystem(REAL **A,const Vector &B,INTEGER Size);
    LinearSystem(const Matrix &A,const Vector &B);
    ~LinearSystem(void);
    void SetMatrix(const Matrix &M);//Set the matrix of the system
    void SetMatrix(REAL **M,INTEGER Size);
    void SetVector(REAL *vec);//Set the vetor of the system
    void SetVector(const Vector &vec);
    void SetAll(const Matrix &M,const Vector &vec);//Set the matrix and the vector of the linear system
    void SetAll(const Matrix &M,REAL *vec);
    void Solve(REAL *X);//Solve the system (result in X)
    void Solve(Vector *X);
    void Solve(Vector &X);
    void Solve(REAL **X);//Backward compatibility
    void LCGM(Vector &X);//Linear Conjugate Gradient Method
};



void QUANC8(REAL (FUN)(REAL X),REAL A,REAL B,REAL ABSERR,REAL RELERR,QUANC8_Result &RESULT);

//tests
void TestAll();
void TestZeroin(void);
void TestUrand(void);
void TestSVD(void);
void TestSpline(void);
void TestDecomp(void);
void TestQuaNC8(void);
void TestFMin(void);
void TestRKF45(void);

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
//SR3u
void TestLagrange(void);
void TestLagrange(REAL (FUN)(REAL X));
void TestLagrange(REAL (FUN)(REAL X),REAL X);//Works Better if X is near [0,1] interval or in it
void TestNewton(void);

class Lagrange
{
public:
    // Lagrange interpolation polynomial by PointsNum points using table of function values //
    static REAL Polynomial(INTEGER PointsNum,REAL *TableF,REAL *TableX,REAL X);//for more info see !ReadMe.rtf
    // Lagrange interpolation polynomial by PointsNum points using function's formula
    static REAL Polynomial(REAL (FUN)(REAL X),INTEGER PointsNum,REAL *TableX,REAL X);//for more info see !ReadMe.rtf
    //Automatic Lagrange interpolation polynomial by PointsNum points using function's formula on [A,B]
    static REAL AutoPolynomial(REAL (FUN)(REAL X),INTEGER PointsNum,REAL A,REAL B,REAL X);//for more info see !ReadMe.rtf
    //It would be great if you will test this function using void TestLagrange(REAL (FUN)(REAL X)) or  function void TestLagrange(REAL (FUN)(REAL X),REAL X)
    static bool IsStable(void);//returns true if Lagrange class is stable
private:
    static REAL _Polynomial(INTEGER PointsNum,REAL *TableF,REAL *TableX,REAL X);//Don't use this
};


class Newton
{
public:
    // Newton interpolation polynomial by PointsNum points using table of function values //
    static REAL Polynomial(INTEGER PointsNum,REAL *TableF,REAL *TableX,REAL X);//for more info see !ReadMe.rtf
    // Newton interpolation polynomial by PointsNum points using function's formula
    static REAL Polynomial(REAL (FUN)(REAL X),INTEGER PointsNum,REAL *TableX,REAL X);//for more info see !ReadMe.rtf
    //Automatic Newton interpolation polynomial by PointsNum points using function's formula on [A,B]
    static REAL AutoPolynomial(REAL (FUN)(REAL X),INTEGER PointsNum,REAL A,REAL B,REAL X);//for more info see !ReadMe.rtf
    static REAL CountDivDif(INTEGER start, INTEGER end,REAL *TableF,REAL *TableX);//divided difference
    static REAL CountDivDif(INTEGER start, INTEGER end,REAL (FUN)(REAL X),REAL *TableX,INTEGER PointsNum);//divided difference using function formula
    static bool IsStable(void);//returns true if Newton class is stable
};

class TableFile
{//Table file for Import Table option in AdvGrapher
private:
    string FName;
public:
    TableFile(void);
    TableFile(const char *FileName);
    TableFile(string FileName);
    void SetFile(const char *FileName);//Sets the file to FileName
    void SetFile(string FileName);
    void Print(INTEGER x,INTEGER F);//Prints a row into the table file
    void Print(INTEGER x,REAL F);
    void Print(REAL x,INTEGER F);
    void Print(REAL x,REAL F);
    void Print(INTEGER x,REAL F,const char *accuracy);//Prints a row into the table file with accuracy digits after the point
    void Print(REAL x,INTEGER F,const char *accuracy);
    void Print(REAL x,REAL F,const char *accuracy);
    void PrintE(INTEGER x,REAL F);//Prints a row in eponential representation into the table file
    void PrintE(REAL x,INTEGER F);
    void PrintE(REAL x,REAL F);
    void PrintE(INTEGER x,REAL F,const char *accuracy);//Prints a row in eponential representation into the table file with accuracy digits after the point
    void PrintE(REAL x,INTEGER F,const char *accuracy);
    void PrintE(REAL x,REAL F,const char *accuracy);
    void PrintA(INTEGER x,REAL F);//Prints a row in shortest (exponential or not) representation into the table file
    void PrintA(REAL x,INTEGER F);
    void PrintA(REAL x,REAL F);
    void PrintA(INTEGER x,REAL F,const char *accuracy);//Prints a row in shortest (exponential or not) representation into
    void PrintA(REAL x,REAL F,const char *accuracy);//the table file with accuracy digits after the point
    void PrintA(REAL x,INTEGER F,const char *accuracy);
    static bool IsStable(void);// returns true if this class is stable
};
namespace FORSYTHE_CONFIG
{
    extern bool SHOWERRORS;
    extern bool SHOWWARNINGS;
    extern bool ERRORREACT;
    extern bool WARNINGSREACT;
}
class FORSYTHE
{
public:
    static const char* Version(void){return FORSYTHE_VERSION;}
    static bool IsStable(){return FORSYTHE_STABLE;}
    static void PrintVersion(void){cout<<"\nFORSYTHE_CPP Version "<<FORSYTHE_VERSION<<"\n";}
    static void _PrintInfo(void)
    {
        static bool RunInfo=true;
        if (RunInfo)
        {
            RunInfo=false;
            FORSYTHE::PrintInfo();
        }
        
    }
    static void PrintInfo(void)
    {
        cout<<"\n/*\n * FORSYTHE (C++ port)\n";
        cout<<" * Version "<<FORSYTHE_VERSION<<" ";
        if(FORSYTHE::IsStable())
            cout<<"(Stable)";
        else
            cout<<"(Unstable)";
        cout<<"\n * \n * Original: by Grigoriev D. I. gr 2081/4 tel 350-22-94";
        cout<<"\n * Adaptation by: Rump Sergey aka SR3u gr 2084-1";
        cout<<"\n * Testing by:\n *   -- Lobov Konstantin aka Anther gr 2084-1\n *   -- Rump Sergey aka SR3u gr 2084-1\n *   -- Borisevich Leonid aka DaZzBouRGh gr 2084-1";
        cout<<"\n * \n * Thank you for using our product!\n * If you have some questions or propositions contact us forsythe.cpp@gmail.com\n * If you've found a bug please send us an e-mail with theme\n * \"FORSYTHE_CPP BUG\" and describe it.\n";
        cout<<" */\n\n";
    }
    static void ErrorsOff(void)//turn off any error messages
    {
        _PrintInfo();
        FORSYTHE_CONFIG::SHOWERRORS=false;
    }
    static void ErrorsOn(void)//turn on any error messages
    {
        _PrintInfo();
        FORSYTHE_CONFIG::SHOWERRORS=true;
    }
    static void ErrorsReactOff(void)//turn off waiting for user reaction to any error messages
    {
        _PrintInfo();
        FORSYTHE_CONFIG::ERRORREACT=false;
    }
    static void ErrorsReactOn(void)//turn on waiting for user reaction to any error messages
    {
        _PrintInfo();
        FORSYTHE_CONFIG::ERRORREACT=true;
    }
    static void WarningsOff(void)//turn off any error messages with danger level <= 100 (they're just warnings)
    {
        _PrintInfo();
        FORSYTHE_CONFIG::SHOWWARNINGS=false;
    }
    static void WarningsOn(void)//turn on any error messages with danger level <= 100 (they're just warnings)
    {
        _PrintInfo();
        FORSYTHE_CONFIG::SHOWWARNINGS=true;
    }
    static void WarningsReactOff(void)//turn off waiting for user reaction to any error messages with danger level <= 100 (they're just warnings)
    {
        _PrintInfo();
        FORSYTHE_CONFIG::WARNINGSREACT=false;
    }
    static void WarningsReactOn(void)//turn on waiting for user reaction to any error messages with danger level <= 100 (they're just warnings)
    {
        FORSYTHE_CONFIG::WARNINGSREACT=true;
    }
    static void Error(string msg,unsigned char danger)//show error message msg and ask if you wish to continue
    //great for warn the user about errors
    //msg - message to show, danger - shows how dangerous this error is (0-255)
    {
        Error(msg.c_str(),danger);
    }
    static void Error(const char *msg,unsigned char danger)//show error message msg and ask if you wish to continue
    //great for warn the user about errors
    //msg - message to show, danger - shows how dangerous this error is (0-255)
    {
        FORSYTHE::_PrintInfo();
        if(danger<25)
            if(!FORSYTHE_CONFIG::SHOWWARNINGS)
                return;
        if(danger>=25)
            if(!FORSYTHE_CONFIG::SHOWERRORS)
                return;
        char choice=0;
        while(choice!='C'||choice!='A'||choice!='c'||choice!='a')
        {
            if(danger>25)
            {cout<<"\nERROR: "<<msg<<"\nDanger level: "<<(int)danger<<"/255\n";}
            else{cout<<"\nWARNING: "<<msg<<"\nDanger level: "<<(int)danger<<"/255\n";}
            if(danger == 255)
            {
                cout<<"OOPS!!!!\n";
                cout<<"Critical Error!!!\n";
                if(FORSYTHE_CONFIG::ERRORREACT)
                    getchar();
                exit(-1024);
            }
            if(danger>=250)
                cout<<"\nIt's strongly recommended for you to stop here!\n";
            if(danger>=200 && danger<250)
                cout<<"That's not good, please stop\n";
            if(danger>100&&danger<200)
                cout<<"\nThis may hurt!\n";
            if(danger>25&&danger<=100)
            {
                cout<<"\nNothing serious, but don't say I didn't warn you!\n";
                if(!FORSYTHE_CONFIG::WARNINGSREACT)
                    return;
            }
            if(danger<=25)
            {
                cout<<"\nIt's not an error\nI just wanted you to know.\n";
                if(!FORSYTHE_CONFIG::WARNINGSREACT)
                return;
            }
            if(!FORSYTHE_CONFIG::ERRORREACT)
                return;
            cout<<"\npress C to continue A to abort: ";
            cin>>choice;
            if(choice=='A'||choice=='a')
            {
                cout<<"Clever choice! Now fix this!!\n";
                exit(-1);
            }
            if(choice=='C'||choice=='c')
            {
                cout<<"\nI don't want to continue, but if you insist...\nAnyway, I've warned you!!\n\n";
                return;
            }
        }
    }
};
Vector LinearConjugateGradientMethod(Matrix A, Vector b);
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// Backward compatibility
REAL LagrangePolynomial(INTEGER PointsNum,REAL *TableF,REAL *TableX,REAL X);//See Lagrange class
REAL LagrangePolynomial(REAL (FUN)(REAL X),INTEGER PointsNum,REAL *TableX,REAL X);//See Lagrange
REAL AutoLagrangePolynomial(REAL (FUN)(REAL X),INTEGER PointsNum,REAL A,REAL B,REAL X);//See Lagrange class
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
typedef FORSYTHE FORSYTHE_CPP;//Backward compatibility
#define FORSYTHE_CPP_STABLE FORSYTHE_STABLE
#define FORSYTHE_CPP_VERSION FORSYTHE_VERSION
#define FORSYTHE_CPP_SPLINE_STABLE FORSYTHE_SPLINE_STABLE
#define FORSYTHE_CPP_DECOMP_STABLE FORSYTHE_DECOMP_STABLE
#define FORSYTHE_CPP_LAGRANGE_STABLE FORSYTHE_LAGRANGE_STABLE
#define FORSYTHE_CPP_QUANC8_RESULT_STABLE FORSYTHE_QUANC8_RESULT_STABLE
#define FORSYTHE_CPP_MATRIX_STABLE FORSYTHE_MATRIX_STABLE
#define FORSYTHE_CPP_NEWTON_STABLE FORSYTHE_NEWTON_STABLE
#define FORSYTHE_CPP_VECTOR_STABLE FORSYTHE_VECTOR_STABLE
#define FORSYTHE_CPP_LINEARSYSTEM_STABLE FORSYTHE_LINEARSYSTEM_STABLE
#define FORSYTHE_CPP_TABLEFILE_STABLE FORSYTHE_TABLEFILE_STABLE
#define _FORSYTHE_CPP_DEBUG_ _FORSYTHE_DEBUG_
#define FORSYTHE_DEBUG _FORSYTHE_DEBUG_

#endif

