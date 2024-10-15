/*******************************************************************************
template <typename T> class vector
{
  
  public:  
    vector(int rows);
    virtual ~vector();
    void generate();
    void print();
    T operator[](int i);
    void append(T i);
    T* to_number() const;
    

  private:
    int _rows;   
    T* _vector;
    int _pos;
};


template <typename T> class matrix
{
  
  public:  
    matrix(int rows, int cols); 
    virtual ~matrix();
    void generate();
    void generate_tridiagonal_matrix_symmetric();
    void print();   
    double at(int i, int j);
    const T* const to_number()const;
    vector<T> *diagonalPrincipal();
    vector<T> *subdiagonalPrincipal();            
    T* to_number();
    
  private:
    int _rows;
    int _cols;
    T* _matrix;
};


class LAPACK
{
  public:  
    LAPACK(int rows)
    void dgemm();
    void dsyevd();
    void dstevd();
    
  private:
   int _rows;
   
};

********************************************************************************/

#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

extern "C" void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc ); 
extern "C" void dsyevd_( char *jobz, char *uplo, int *n, double *A, int *lda, double *w, double *work, int *ldwork ,int *lwork, int *liwork,int *info); 
extern "C" void dstevd_(char *jobz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *ldwork, double *iwork, int *liwork, int *info);

using namespace std;

/*vector*/
template <typename T> class vector
{
  
  public:  
    vector(int rows)
    {
      _pos = 0;
      _vector = new T[rows]; 
      _rows = rows;  
    }
    
    virtual ~vector()
    { 
     delete[] _vector; 
    }
    
    void generate()
    {
     for(int i = 0; i < _rows; ++i)
	     _vector[i] = rand()%(100-1)*1;
    
    }
   
    void print()
    {
    	for(int i = 0; i < _rows; ++i)
         cout << _vector[i] << '\t';
  
      cout << '\n';
    
    }
    
    T operator[](int i)
    {
      return _vector[i];
    }
 
      
    void append(T i)
    {
       _vector[_pos] = i;
       _pos += 1;     
    }
    
    T* to_number() const
    {
      return _vector;
    } 
   

  private:
    int _rows;
       
    T* _vector;
    int _pos;
};


/*matrix*/
template <typename T> class matrix
{
  
  public:  
    matrix(int rows, int cols)
    {
      _matrix = new T[rows*cols]; 
      _rows = rows;
      _cols  = cols;
    }
    
    virtual ~matrix()
    { 
     delete[] _matrix; 
    }
    
    void generate()
    {
     for(int i = 0; i < _rows; ++i)
	    for(int j = 0; j < _cols; ++j)
	     _matrix[i + j*_rows] = rand()%(100-1)*1;
    
    }/*generate*/
 
    void generate_tridiagonal_matrix_symmetric()
    {
      for(int i = 0; i < _rows; ++i)
	    { 
       for(int j = 0; j < _cols; ++j)
	     {   
                           
         if(i==j)
           _matrix[i + j*_rows] = rand()%(100-1)*1;
          
           else if(i==(j-1)) 
               _matrix[i + j*_rows] = rand()%(100-1)*1;
             
                else if((i-1)==j) 
                   _matrix[i + j*_rows] = _matrix[j + i*_rows];
          
          else
             _matrix[i + j*_rows] = 0;       
       }    
      }
    
    }/*generate_tridiagonal_matrix_symmetric*/
 
   
    void print()
    {
      for(int i = 0; i < _rows; ++i)
      {
	     for(int j = 0; j < _cols; ++j)
	       cout << _matrix[i+ j*_rows] << '\t';
    
      cout << '\n';
      }
    cout << '\n';
    }/*print*/
       
       
    double at(int i, int j)
    {
      return _matrix[i+ j*_rows];
    }/*at*/
    
    const T* const to_number()const
    {
      return _matrix;
    }/*to_number*/
    
    
    vector<T> *diagonalPrincipal()
    {
       vector<T> *v = new vector<T>(_rows);
    
       for(int i = 0; i < _rows; ++i)
       {
         v->append( at(i, i));
       }
       
       return v;
    }/*diagonalPrincipal*/  
    
    
    vector<T> *subdiagonalPrincipal()
    {
       vector<T> *s = new vector<T>(_rows-1);
    
       for(int i = 1; i < _rows; ++i)
         s->append( at(i, i-1));
       
       
       return s;
    }/*subdiagonalPrincipal*/
          
          
    T* to_number()
    {
     return _matrix;
    }/*to_number*/


  private:
    int _rows;
    int _cols;
    T* _matrix;
};



class LAPACK
{
  public:  
    LAPACK(int rows)
    {
     _rows = rows;
    }
      
    void dgemm()
    {
     int m = _rows;
     int n = _rows;
     int k = _rows;
     int lda = _rows; 
     int ldb = _rows; 
     int ldc = _rows; 
     double alpha = 1.;
     double beta =  0.;
     char transa ='N';
     char transb ='N';
      
      matrix<double> *A = new matrix<double>(_rows, _rows);
      matrix<double> *B = new matrix<double>(_rows, _rows);
      matrix<double> *C = new matrix<double>(_rows, _rows);
 
      A->generate();
      A->print(); 
   
      B->generate();
      B->print(); 
  
    /*--------------------------------------------------------------------------------------------------------------------------------*/
      dgemm_(&transa, &transb, &m, &n, &k, &alpha, A->to_number(), &lda, B->to_number(), &ldb, &beta, C->to_number(), &ldc);
    /*--------------------------------------------------------------------------------------------------------------------------------*/
   
      C->print(); 
   
      delete A;
      delete B;
      delete C;  
    
    }/*dgemm*/
  

    void dsyevd()
    {
      char JOB = 'N';          
      char UPLO = 'L';         
      int N = _rows;
      int LDA = _rows;
      int LWORK;
      int LIWORK; 
      int info;
 
      matrix<double> *A = new matrix<double>(_rows, _rows);
      vector<double> *w = new vector<double>(_rows);
  
      A->generate();
      A->print(); 
      
      LWORK = (2*N + 1)*10;     
      LIWORK = 1;
      
      vector<double> *WORK  = new vector<double>(LWORK);
      vector<int>   *IWORK  = new vector<int>(LIWORK);
  
    /*--------------------------------------------------------------------------------------------------------------------------------*/
         dsyevd_(&JOB, &UPLO, &N, A->to_number(), &LDA, w->to_number(), WORK->to_number(), &LWORK, IWORK->to_number(), &LIWORK, &info); 
    /*--------------------------------------------------------------------------------------------------------------------------------*/
  
      w->print();
    
      delete A;
      delete w;
      delete WORK;
      delete IWORK;    
    
    }/*dsyevd*/
  
  
    void dstevd()
    {
     char JOBZ = 'N';   
     int N = _rows;
     int LDZ = _rows;
     int info;
    
     int  LWORK = 1;
     int  LIWORK = 1; 
     
     vector<double> *Z     = new vector<double>(LDZ);
     vector<double> *WORK  = new vector<double>(LWORK);
     vector<double> *IWORK = new vector<double>(LIWORK);
      
     matrix<double> *A = new matrix<double>(N, N);
     vector<double> *D = NULL;
     vector<double> *S = NULL;
   
     A->generate_tridiagonal_matrix_symmetric();
     A->print(); 

     D = A->diagonalPrincipal();
     D->print();

     S = A->subdiagonalPrincipal();
     S->print();

    /*------------------------------------------------------------------------------------------------------------------------------------------*/  
        dstevd_(&JOBZ, &N, D->to_number(), S->to_number() , Z->to_number(), &LDZ, WORK->to_number(), &LWORK, IWORK->to_number(), &LIWORK, &info); 
    /*-------------------------------------------------------------------------------------------------------------------------------------------*/
  
     D->print();
     
     delete Z;
     delete WORK;
     delete IWORK; 
     delete A;
     delete D;
     delete S;
    
    }/*dstevd*/ 
    
  private:
   int _rows;
   
};
