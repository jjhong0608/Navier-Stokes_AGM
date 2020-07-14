#ifndef HEADER_H
#define HEADER_H

#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>
#include <limits>
#include <functional>
#include <unordered_map>

// pi value
#ifndef PI
#define PI M_PI
#endif

// zero value (0.0)
#ifndef ZeroValue
#define ZeroValue 0.0000000000000000E+00
#endif

// unit value (1.0)
#ifndef UnitValue
#define UnitValue 1.0000000000000000E+00
#endif

// marginal zero
#ifndef NearZero
#define NearZero 5.0000000000000000E-14
#endif

// 반올림 margin
#ifndef Margin
#define Margin 10
#endif

// dimensions
#ifndef D2
#define D2 2
#endif

#ifndef D3
#define D3 3
#endif

#ifndef SIGDIGITS
#define SIGDIGITS 8
#endif

#ifndef True
#define True true
#endif

#ifndef False
#define False false
#endif

using namespace std;

int SolverType;
/* Solver Type               */
/* Elliptic solver     : 0   */
/* Heat transfer solver: 1   */

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              Declaration                               |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Vector;
class CoLMaT;
class OPRT;
class Rot;
class Cell;
class Thing;
class Triangle;
class Line2D;
class Line3D;
class Section;
class Region;
class AxialLine;
class Thing_AxialLine;
class ControlData;
class AxialData;
class Coordinate;
class Point;

struct _xData;
struct _yData;

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              Vector Class                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Vector {
public:
  Vector ();
  Vector (int);
  Vector (double, double);
  Vector (double, double, double);
  Vector (const Vector &);
  ~Vector ();

  Vector & setvalue (int, double);
  Vector & setvector (int, double *);
  Vector & setvector (int, ...);
  Vector & setvector (Vector &);
  Vector & setzerovector ();
  Vector & setzerovector (int);
  Vector & setonevector ();
  Vector & setonevector (int);
  Vector & setunit ();
  Vector & setEk (int, int);
  Vector & cut (int);
  double   getvalue (int k) const;
  int      getdim () const {return dim;}
  Vector   getunit ();

  double   dot (Vector &);
  double   norm ();
  double   triple (Vector &, Vector &);
  Vector   cross (Vector &);
  Vector & add (const Vector &);
  Vector & scalar (double);

  Vector & showcontents (const char *);
  Vector & showcontents (const char *, const char *);
  Vector & print (const char *);
  Vector & print (FILE *, const char *);
  Vector & print (ofstream & , const char *);
  Vector & gnuplot (FILE *, const char *);

  double &  operator[] (int);
  Vector    operator+ (Vector &);
  Vector    operator- (Vector &);
  Vector    operator* (double);
  double    operator* (Vector &);

  Vector &  operator= (const Vector &); // assignment operator
  Vector &  operator+= (const Vector &); // assignment operator

private:
  int      dim;
  double  *comp;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                (CoLMat) Matrix of column vectors Class                 |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class  CoLMaT
{
  public :
  CoLMaT ();
  CoLMaT (int);
  CoLMaT (int, int);
  CoLMaT (const CoLMaT &);
  ~CoLMaT ();

  int       getcoldim () const {return col_dim;}
  int       getrowdim () const {return vec_n;}
  Vector *  getmatrix () const {return colvec;}
  Vector &  getcolvector (int) const;
  double    getvalue (int, int);
  CoLMaT &  setvalue (int, int, double);
  CoLMaT &  reset ();
  CoLMaT &  setInxn (int);
  CoLMaT &  setOnxn (int);
  CoLMaT &  cut (int);
  CoLMaT &  putvector (int, Vector &);
  CoLMaT &  add (CoLMaT &);
  CoLMaT &  scalar (double);
  CoLMaT &  rowmatrix (Vector &);
  CoLMaT &  colmatrix (Vector &);
  CoLMaT &  rankone (Vector &, Vector &);

  CoLMaT &  showcontents (const char *);
  CoLMaT &  showsize (const char *);

  Vector &  operator[] (int);
  CoLMaT    operator+ (CoLMaT &);
  CoLMaT    operator* (CoLMaT &);
  Vector    operator* (Vector &);

  CoLMaT &  operator= (const CoLMaT &);

  private :
  int     col_dim;
  int     vec_n;
  Vector *colvec;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|           (OPRT) Operators of projectors and Rotators Class            |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class OPRT
{
  public :
  OPRT ();
  OPRT (int);
  OPRT (const OPRT &);
  ~OPRT ();

  OPRT &    setvalue (int, int, double);
  double    getvalue (int, int) const;
  CoLMaT *  getmatrix () const {return mat;}
  int       getdim () const {return dim;}
  OPRT &    setmatrix (int);

  OPRT &    OP_onto (Vector &);
  OPRT &    OP_onto (Vector &,  Vector &);
  OPRT &    OP_onto (int, Vector *);
  OPRT &    COP_onto (Vector &);
  OPRT &    Omega (Vector &);

  Vector    acting (Vector &);
  OPRT &    add (const OPRT &);

  OPRT &    showcontents (const char *);

  Vector    operator* (Vector &);
  OPRT &    operator= (const OPRT &);

  private :
  int      dim;
  CoLMaT  *mat;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                     (ROT) Rotators (3D only) Class                      |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Rot
{
  public :
  Rot ();
  Rot (Vector &);
  Rot (const Rot &);
  ~Rot ();

  Vector    getRotAxis () const {return omega;}
  Rot &     setRotAxis (Vector &);
  Rot &     setOPRT (Vector &);
  Vector    acting (Vector &, double);
  OPRT      getPw () const {return Pw;}
  OPRT      getImPw () const {return ImPw;}
  OPRT      getOmg () const {return Omg;}
  OPRT      getOPRT (int) const;

  Rot &     showcontents (const char *);

  private :
  Vector    omega;
  OPRT      Pw, ImPw, Omg;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                 (Cell) Unit cell of linked list Class                  |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class   Cell
{
  public :
  Cell ();
  Cell (int);
  Cell (Vector &);
  Cell (Vector &, char, double, char);
  Cell (Vector &, char, double, char, Vector &);
  Cell (int, Vector &);
  Cell (const Cell &);
  ~Cell ();

  int         Index () const {return index;}
  Vector      Item () const {return vec;}
  double      Item (int offs) {return vec[offs];}
  Vector      Normal () const {return normal;}
  double      Normal (int offs) {return normal[offs];}
  Cell*       Next () const {return nextcell;}
  char        Cond () const {return cond;}
  char        Status () const {return status;}
  double      Value () const {return value;}
  AxialLine * Xaxialline () const {return xaxialline;}
  AxialLine * Yaxialline () const {return yaxialline;}

  Cell &   SetIndex (int off){index = off; return *this;}
  Cell &   SetItem (Vector &);
  Cell &   SetItem (int, Vector &);
  Cell &   SetNormal (Vector &);
  Cell &   SetNormal (int, Vector &);
  Cell &   SetNext (Cell*);
  Cell &   SetCond (char);
  Cell &   SetStatus (char);
  Cell &   SetValue (double);
  Cell &   SetAxialLine (AxialLine*, AxialLine*);
  Cell &   AppendCell (int);
  Cell &   AppendCell (Vector &);
  Cell &   AppendCell (Vector &, char, double, char);
  Cell &   AppendCell (Vector &, char, double, char, Vector &);
  Cell &   AppendCell (int, Vector &);
  Cell &   AppendCell (int, Vector &, char, double, char);
  Cell &   AppendCell (int, Vector &, char, double, char, Vector &);
  Cell &   AppendCell (int, Vector &, Cell*, Cell*);

  Cell &   showcontents (const char *);

  private :
  int        index;
  AxialLine *xaxialline;
  AxialLine *yaxialline;

  Vector     vec;
  Vector     normal;
  char       cond;
  char       status;
  double     value;
  Cell      *nextcell;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|               (Thing) Master of Cells (Linked List) Class               |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class   Thing
{
  public :
  Thing ();
  Thing (const Thing &);
  ~Thing ();

  int      Howmany () const {return ncell;}
  Cell  *  Head () const {return  headcell;}
  Cell  *  Run () const {return run;}
  Thing &  MakeCell (int);
  Thing &  MakeCell (Vector &);
  Thing &  MakeCell (int, Vector &);
  Thing &  MakeCell (Vector &, char, double, char);
  Thing &  MakeCell (int, Vector &, char, double, char);
  Thing &  MakeCell (Vector &, char, double, char, Vector &);
  Thing &  MakeCell (int, Vector &, char, double, char,Vector &);
  Thing &  RemoveCell ();
  Thing &  RemoveAllCell ();
  Thing &  TakeoffHeadCell ();
  Thing &  GoPost ();
  Thing &  GoPost (int);
  Thing &  GoPre ();
  Thing &  GoPre (int);
  Thing &  GoHead ();
  Thing &  GoTail ();
  Thing &  MoveRun (int);
  Thing &  MoveRun (Cell *);
  Thing &  SetHead (Cell *);
  Thing &  SetNcell (int nc) {ncell = nc; return *this;}

  int      CountCells ();
  int      Numbering ();
  int      Numbering (int);
  Cell *   WhichCell (int);
  Thing &  HeadChange ();
  Thing &  HeadChange (int);
  Thing &  Swapping ();
  Thing &  Swapping (Cell *);
  int      Ordering (char *, int);
  Thing &  Sorting (const char *);
  Thing &  Sorting (const char *, int);
  Thing &  Sorting (char *, char, int);

  Thing &  Import (Thing &, char);
  Thing &  Import (Thing &, int, char);

  int      OnRun ();
  Thing &  showcontents (const char *);

  // operators
  Thing &  operator= (Thing &);
  Cell  &  operator[] (int);

  private :
  int      ncell;
  Cell    *headcell;
  Cell    *run;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                          Triangle in 3D Class                          |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Triangle
{
  public :
  // Constructors and Destructors

  Triangle ();
  Triangle (Vector & , Vector &, Vector &);
  Triangle (Vector & , Vector &, Vector &, char, double);
  ~Triangle ();

  //  Access Functions

  double      GetArea (){ return Area; }
  Vector &    GetVertex (int offset){ return Vertex[offset]; }
  Vector &    GetNormal (){ return Normal; }
  OPRT &      GetOPRT (){ return CP_n; }
  char        GetBDcondition () {return cond;}
  double      GetBDvalue () {return value;}
  Triangle &  SetIndex (int offs) { index = offs; return *this; }
  Triangle &  SetArea ();
  Triangle &  SetNormal ();
  Triangle &  SetNormal (Vector &);
  Triangle &  SetTriangle (Vector &, Vector &, Vector &);
  Triangle &  SetTriangle (Vector & , Vector &, Vector &, char, double);
  Triangle &  SetBD (char, double);
  double      Distance (Vector &);
  Triangle &  PutNormal (Vector &);
  Triangle &  PutVertex (int, Vector &);
  Triangle &  SetOPRT (Vector &);
  Triangle &  SetOPRT ();
  int         Inside (Vector &);

  //  Methods

  Triangle & CalculateProperties ();
  Vector     MapFromReference (double, double);
  Triangle & showcontents (const char *);

  // Operators

  Vector &   operator[] (int);

  private :
  int     index;
  double  Area;
  Vector  Vertex[D3];
  Vector  Normal;
  OPRT    CP_n;
  double  value;
  char    cond;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              Line2D Class                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Line2D
{
  public :
  // Constructors and Destructors

  Line2D ();
  Line2D (Vector & , Vector &);
  Line2D (double, double, double, double);
  Line2D (Vector & , Vector &, char, double);
  Line2D (double, double, double, double, char, double);
  ~Line2D ();

  //  Access Functions

  Vector &   GetStart (){ return start; }
  Vector &   GetEnd (){ return end; }
  Vector &   GetNormal (){ return Normal; }
  Vector &   GetTangent (){ return Tangent; }
  char       GetBDcondition () {return cond;}
  double     GetBDvalue () {return value;}
  double     Distance (Vector &);
  double     GetLength () { return Len; }
  Line2D &   SetLength ();
  Line2D &   SetNormal ();
  Line2D &   SetNormal (Vector &);
  Line2D &   SetNormal (double, double);
  Line2D &   SetLine (Vector &, Vector &);
  Line2D &   SetLine (double, double, double, double);
  Line2D &   SetLine (Vector & , Vector &, char, double);
  Line2D &   SetLine (double, double, double, double, char, double);
  Line2D &   SetOPRT (Vector &);
  Line2D &   SetOPRT ();
  Line2D &   SetBD (char, double);
  OPRT &     GetOPRT () { return CP_q; }
  int        IsCross (Line2D &, Vector &, double &);
  // int        ScanCross (int, Line2D *, Vector** &);
  int        ScanCross (int, Line2D *, Thing &, double);

  int         Inside (Vector &);

  //  Methods

  Line2D &   CalculateProperties ();
  Vector     MapFromReference (double);
  Line2D &   showcontents (const char *);

  // Operators

  Vector &   operator[] (int);

  private :
  Vector  start, end;
  Vector  Normal, Tangent;
  OPRT    CP_q;
  double  Len;
  double  value;
  char    cond;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              Line3D Class                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Line3D
{
  public :
  // Constructors and Destructors

  Line3D ();
  Line3D (Vector & , Vector &);
  Line3D (double, double, double, double, double, double);
  ~Line3D ();

  //  Access Functions

  Vector &   GetStart (){ return start; }
  Vector &   GetEnd (){ return end; }
  Vector &   GetNormal (){ return Normal; }
  Vector &   GetTangent (){ return Tangent; }
  double     Distance (Vector &);
  double     GetLength () { return Len; }
  Line3D &   SetLength ();
  Line3D &   SetNormal ();
  Line3D &   SetLine (Vector &, Vector &);
  Line3D &   SetLine (double, double, double, double, double, double);
  Line3D &   SetOPRT (Vector &);
  Line3D &   SetOPRT ();
  OPRT &     GetOPRT () { return CP_q; }
  int        IsCross (Triangle &, Vector &, double &);
  // int        ScanCross (int, Triangle *, Vector **&);
  int        ScanCross (int, Triangle *, Thing &);

  int         Inside (Vector &);

  //  Methods

  Line3D &   CalculateProperties ();
  Vector     MapFromReference (double);
  Line3D &   showcontents (const char *);

  // Operators

  Vector &   operator[] (int);

  private :
  Vector     start, end;
  Vector     Normal, Tangent;
  OPRT       CP_q;
  double     Len;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                             Section Class                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Section {
private:
  string    name; // section name
  int       nE;   // the number of elemens
  Line2D   *SL2D; // 2D line element
  Triangle *SL3D; // 3D line element


public:
  Section ();
  Section (string);
  Section (string, int);
  Section (string, int, Line2D*);
  Section (string, int, Triangle*);
  virtual ~Section ();

  Section & setName (string nm) {name = nm; return *this;}
  string    getName () {return name;}

  Section & setEltnum (int n) {nE = n; return *this;}
  int       getEltnum () {return nE;}

  Section   & setElt (Line2D *ln2d) {SL2D = ln2d; return *this;}
  Section   & setElt (Triangle *tr3d) {SL3D = tr3d; return *this;}
  Line2D    * get2DElt (int n) {return &SL2D[n];}
  Triangle  * get3DElt (int n) {return &SL3D[n];}

  Section & setSection (string, int);
  Section & setSection (string, int, Line2D*);
  Section & setSection (string, int, Triangle*);

  Section & showcontents (const char*);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                             Region Class                               |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Region {
private:
  string    name; // region name
  int       nS;   // the number of sections
  Section **sctn; // sections
  double   *sign; // section sign
  Vector    gridinfo; // grid informataion
  Line2D   *xgrid2d; // 2D background linesegment along x-axis
  Line2D   *ygrid2d; // 2D background linesegment along y-axis
  Line3D   *xgrid3d; // 3D background linesegment along x-axis
  Line3D   *ygrid3d; // 3D background linesegment along y-axis
  Line3D   *zgrid3d; // 3D background linesegment along z-axis
  int       nXgrid; // the number of xgrid
  int       nYgrid; // the number of ygrid
  int       nZgrid; // the number of zgrid

public:
  Region ();
  Region (string);
  Region (string, int);
  Region (string, int, Vector*);
  virtual ~Region ();

  Region & setName (string nm) {name = nm; return *this;}
  string   getName () {return name;}

  Region & setGridinfo (Vector* src) {gridinfo = *src; return *this;}
  Vector * getGridinfo () {return &gridinfo;}

  Region & setSectionnum (int n) {nS = n; return *this;}
  int      getSectionnum () {return nS;}

  Region  &  setSection (Section*);
  Region  &  setSection (Section*, int);
  Section *  getSection (int n) {return sctn[n];}

  Region & setSign (string, int);
  double   getSign (int n) {return sign[n];}

  Region & make2Dgrid ();
  Region & make3Dgrid ();

  int getXgridnum () {return nXgrid;}
  int getYgridnum () {return nYgrid;}
  int getZgridnum () {return nZgrid;}

  Line2D * get2Dxgrid (int n) {return &xgrid2d[n];}
  Line2D * get2Dygrid (int n) {return &ygrid2d[n];}

  Line3D * get3Dxgrid (int n) {return &xgrid3d[n];}
  Line3D * get3Dygrid (int n) {return &ygrid3d[n];}
  Line3D * get3Dzgrid (int n) {return &zgrid3d[n];}

  Region & setRegion (string, int);
  Region & setRegion (string, int, Vector*);

  Region & showcontents (const char*);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                            Axialline Class                             |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class AxialLine {
private:
  Cell      *start_pt;
  Cell      *end_pt;
  int        index;
  AxialLine *nextaxialline;
  Vector     center_pt;

public:
  AxialLine ();
  AxialLine (Cell*, Cell*);
  AxialLine (int, Cell*, Cell*);
  virtual ~AxialLine ();

  int         Index  () const {return index;}
  Cell      * Start  () const {return start_pt;}
  Cell      * End    () const {return end_pt;}
  AxialLine * Next   () const {return nextaxialline;}
  Vector    * Center ()       {return &center_pt;}

  AxialLine & SetIndex (int);
  AxialLine & SetNext (AxialLine*);
  AxialLine & setpoints (Cell*, Cell*);

  AxialLine & AppendAxialLine (Cell*, Cell*);
  AxialLine & AppendAxialLine (int, Cell*, Cell*);

  AxialLine & CalcCenterPt ();

  AxialLine & showcontents (const char*);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                         Axialline Linked list                          |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Thing_AxialLine {
private:
  int        nline;
  AxialLine *headline;
  AxialLine *run;

public:
  Thing_AxialLine ();
  Thing_AxialLine (const Thing_AxialLine&);
  virtual ~Thing_AxialLine ();

  int Howmany () const {return nline;}
  AxialLine * Head () const {return headline;}
  AxialLine * Run () const {return run;}

  Thing_AxialLine & MakeLine (Cell*, Cell*);
  Thing_AxialLine & MakeLine (int, Cell*, Cell*);
  Thing_AxialLine & RemoveLine ();
  Thing_AxialLine & RemoveAllLine ();
  Thing_AxialLine & TakeoffHeadLine ();
  Thing_AxialLine & GoPost ();
  Thing_AxialLine & GoPost (int);
  Thing_AxialLine & GoPre ();
  Thing_AxialLine & GoPre (int);
  Thing_AxialLine & GoHead ();
  Thing_AxialLine & GoTail ();
  Thing_AxialLine & MoveRun (int);
  Thing_AxialLine & MoveRun (AxialLine*);
  Thing_AxialLine & SetHead (AxialLine*);
  Thing_AxialLine & SetNline (int nl) {nline = nl; return *this;}

  int         CountLines ();
  int         Numbering ();
  int         Numbering (int);
  AxialLine * WhichLine (int);

  Thing_AxialLine & Import (Thing_AxialLine&, char);
  Thing_AxialLine & Import (Thing_AxialLine&, int, char);

  Thing_AxialLine &  showcontents (const char *);

  // operators
  Thing_AxialLine & operator= (Thing_AxialLine&);
  AxialLine       & operator[] (int);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              ControlData                               |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class ControlData {
private:
  double   dt              ; //time interval
  double   initialTime     ; //initial time
  double   terminalTime    ; //terminal time
  double   t               ; //Present time
  int      timeStep        ; //time step
  int      rest_out_intvl  ; //out interval
  string   AxialFile       ; //axial data file name
  string   output_sol      ; //output solution file name
  string   output_del      ; //output derivative file name
  string   boundary        ;
  string   xaxial          ;
  string   yaxial          ;

public:
  ControlData ();
  virtual ~ControlData ();

  double Dt ()           const {return this->dt;}
  double InitialTime ()  const {return this->initialTime;}
  double TerminalTime () const {return this->terminalTime;}
  double T ()            const {return this->t;}

  int TimeStep ()       const {return this->timeStep;}
  int Rest_Out_Intvl () const {return this->rest_out_intvl;}

  string Axialfile ()  const {return this->AxialFile;};
  string Output_Sol () const {return this->output_sol;}
  string Output_Del () const {return this->output_del;}
  string Boundary () const {return this->boundary;}
  string Xaxial () const {return this->xaxial;}
  string Yaxial () const {return this->yaxial;}

  // control data 읽기
  ControlData & LoadCtrlData (string, string, double, double, int, double);
  // control data 출력
  ControlData & ShowCtrlData ();
};

ControlData::ControlData () {}
ControlData::~ControlData () {}

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                               AxialData                                |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class AxialData {
private:
  int      nx           ; // total number of x-axial line
  int      ny           ; // total number of y-axial line
  int      bd_pts_num   ; // total number of the boundary poins
  int      xaxial_num   ; // total number of x-axial lines
  int      yaxial_num   ; // total number of y-axial lines
  int      xxaxial_num  ; // size of xxaxial
  int      yyaxial_num  ; // size of yyaxial
  int      in_pts_num   ; // total number of the cross points
  int      phi_pts_num  ; // total number of the cross points except interface points
  int      pts_num      ; // total number of points

  int    **EWNS_index   ; // EWNS index of the point    [ 4 * n + i ] ( 6 * in_pts_num   )
  int     *xaxial_index ; // index of points along x-axial line
  int     *yaxial_index ; // index of points along y-axial line
  int     *xxaxial_index; // first number of each x-axial line
  int     *yyaxial_index; // first number of each y-axial line
  int     *ptsTOin_pts  ; // all points to inner points
  int     *in_ptsTOpts  ; // inner points to all points
  int     *ptsTOphi_pts ; // all points to inner points
  int     *phi_ptsTOpts ; // inner points to all points
  int    **axial_index  ; // index of the axial line      [ 3 * n + i ] ( 3 * nx * ny * nz )

  double **pts          ; // coordinate of the points     [ 2 * n + i ] ( 2 * pts_num      )
  double **xaxial       ; // x-axial lines                [ 4 * l + i ] ( 4 * nx           )
  double **yaxial       ; // y-axial lines                [ 4 * l + i ] ( 4 * ny           )
  double  *b_u          ; // bd value u      [      n + i ] (      1 * in_pts_num )
  double  *mp_u         ; // material property
  double **normal       ; // normal vector

  char    *bc_u         ; // bd condition u  [      n + i ] (      1 * in_pts_num )

public:
  AxialData ();
  virtual ~AxialData ();

  double Pts (int, char);
  double XYaxial (char, int, int);
  double Boundaryvalue (int);
  double MaterialProperty (int);
  double Normal (int, char);

  char Boundarycondition (int);

  int Pts_Num () const {return pts_num;}
  int In_Pts_Num () const {return in_pts_num;}
  int Phi_Pts_Num () const {return phi_pts_num;}
  int XYaxial_Num (char);
  int XXYYaxial_Num (char);
  int EWNS_Index (int, char);
  int XXYYaxial_Index (char, int);
  int XYaxial_Index (char, int);
  int PtsTOPts (char, int);
  int Axial_Index (int, char);

  AxialData & SetPtsTOpts (char, int, int);
  AxialData & AllocatePhipts (int);
  AxialData & SortEWNS ();

  AxialData & ExportAxialData (ControlData*);
  AxialData & LoadAxialData (string);
  AxialData & AssignBoundaryValue ();
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                             Coordinate                                 |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Coordinate {
private:
  int dim;
  double *comp;

public:
  Coordinate ();
  Coordinate (int);
  Coordinate (double, double);
  Coordinate (double, double ,double);
  Coordinate (const Coordinate &);
  virtual ~Coordinate ();

  Coordinate & SetCoordinate (char, double);
  Coordinate & SetCoordinate (double, double);
  Coordinate & SetCoordinate (double, double, double);
  Coordinate & SetCoordinate (Coordinate &);

  int Dimension () const {return this->dim;}

  double & Value (char) const;
  double Distance (const Coordinate &);

  double     & operator[] (const int);
  double     & operator[] (const int) const;
  Coordinate & operator=  (const Coordinate &);
  double       operator-  (const Coordinate &);
  bool         operator== (const Coordinate &);

  Coordinate & ShowContents (const char *, const char *);
  Coordinate & ShowContents (const char *);
};


/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                                 Point                                  |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Point {
private:
  Coordinate coord;
  Coordinate min_coord;
  Coordinate max_coord;
  Coordinate normal;

  double value;

  Point *u;
  Point *ux, *uy;
  Point *intg;
  Point *phi;
  double f;

  string mark;

  double b_u;
  double mp_u;

  double terminalTime;
  double presentTime;
  double dt;

  int timeStep;

  Point *pre, *old;
  Point *hat;

  Point* E;
  Point* EN;
  Point* ES;
  Point* W;
  Point* WN;
  Point* WS;
  Point* N;
  Point* NE;
  Point* NW;
  Point* S;
  Point* SE;
  Point* SW;

  Point* e;
  Point* en;
  Point* es;
  Point* w;
  Point* wn;
  Point* ws;
  Point* n;
  Point* ne;
  Point* nw;
  Point* s;
  Point* se;
  Point* sw;

  Point* uVelocity;
  Point* vVelocity;
  Point* pressure;

  int    index;
  int    mtrxindex;
  int    axis[2];

  char   condition;

public:
  Point ();
  virtual ~Point ();
  Coordinate & Coord () {return this->coord;}
  Coordinate & MinCoord () {return this->min_coord;}
  Coordinate & MaxCoord () {return this->max_coord;}
  Coordinate & Normal () {return this->normal;}
  double MinMaxCoordinate (char, char);
  double Value () const {return this->value;}
  double F () const {return this->f;}
  double Boundaryvalue () const {return this->b_u;}
  double MaterialProperty ();
  double TerminalTime () const {return this->terminalTime;}
  double Time () const {return this->presentTime;}
  double Dt () const {return this->dt;}
  double CalcRightHandSide ();

  int Index () {return this->index;}
  int Mtrx_Index () {return this->mtrxindex;}
  int TimeStep () const {return this->timeStep;}

  char Condition () {return this->condition;}

  string & Mark () {return this->mark;}

  int Axis (char);

  Point * Diff (char);
  Point * Intg () const {return this->intg;}
  Point * Phi () const {return this->phi;}
  Point * Pre () const {return this->pre;}
  Point * Old () const {return this->old;}
  Point * Hat () const {return this->hat;}

  Point * EWNS (char, char);
  Point * EWNS2nd (char, char);

  Point * Velocity (char);
  Point * Pressure () const {return pressure;}

  Point & SetCoordinate (char, double);
  Point & SetCoordinate (double, double);
  Point & SetMinMaxCoordinate (char, char, double);
  Point & SetIndex (int);
  Point & SetMtrx_Index (int);
  Point & SetCondition (char);
  Point & SetBoundaryvalue (double);
  Point & SetMinMax ();
  Point & SetInterfaceMinMax ();
  Point & SetValue (double);
  Point & SetPre (Point*);
  Point & SetOld (Point*);
  Point & SetHat (Point*);
  Point & SetDiff (char, Point*);
  Point & SetIntg (Point*);
  Point & SetPhi (Point*);
  Point & SetF (double);
  Point & SetVelocity (char, Point*);
  Point & SetPressure (Point*);
  Point & SetMark (string);
  Point & SetMaterialProperty (double);
  Point & SetTerminalTime (double);
  Point & SetTime (double);
  Point & SetTimeStep (int);
  Point & SetDt (double);
  Point & SetNormal (char, double);
  Point & SetNormal (Coordinate);
  Point & FindAxialElement (AxialData*, Point*);
  Point & Find2ndAxialElement (AxialData*, Point*);
  Point & FindAxialElementIntp (AxialData*, Point*);
  Point & FindBoundaryElement ();
  Point & IsInterface ();

  bool NextTime ();

  int FindEastAxialLine(AxialData*);
  int FindWestAxialLine(AxialData*);
  int FindNorthAxialLine(AxialData*);
  int FindSouthAxialLine(AxialData*);

  int FindEast2ndAxialLine(AxialData*);
  int FindWest2ndAxialLine(AxialData*);
  int FindNorth2ndAxialLine(AxialData*);
  int FindSouth2ndAxialLine(AxialData*);

  int FindVerticalPoints(AxialData*, int, char);
  int FindHorizontalPoints(AxialData*, int, char);

  void Find_extra_point (FILE*, FILE*, FILE*, FILE*, AxialData*, Point*);
  Point & SetEWNS (char, char, Point*);
  Point & SetEWNS2nd (char, char, Point*);

  bool IsBoundary ();
  bool IsBoundary (const char);

  Point & ExportEWNS (FILE*);
  Point & ExportEWNS2nd (FILE*);

  double PrintDebuggingData (const string, _xData*, _yData*, bool);
  Point & PrintDebuggingData (const string, const char*);
  // Point & CheckRepresentationFormula (Point*, _xData*, _yData*);
  double  CheckRepresentationFormula (Point*, _xData*, _yData*, bool);

  bool operator== (const int);
  bool operator== (const char);
  bool operator== (const Coordinate&);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                                 xyData                                 |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

typedef struct _xData {
  double F;
  double Cu, Cphi, Cf, Cdu;
  double Eu, ENu, ESu;
  double Wu, WNu, WSu;
  double Edu, ENdu, ESdu;
  double Wdu, WNdu, WSdu;
  double Ephi, ENphi, ESphi;
  double Wphi, WNphi, WSphi;
  double Ef, ENf, ESf;
  double Wf, WNf, WSf;
}xData;

typedef struct _yData {
  double F;
  double Cu, Cphi, Cf, Cdu;
  double Nu, NEu, NWu;
  double Su, SEu, SWu;
  double Ndu, NEdu, NWdu;
  double Sdu, SEdu, SWdu;
  double Nphi, NEphi, NWphi;
  double Sphi, SEphi, SWphi;
  double Nf, NEf, NWf;
  double Sf, SEf, SWf;
}yData;

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                             MatrixProcess                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class MatrixProcess {
private:
  double *rb      ; // vector value    [   3*n + i] (      3 * in_pts_num )
  int    *ia      ;
  int    *ja      ;
  double *ent     ;
  int     ent_num ; // total number of the matrix elements

  Point **arrInt;
  double *arrEnt;
  int    *uniqInt;
  double *uniqEnt;
  int    *rowsInt;
  double *rowsEnt;
  int    *sortInt;
  int     Int_num;
  int     ja_num;
  int     matrixSize;

public:
  MatrixProcess (AxialData*);
  virtual ~MatrixProcess ();

  MatrixProcess & initialization (AxialData*, int);
  MatrixProcess & countEnt_num (int, AxialData*, Point*, xData* ,yData*);
  MatrixProcess & MakeMatrixSystem (int, AxialData*, Point*, xData* ,yData*);
  MatrixProcess & calcMatrix (ControlData*, AxialData*, Point*);
  MatrixProcess & ExportMatrixData (AxialData*);
  MatrixProcess & SettingAzimuth (Point*, AxialData*, const int);
  MatrixProcess & printPointDataError (int, Point*);
  MatrixProcess & SettingArray (int, AxialData*, Point*, xData*, yData*);
  MatrixProcess & PrintDebuggingData (AxialData*, Point*, xData* ,yData*);
  MatrixProcess & PrintDebuggingData (AxialData*, Point*, xData* ,yData*, double*);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                         User-Defined functions                         |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

int    IsIn( char, const char * );
int    copyupto( char *, char *, const char * );
int    duplicate( char, const char * );
int    IsRealnum( char * );
int    IsIntnum( char * );
int    all_IsIn( char *, const char * );
char   readingblanks( ifstream &, const char * );
char   readingupto( ifstream &, const char * );
int    getword( char *, ifstream & );
int    getword( char *, ifstream &, const char * );
int    getword( char *, ifstream &, const char *, const char * );
int    getwordupto( char *, ifstream &, const char * );
int    inverseindex( int, int *, int );
int    rmtailblanks( char *, const char * );
void   Make_name( char *, const char *, char * );
void   Make_name( char *, const char * );
void   Make_name( int, char *, char * );
int    CellMemoryCheck( Cell * );
int    AxialLineMemoryCheck (AxialLine *);

/* EllipticSolver */
void EllipticSolver (ControlData*, AxialData*, Point*, string);

/* util.h */
//double형의 두 실수가 같은지 비교 (NearZero값으로 비교)
bool IsEqualDouble (double, double);
// double형의 두 실수가 같은지 비교 (tolerance가 주어진 비교)
bool IsEqualDouble (double, double, double);
// integer set인 memberSet의 0 ~ range안에 target이 있으면 있으면 위치를 return 없으면 -1을 return
int is_member (int, int, int);
// double set인 memberSet의 0 ~ range안에 target이 있으면 있으면 위치를 return 없으면 -1을 return
int is_member (double, double, int);
// 두 점 사이의 거리를 측정
double point_distance (double, double, double, double);
// 반올림을 하는 모듈
double round (double, int);
// 4pt Gaussian quadrature rule integration
double gauss_quadrature (std::function<double (double)>, double, double);
// 해의 relative L2-error를 계산해 주는 모듈
double Calc_Error (AxialData*, Point*);
// 미분의 relative L2-error를 계산해 주는 모듈
double Calc_Derivative_Error (AxialData*, Point*);
// Unbounded domain에서 undetermined coefficient를 계산하는 모듈
double Calc_Undetermined_Coefficient_in_Unbounded_Domian (AxialData*, Point*);
// 에러메시지 출력
void PrintError (const char*);
// 점의 정보를 출력
void PrintPtsInfo (Point*);
// xdat의 정보를 출력
void PrintXdata (xData*);
// ydat의 정보를 출력
void PrintYdata (yData*);
// 점의 오른쪽, 왼쪽, 위쪽, 아래쪽을 대입하는 모듈
void SetAllEWNS (AxialData*, Point*, Point*, Point*);
// 점의 두번재 오른쪽, 왼쪽, 위쪽, 아래쪽을 대입하는 모듈
void SetAllEWNS2nd (AxialData*, Point*, Point*, Point*);
//
Point * HeadVelocity (Point*);

/* Greensfunction.h */
// 그린함수
double greens_function      (double, double, double, double, double, double, int, int, double, double);
// 그린함수의 첫번째 변수로의 미분
double greens_function_t    (double, double, double, double, double, double, int, int, double, double);
// 그린함수의 두번째 변수로의 미분
double greens_function_tau  (double, double, double, double, double, double, int, int, double, double);
// 그린함수의 첫번째, 두뻔째 변수로의 미분
double greens_function_ttau (double, double, double, double, double, double, int, int, double, double);
// 그린함수와 basis function의 적분
double greens_integral     (int, double, double, double, double, double, int, int, double, double);
// 그린함수의 두번째변수로의 미분과 basis function의 적분
double greens_integral_tau (int, double, double, double, double, double, int, int, double, double);

/* CalcRepresenSol.h */
// 해의 표현식의 계수를 계산
void CalcRepresenCoef (Point*, xData*, yData*, bool);
// Neumann 경계점에서의 해의 표현식의 계수를 계산
void CalcRepresenCoefAtNeumann (Point*, xData*, yData*);
// 경계근처에서의 해의 표현식을 처리
void TransposeBoundaryData (Point*, xData*, yData*, bool);
// 경계에서의 phi값을 부여
void AssignPhivalue (AxialData*, Point*);
// 경계에서의 phi값을 부여
int AssignPhi (Point*);
// Neumann 경계와 Interface가 만나는지 확인
bool SetInterfaceBoundary (Point*, char);
// 우변의 항 f를 계산
void RightHandSide (Point*, xData*, yData*);
void RightHandSideDirichlet (Point*, Point*, xData*, yData*);
void RightHandSideNeumann (Point*, xData*, yData*);
// 축선의 양끝점중 한 개의 점의 경계조건이 Infinity 혹은 Singularity일 때, 정확한 적분을 계산
void ExactIntegration(char, xData*, yData*, int, int, int, double, double, double, double, double, double);
// 해의 표현식의 계수의 확인
double CheckRepresentationFormula (Point*, xData*, yData*, int, double, double, double, double, double, double, int, int);
// y-축선에서의 u의 계수를 계산
double CalcVerticalUCoefficient (char, double, double, double, double, int, double);
// x-축선에서의 u의 계수를 계산
double CalcHorizontalUCoefficient (char, double, double, double, double, int, double);
// y-축선에서의 phi의 계수를 계산
double CalcVerticalPHICoefficient (char, double, double, double, double, int, double);
// x-축선에서의 phi의 계수를 계산
double CalcHorizontalPHICoefficient (char, double, double, double, double, int, double);
// y-축선에서의 f의 값을 계산
double CalcVerticalFCoefficient (char, double, double, double, double, int, double);
// x-축선에서의 f의 값을 계산
double CalcHorizontalFCoefficient (char, double, double, double, double, int, double);
// y-축선에서의 u의 계수를 계산
double CalcVerticalDiffUCoefficient (char, double, double, double, double, int, double);
// x-축선에서의 u의 계수를 계산
double CalcHorizontalDiffUCoefficient (char, double, double, double, double, int, double);
// y-축선에서의 phi의 계수를 계산
double CalcVerticalDiffPHICoefficient (char, double, double, double, double, int, double);
// x-축선에서의 phi의 계수를 계산
double CalcHorizontalDiffPHICoefficient (char, double, double, double, double, int, double);
// y-축선에서의 f의 값을 계산
double CalcVerticalDiffFCoefficient (char, double, double, double, double, int, double);
// x-축선에서의 f의 값을 계산
double CalcHorizontalDiffFCoefficient (char, double, double, double, double, int, double);
//
void SettingMaterialProperties (Point*, double*, double*, double*, double*);
//
void InitializationCoef (xData*, yData*, string);
//
void SettingGreenShape (Point*, int*, int*);
//
void SettingCoefficient (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double, const int, const bool);
//
void SettingNeumannCoefficient (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double, const bool);
//
void SettingDirichletCoefficient (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double, const int, const bool);
//
void SettingDiffCoefficient (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double);
//
void SettingDiffDiffCoefficient (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double);
//
void SettingTimeIntegration (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double, const bool);
//
void SettingDiffTimeIntegration (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double);
//
void SettingNeumannTimeIntegration (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double, const bool);
//
void SettingTimeIntegrationRightHand (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double, const bool);
//
void SettingDiffTimeIntegrationRightHand (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double);
//
void SettingDiffDiffTimeIntegrationRightHand (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double);
//
void SettingNeumannTimeIntegrationRightHand (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double, const bool);
//
void SettingNeumannTimeIntegrationRightHandCoefficient (Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double, const bool);
//
void SettingDirichletTimeIntegrationRightHand (Point*, Point*, xData*, yData*, const double, const double, const double, const double, const double, const double, const int, const int, const double, const double, const double, const double, const bool);
//
void CalcInterfaceCoef (Point*, xData*, yData*, const char*, const double, const double, const double, const double, const double);
//
void CalcNeumannCoef (Point*, xData*, yData*, const char*, const double, const double, const double, const double, const double);
//
void ApproximateInterfacePt (Point*, xData*, yData*, double, double ,double, double, double ,double ,double, double, double, double, bool);
//
void ApproximateDiffDiffPt (Point*, xData*, yData*, double, double ,double, double, double ,double ,double, double, double, double, bool);
//
void TransposeNeumannBoundaryData (Point*, xData*, yData*);
//
void TransposeOtherBoundaryData (Point*, xData*, yData*);
/* CalcAtNeumannPt */
// Neumann 경계조건이 주어진 점에서의 해의 값을 계산하는 모듈
double CalcAtNeumannpt (char, Point*, xData*, yData*);
// Neumann 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산하는 모듈
void CalcNeumannCoef (Point*, xData*, yData*);
// Infinity 경계조건이 주어진 점에서의 해의 값을 계산하는 모듈
double CalcAtInfinitypt (char, Point*, xData*, yData*);
// Infinity 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산하는 모듈
void CalcInfinityCoef (Point*, xData*, yData*);


/* 경계점에서의 미분의 값을 계산하는 모듈 */
double CalcDiffAtBoundary (char, Point*, xData*, yData*);

/* 경계에서의 미분값을 계산하기 위한 모듈 */
void CalcBoundaryDiffCoef (Point*, xData*, yData*);


// 축선의 양끝점중 한 개의 점의 경계조건이 Infinity 혹은 Singularity일 때, 정확한 미분의 적분을 계산
void ExactDiffIntegration(char, xData*, yData*, int, int, int, double, double, double, double, double, double);

double greens_function (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2);
double greens_function_t (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2);
double greens_function_tau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2);
double greens_function_ttau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2);
double greens_integral (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2);
double greens_integral_t (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2);
double greens_integral_tau (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2);
double greens_integral_ttau (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2);
double greens_coefficient_t (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2);
double greens_coefficient_ttau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2);

#endif
