#ifndef STUFF_H
#define STUFF_H

#include "Basic++.hpp"

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              Line2D Class                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

Line2D::Line2D()
:CP_q(D2)
{
  start.setvector(D2, ZeroValue,ZeroValue);
  end.setvector(D2, ZeroValue,ZeroValue);

  this->CalculateProperties();
}

Line2D::Line2D(Vector & sv, Vector & ev)
:CP_q(D2)
{
  start = sv;
  end = ev;

  this->CalculateProperties();
}

Line2D::Line2D(Vector & sv, Vector & ev, char bdc, double bdv)
:CP_q(D2)
{
  start = sv;
  end = ev;
  cond = bdc;
  value = bdv;

  this->CalculateProperties();
}

Line2D::Line2D(double xs, double ys, double xe, double ye)
:CP_q(D2)
{
  start.setvector(D2, xs, ys);
  end.setvector(D2, xe, ye);

  this->CalculateProperties();
}

Line2D::Line2D(double xs, double ys, double xe, double ye, char bdc, double bdv)
:CP_q(D2)
{
  start.setvector(D2, xs, ys);
  end.setvector(D2, xe, ye);
  cond = bdc;
  value = bdv;

  this->CalculateProperties();
}

Line2D::~Line2D(){}

Line2D &  Line2D::SetLength()
{
  Len = (end-start).norm();

  return  *this;
}

Line2D &  Line2D::SetOPRT(Vector & qvec)
{
  if(qvec.norm() < NearZero) CP_q.getmatrix()->setInxn(D2);
  else                       CP_q.COP_onto(qvec);

  return *this;
}

Line2D &  Line2D::SetOPRT()
{
  this->SetOPRT(this->GetTangent());

  return *this;
}

//  Line2D-this(master) <---- Line2D(slave)
//  <--(projecting onto the normal surface to the master line)

int  Line2D::IsCross(Line2D & slave, Vector & vecCROSS, double & cptPOS)
{
  Vector   ref=this->GetStart();
  Vector   sL=slave.GetStart()-ref;
  Vector   eL=slave.GetEnd()-ref;
  Vector   tau(D2);
  Vector   ps(D2), pe(D2);
  double   s_len, e_len, t=ZeroValue;

  // sL= Xs-MXs: MXs is master start point
  // eL= Xe-MXs

  // projecting sL and eL onto the plane perpendicular to the master Line2D-<this>
  ps = CP_q.acting(sL);
  pe = CP_q.acting(eL);

  s_len = ps.norm();
  e_len = pe.norm();

  //     if((pe-ps).norm() < NearZero) { cout << "/"; return  0; } // parallel
  if((pe-ps).norm() < NearZero) { return  0; } // parallel

  //     if(ps*pe > NearZero) { cout << "|"; return 0; } // on one side
  if(ps*pe > NearZero) { return 0; } // on one side
  // if(ps*pe > -NearZero) { return 0; } // on one side // THIS IS TEST BY JUNHONG, JO

  vecCROSS.setzerovector(D2); // set zero vector(important)

  // Across anyway and find the cross point
  sL = slave.GetStart();
  eL = slave.GetEnd();
  vecCROSS.add(sL.scalar(e_len/(s_len+e_len)));
  vecCROSS.add(eL.scalar(s_len/(s_len+e_len)));

  // Meet or not
  tau=this->GetTangent();
  tau.scalar(UnitValue/this->GetLength());
  t = (vecCROSS-ref)*tau;

  //     if(t*(t-UnitValue) > ZeroValue) { cout << "-"; return 0;} // meet away
  if(t*(t-UnitValue) > ZeroValue) { return 0;} // meet away

  cptPOS = t;

  // cout << "X"; // by junhong, jo

  return 1;
}

/*
int  Line2D::ScanCross(int nSlv, Line2D *slave, Vector** & crosspts)
{
int      i, k, count, nFind;
char     tstr[128];
Vector   *vtmp;
double   *tindex=NULL;
double   **tp=NULL;
double   *ttmp;

tindex = new double[nSlv];
if(!tindex) { cout << "not memory for [tindex]..." << endl; exit(0); }

cout << "<" << nSlv << ">";
for(nFind = 0, i = 0 ; i < nSlv ; i++) {
k = this->IsCross(slave[i], crosspts[nFind], tindex[nFind]);
nFind += k;
}
cout << "<" << nFind << ">" << endl;

// sorting
tp = new double*[nFind];
if(!tp) { cout << "not memory for [tp**]..." << endl; exit(0); }

for(i = 0 ; i < nFind ; i++) tp[i] = tindex+i;

// using pointer changes
for(i = 0 ; i < nFind ; i++)
for(k = 1 ; k < nFind-i ; k++)
if(*(tp[k-1]) > *(tp[k])) {
ttmp = tp[k-1]; tp[k-1] = tp[k]; tp[k] = ttmp;
vtmp  = crosspts[k-1];
crosspts[k-1] = crosspts[k];
crosspts[k] = vtmp;
}

// eliminate duplicated
//       for(k = 1 ; k < nFind ; k++)
//         if(fabs(*(tp[k]) - *(tp[k-1])) < NearZero) {
//                 delete crosspts[k-1];
//                 crosspts[k-1] = NULL;
//         }

// free memory
for(i = 0 ; i < nFind ; i++) tp[i]=NULL;
delete [] tp;
delete [] tindex;

return nFind;
}
*/

int  Line2D::ScanCross(int nSlv, Line2D *slave, Thing & crsPT, double sign)
{
  int      i, j, nFind;
  Vector   crsvec(D2); // crsvec = [ x y ]
  Vector   omnivec(D2+1); // omnivec = [ x y : t]
  Vector  &cpt=crsvec;
  double   tPOS;

  // cout << "<" << nSlv << ">"; // by junhong, jo
  for(nFind = 0, i = 0 ; i < nSlv ; i++)
  if(this->IsCross(slave[i], cpt, tPOS)) {
    // make omnivec = [ crsvec : tPOS ]
    for(j = 0 ; j < D2 ; j++) omnivec[j] = round(crsvec[j], Margin);
    omnivec[D2] = tPOS;

    // save [ crsvec : tPOS ] to crsPT
    if      (this->GetTangent() * slave[i].GetNormal() * sign < ZeroValue) crsPT.MakeCell(omnivec, slave[i].GetBDcondition(), slave[i].GetBDvalue(), 'I', slave[i].GetNormal());
    else if (this->GetTangent() * slave[i].GetNormal() * sign > ZeroValue) crsPT.MakeCell(omnivec, slave[i].GetBDcondition(), slave[i].GetBDvalue(), 'O', slave[i].GetNormal());

    nFind++;
  }

  // cout << "<" << nFind << ">" << endl; // by junhong, jo

  return nFind;
}

Line2D &  Line2D::SetNormal()
{
  if(Len > NearZero) {
    Tangent = (end - start)*(UnitValue/Len);
    Normal.setvector(D2, Tangent[1], -Tangent[0]);
  } else {
    Tangent.setvector(D2, ZeroValue,ZeroValue);
    Normal.setvector(D2, ZeroValue,ZeroValue);
  }

  return  *this;
}

Line2D &  Line2D::SetNormal(Vector &src)
{
  Tangent = (end - start)*(UnitValue/Len);
  Normal.setvector(src);

  return  *this;
}

Line2D &  Line2D::SetNormal(double nx, double ny)
{
  Tangent = (end - start)*(UnitValue/Len);
  Normal.setvector(D2, nx, ny);

  return  *this;
}


int   Line2D::Inside(Vector & vec)
{
  double  s_len, e_len;

  if(vec.getdim() != D2) {
    cout << "Inconsistant vector dimension in LINE2D::Inside()..." << endl;
    exit(0);
  }

  s_len=(this->GetStart()-vec).norm();
  e_len=(this->GetEnd()-vec).norm();

  if((s_len+e_len)-this->GetLength() > NearZero) return 0;

  return  1;
}

Line2D &  Line2D::CalculateProperties()
{
  this->SetLength();
  this->SetNormal();
  this->SetOPRT();

  return  *this;
}

Vector   Line2D::MapFromReference(double t)
{
  Vector  tv(D2);

  if(t < ZeroValue || t > UnitValue) {
    cout << "Out of range in t-variable!\n";
    exit(-1);
  }

  tv.setvalue(0, (UnitValue-t)*start[0] + t*end[0]);
  tv.setvalue(1, (UnitValue-t)*start[1] + t*end[1]);

  return tv;
}

double    Line2D::Distance(Vector & vec)
{
  Vector  pvec=vec-this->GetStart();

  pvec = this->GetOPRT().acting(pvec);

  return  pvec.norm();
}

Line2D & Line2D::SetLine(Vector & spt, Vector & ept)
{
  start.setvector(spt);
  end.setvector(ept);

  this->CalculateProperties();

  return *this;
}

Line2D & Line2D::SetLine(double xs, double ys, double xe, double ye)
{
  start.setvector(D2, xs, ys);
  end.setvector(D2, xe, ye);

  this->CalculateProperties();

  return  *this;
}

Line2D & Line2D::SetLine(Vector & spt, Vector & ept, char bdc, double bdv)
{
  start.setvector(spt);
  end.setvector(ept);

  cond = bdc;
  value = bdv;

  this->CalculateProperties();

  return *this;
}

Line2D & Line2D::SetLine(double xs, double ys, double xe, double ye, char bdc, double bdv)
{
  start.setvector(D2, xs, ys);
  end.setvector(D2, xe, ye);

  cond = bdc;
  value = bdv;

  this->CalculateProperties();

  return  *this;
}

Line2D & Line2D::SetBD(char bdc, double bdv)
{
  cond = bdc;
  value = bdv;

  this->CalculateProperties();

  return  *this;
}

Vector  &  Line2D::operator[](int offset)
{
  if(offset == 0) return start;
  else              return end;
}

Line2D &  Line2D::showcontents(const char *str)
{
  cout << str;
  this->GetStart().showcontents(" "," ");
  cout << "+--|" << Len << "|-->";
  this->GetEnd().showcontents(" ","::");

  this->GetTangent().showcontents("<T>"," ");
  this->GetNormal().showcontents("<N>");

  return *this;
}

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                          Triangle in 3D Class                          |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

Triangle::Triangle()
:index(-1),CP_n(D3)
{
  int   i;

  for(i = 0 ; i < D3 ; i++)
  Vertex[i].setvector(D3, ZeroValue, ZeroValue, ZeroValue);

  this->CalculateProperties();
}

Triangle::Triangle(Vector & v1, Vector & v2, Vector & v3)
:index(-1),CP_n(D3)
{
  Vertex[0] = v1;
  Vertex[1] = v2;
  Vertex[D3-1] = v3;

  this->CalculateProperties();
}

Triangle::Triangle(Vector & v1, Vector & v2, Vector & v3, char bdc, double bdv)
:index(-1),CP_n(D3)
{
  Vertex[0] = v1;
  Vertex[1] = v2;
  Vertex[D3-1] = v3;

  cond = bdc;
  value = bdv;

  this->CalculateProperties();
}

Triangle::~Triangle()
{}

  Triangle  &  Triangle::SetArea()
  {
    Vector  u=Vertex[1]-Vertex[0];
    Vector  v=Vertex[D3-1]-Vertex[0];

    Area = (u.cross(v).norm())/2.0;

    return  *this;
  }

  Triangle & Triangle::SetTriangle(Vector & v1, Vector & v2, Vector & v3)
  {
    Vertex[0] = v1;
    Vertex[1] = v2;
    Vertex[D3-1] = v3;

    this->CalculateProperties();

    return  *this;
  }

  Triangle & Triangle::SetTriangle(Vector & v1, Vector & v2, Vector & v3, char bdc, double bdv)
  {
    Vertex[0] = v1;
    Vertex[1] = v2;
    Vertex[D3-1] = v3;

    cond = bdc;
    value = bdv;

    this->CalculateProperties();

    return  *this;
  }

  Triangle & Triangle::SetBD(char bdc, double bdv)
  {

    cond = bdc;
    value = bdv;

    this->CalculateProperties();

    return  *this;
  }

  Triangle & Triangle::PutVertex(int offset, Vector & vec)
  {
    if(vec.getdim() != D3) {
      cout << "Dimenssion of [vec] is not 3..." << endl;
      exit(0);
    }

    Vertex[offset] = vec;

    return  *this;
  }

  Triangle  &  Triangle::SetNormal()
  {
    Vector  u=Vertex[1]-Vertex[0];
    Vector  v=Vertex[D3-1]-Vertex[0];

    if(Area > NearZero) {
      Normal=u.cross(v);
      Normal.setunit();
    } else  Normal.setvector(D3, ZeroValue, ZeroValue, ZeroValue);

    return  *this;
  }

  Triangle  &  Triangle::SetNormal(Vector &src)
  {
    Normal.setvector(src);

    return  *this;
  }

  Triangle  &  Triangle::PutNormal(Vector & src)
  {
    if(src.getdim() != D3) {
      cout << "Dimenssion of [src] is not 3..." << endl;
      exit(0);
    }

    Normal = src;

    return  *this;
  }

  Triangle &  Triangle::SetOPRT(Vector & qvec)
  {
    if(qvec.norm() < NearZero) CP_n.getmatrix()->setInxn(D3);
    else                         CP_n.COP_onto(qvec);

    return *this;
  }

  Triangle &  Triangle::SetOPRT()
  {
    this->SetOPRT(this->GetNormal());

    return *this;
  }

  Triangle &  Triangle::CalculateProperties()
  {
    this->SetArea();
    this->SetNormal();
    this->SetOPRT();

    return  *this;
  }

  int   Triangle::Inside(Vector & vec)
  {
    int     i;
    double  areaSide;
    Vector  Stri[D3];

    if(vec.getdim() != D3) {
      cout << "Inconsistant vector dimension in Triangle::Inside()..." << endl;
      exit(0);
    }

    for(i = 0 ; i < D3 ; i++) Stri[i] = this->GetVertex(i)-vec;

    for(areaSide = ZeroValue, i = 0 ; i < D3 ; i++)
    areaSide += Stri[i].cross(Stri[(i+1)%D3]).norm()/2.0;

    if(areaSide-this->GetArea() > NearZero) return 0;

    return  1;
  }

  Triangle & Triangle::showcontents(const char * str)
  {
    cout << str;
    cout << "(" << index << ")" << "[AREA =" << Area << "]" << endl;
    Vertex[0].showcontents("X[0]");
    Vertex[1].showcontents("X[1]");
    Vertex[2].showcontents("X[2]");

    Normal.showcontents("N = ");

    return *this;
  }

  Vector  &  Triangle::operator[](int offset)
  {
    if(offset == 0 || offset == 1) return Vertex[offset];
    else                             return Vertex[D3-1];
  }

  /****************************************************************************/
  /*+------------------------------------------------------------------------+*/
  /*|                              Line3D Class                              |*/
  /*+------------------------------------------------------------------------+*/
  /****************************************************************************/

  Line3D::Line3D()
  :CP_q(D3)
  {
    start.setvector(D3, ZeroValue,ZeroValue, ZeroValue);
    end.setvector(D3, ZeroValue,ZeroValue, ZeroValue);

    this->CalculateProperties();
  }

  Line3D::Line3D(Vector & sv, Vector & ev)
  :CP_q(D3)
  {
    start = sv;
    end = ev;

    this->CalculateProperties();
  }

  Line3D::Line3D(double xs, double ys, double zs, double xe, double ye, double ze)
  :CP_q(D3)
  {
    start.setvector(D3, xs, ys, zs);
    end.setvector(D3, xe, ye, ze);

    this->CalculateProperties();
  }

  Line3D::~Line3D(){}

  Line3D &  Line3D::SetLength()
  {
    Len = (end-start).norm();

    return  *this;
  }

  Line3D &  Line3D::SetOPRT(Vector & qvec)
  {
    if(qvec.norm() < NearZero) CP_q.getmatrix()->setInxn(D3);
    else                         CP_q.COP_onto(qvec);

    return *this;
  }

  Line3D &  Line3D::SetOPRT()
  {
    this->SetOPRT(this->GetTangent());

    return *this;
  }

  //  Line3D-this(master) <---- Line3D(slave)
  //  <--(projecting onto the normal surface to the master line)

  int  Line3D::IsCross(Triangle & slave, Vector & vecCROSS, double & cptPOS)
  {
    Vector   ref=this->GetStart();
    Vector   vtx[D3], pvtx[D3];
    Vector   tau(D3);
    Triangle Ptri;
    double   areaPtri, t=ZeroValue;
    double   Parea[D3];
    int      i;

    // vtx[i] = X(i)-MXs: MXs is master start point
    for(i = 0 ; i < D3 ; i++) vtx[i] = slave.GetVertex(i)-ref;

    // projecting T-vertices onto the plane perpendicular to the master Line3D-<this>
    for(i = 0 ; i < D3 ; i++) pvtx[i] = CP_q.acting(vtx[i]);

    // projected triangle to the master Line3D-<this>
    Ptri.SetTriangle(pvtx[0], pvtx[1], pvtx[2]);

    //     if(Ptri.GetArea() < NearZero) { cout << "/"; return  0; } // parallel
    if(Ptri.GetArea() < NearZero) { return  0; } // parallel

    // Area of the side three surfaces by the projected triangle and the origin
    for(areaPtri = ZeroValue, i = 0 ; i < D3 ; i++)
    areaPtri += Parea[(i+2)%D3] = Ptri[i].cross(Ptri[(i+1)%D3]).norm()/2.0;

    //     if(areaPtri-Ptri.GetArea() > NearZero) { cout << "|"; return 0;} // Away from M-Line
    if(areaPtri-Ptri.GetArea() > NearZero) { return 0;} // Away from M-Line

    vecCROSS.setzerovector(D3);  // set zero vector(important)

    //  Find the cross point inside the triangle
    for(i = 0 ; i < D3 ; i++) {
      vtx[i] = slave.GetVertex(i);
      vecCROSS.add(vtx[i].scalar(Parea[i]/Ptri.GetArea()));
    }
    // Meet or not
    tau=this->GetTangent();
    tau.scalar(UnitValue/this->GetLength());
    t = (vecCROSS-ref)*tau;

    //     if(t*(t-UnitValue) > ZeroValue) { cout << "-"; return 0;} // meet away
    if(t*(t-UnitValue) > ZeroValue) { return 0;} // meet away

    cptPOS = t;

    cout << "X";

    return 1;
  }

  /*
  int  Line3D::ScanCross(int nSlv, Triangle *slave, Vector** & crosspts)
  {
  int      i, k, count, nFind;
  char     tstr[128];
  Vector   *vtmp;
  double   *tindex=NULL;
  double   **tp=NULL;
  double   *ttmp;

  tindex = new double[nSlv];
  if(!tindex) { cout << "not memory for [tindex]..." << endl; exit(0); }

  cout << "<" << nSlv << ">";
  for(nFind = 0, i = 0 ; i < nSlv ; i++) {
  k = this->IsCross(slave[i], crosspts[nFind], tindex[nFind]);
  nFind += k;
}
cout << "<" << nFind << ">" << endl;

// sorting
tp = new double*[nFind];
if(!tp) { cout << "not memory for [tp**]..." << endl; exit(0); }

for(i = 0 ; i < nFind ; i++) tp[i] = tindex+i;

// using pointer changes
for(i = 0 ; i < nFind ; i++)
for(k = 1 ; k < nFind-i ; k++)
if(*(tp[k-1]) > *(tp[k])) {
ttmp = tp[k-1]; tp[k-1] = tp[k]; tp[k] = ttmp;
vtmp  = crosspts[k-1];
crosspts[k-1] = crosspts[k];
crosspts[k] = vtmp;
}

// eliminate duplicated
//       for(k = 1 ; k < nFind ; k++)
//         if(fabs(*(tp[k]) - *(tp[k-1])) < NearZero) {
//               delete  crosspts[k-1];
//               crosspts[k-1] = NULL;
//         }

// free memory
for(i = 0 ; i < nFind ; i++) tp[i]=NULL;
delete [] tp;
delete [] tindex;

return nFind;
}
*/

int  Line3D::ScanCross(int nSlv, Triangle *slave, Thing  & crsPT)
{
  int      i, j, nFind;
  Vector   crsvec(D3); // crsvec = [ x y z ]
  Vector   omnivec(D3+1); // omnivec = [ crsvec : t ]
  Vector  &cpt=crsvec;
  double   tPOS;

  crsPT.GoTail(); // place [run] to the tail cell

  cout << "<" << nSlv << ">";
  for(nFind = 0, i = 0 ; i < nSlv ; i++)
  if(this->IsCross(slave[i], cpt, tPOS)) {
    // make omnivec = [ crsvec : tPOS ]
    for(j = 0 ; j < D3 ; j++) omnivec[j] = cpt[j];
    omnivec[D3] = tPOS;

    // save omnivec to crsPT
    crsPT.MakeCell(omnivec);

    nFind++;
  }

  cout << "<" << nFind << ">" << endl;

  return nFind;
}

Line3D &  Line3D::SetNormal()
{
  if(Len > NearZero) {
    Tangent = (end - start)*(UnitValue/Len);
    Normal.setvector(D3, ZeroValue, ZeroValue, ZeroValue);
  } else {
    Tangent.setvector(D3, ZeroValue, ZeroValue, ZeroValue);
    Normal.setvector(D3, ZeroValue, ZeroValue, ZeroValue);
  }

  return  *this;
}

Line3D &  Line3D::CalculateProperties()
{
  this->SetLength();
  this->SetNormal();
  this->SetOPRT();

  return  *this;
}

int   Line3D::Inside(Vector & vec)
{
  double  s_len, e_len;

  if(vec.getdim() != D3) {
    cout << "Inconsistant vector dimension in LINE3D::Inside()..." << endl;
    exit(0);
  }

  s_len=(this->GetStart()-vec).norm();
  e_len=(this->GetEnd()-vec).norm();

  if((s_len+e_len)-this->GetLength() > NearZero) return 0;

  return  1;
}

Vector   Line3D::MapFromReference(double t)
{
  Vector  tv(D3);

  if(t < ZeroValue || t > UnitValue) {
    cout << "Out of range in t-variable!\n";
    exit(-1);
  }

  tv.setvalue(0, (UnitValue-t)*start[0] + t*end[0]);
  tv.setvalue(1, (UnitValue-t)*start[1] + t*end[1]);

  return tv;
}

double    Line3D::Distance(Vector & vec)
{
  Vector  pvec = vec - this->GetStart();

  pvec = this->GetOPRT().acting(pvec);

  return  pvec.norm();
}

Line3D & Line3D::SetLine(Vector & spt, Vector & ept)
{
  start.setvector(spt);
  end.setvector(ept);

  this->CalculateProperties();

  return *this;
}

Line3D & Line3D::SetLine(double xs, double ys, double zs, double xe, double ye, double ze)
{
  start.setvector(D3, xs, ys, zs);
  end.setvector(D3, xe, ye, ze);

  this->CalculateProperties();

  return  *this;
}

Vector  &  Line3D::operator[](int offset)
{
  if(offset == 0) return start;
  else              return end;
}

Line3D &  Line3D::showcontents(const char *str)
{
  cout << str;
  this->GetStart().showcontents(" "," ");
  cout << "+--|" << Len << "|-->";
  this->GetEnd().showcontents(" ","::");

  this->GetTangent().showcontents("<T>"," ");
  this->GetNormal().showcontents("<N>");

  return *this;
}

#endif
