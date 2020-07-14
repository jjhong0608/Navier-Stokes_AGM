#ifndef BASIC_H
#define BASIC_H

#include "Header.hpp"

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              Vector Class                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

Vector::Vector()
:dim(0)
{
  comp = NULL;
}

Vector::Vector(int dms)
:dim(dms)
{
  int  i;

  comp = new double[dms];
  if(!comp) {
    cout << "no memory for " << dms << "D vector..." << endl;
    exit(1);
  }

  for(i = 0 ; i < dim ; i++) comp[i] = ZeroValue;
}

Vector::Vector(double x, double y)
:dim(D2)
{
  comp = new double[D2];
  if(!comp) {
    cout << "2D vector constructor(no memory)..." << endl;
    exit(1);
  }

  comp[0] = x;
  comp[1] = y;
}

Vector::Vector(double x, double y, double z)
:dim(D3)
{
  comp = new double[D3];
  if(!comp) {
    cout << "3D vector constructor(no memory)..." << endl;
    exit(1);
  }

  comp[0] = x;
  comp[1] = y;
  comp[2] = z;
}

Vector::Vector(const Vector & src)
{
  int  i;

  dim = src.getdim();

  comp = new double[dim];
  if(!comp){
    cout << dim << "D vector  memory reservation fails at Vector constructor." << endl;
    exit(1);
  }

  for(i = 0 ; i < dim ; i++) comp[i] = src.getvalue(i);
}

Vector::~Vector()
{
  if(comp) delete [] comp;
}

double  Vector::getvalue(int k) const
{
  if(k < dim)  return comp[k];
  else           return ZeroValue;
}

Vector &  Vector::setvalue(int k, double x)
{
  if(k >= 0 && k < dim) comp[k] = x;

  return  *this;
}

Vector &  Vector::setvector(int nitem, ...)
{
  va_list  arguments;
  double   value;
  int  i;

  if(nitem <= 0) {
    cout << "Check out that 1-st argument is positive integer..." << endl;

    if(comp) delete [] comp;

    dim = 0; comp = NULL;

    return *this;
  }

  if(comp != NULL && this->getdim() != nitem) delete [] comp;

  dim = nitem; // reset new dimension
  comp = new double[nitem];

  va_start(arguments, nitem);
  for(i = 0 ; i < nitem ; i++){
    value = va_arg(arguments, double);
    this->setvalue(i, value);
  }
  va_end(arguments);

  return  *this;
}

Vector &  Vector::setvector(int dms, double *x)
{
  int   i;

  if(this->getdim() < dms) {
    if(comp) delete [] comp;

    dim = dms;
    comp = new double[dms];
  }

  for(i = 0 ; i < dms ; i++) this->setvalue(i, x[i]);
  for(i = dms ; i < dim ; i++) this->setvalue(i, ZeroValue);

  return  *this;
}

Vector &  Vector::setvector(Vector & src)
{
  int    i;

  if(this->getdim() != src.getdim()){
    if(comp) delete [] comp;

    dim = src.getdim();
    comp = new double[dim];
  }

  for(i = 0 ; i < dim ; i++) comp[i] = src.getvalue(i);

  return  *this;
}

Vector &  Vector::setunit()
{
  double    alpha;

  alpha = this->norm();
  if(alpha < NearZero) {
    cout << "Almost zero vector in [Vector::setunit]..." << endl;
    exit(1);
  }

  this->scalar(UnitValue/alpha);

  return  *this;
}

Vector  Vector::getunit()
{
  double    alpha;
  Vector    tvec;

  tvec = *this;

  alpha = tvec.norm();
  if(alpha < NearZero) {
    cout << "Almost zero vector in [Vector::getunit]..." << endl;
    exit(1);
  }
;
  tvec.scalar(UnitValue/alpha);

  return  tvec;
}

Vector &  Vector::cut(int dim_new)
{
  int      i, dim_org=this->getdim();
  Vector   *tmp;

  if(dim_new <= 0) {
    if(comp) delete [] comp;

    dim = 0;
    comp = NULL;

    return *this;
  } else if(dim_new == dim_org) return *this;

  tmp = new Vector(*this); // copy this to tmp

  if(comp) delete [] comp;

  comp = new double[dim_new];
  dim = dim_new;

  for(i = 0 ; i < dim_new ; i++)
  this->setvalue(i, tmp->getvalue(i)); // restore tmp to this

  delete tmp;

  return  *this;
}

Vector &  Vector::setzerovector()
{
  int    i;

  for(i = 0 ; i < this->getdim() ; i++) comp[i] = ZeroValue;

  return  *this;
}

Vector &  Vector::setzerovector(int n)
{
  int    i;

  if(this->getdim() != n){
    if(comp) delete [] comp;

    dim = n;
    comp = new double[n];
  }

  for(i = 0 ; i < this->getdim() ; i++) comp[i] = ZeroValue;

  return  *this;
}

Vector &  Vector::setonevector()
{
  int    i;

  for(i = 0 ; i < this->getdim() ; i++) comp[i] = UnitValue;

  return  *this;
}

Vector &  Vector::setonevector(int n)
{
  int    i;

  if(this->getdim() != n){
    if(comp) delete [] comp;

    dim = n;
    comp = new double[n];
  }

  for(i = 0 ; i < this->getdim() ; i++) comp[i] = UnitValue;

  return  *this;
}

Vector &  Vector::setEk(int ndim, int off)
{
  int    i;

  if(this->getdim() != ndim){
    if(comp) delete [] comp;

    dim = ndim;
    comp = new double[ndim];
  }

  for(i = 0 ; i < ndim ; i++) comp[i] = ZeroValue;

  comp[off] = UnitValue;

  return  *this;
}

Vector & Vector::scalar(double sc)
{
  int      i;
  double    crd;

  for(i = 0 ; i < dim ; i++) {
    crd = this->getvalue(i);
    this->setvalue(i, sc*crd);
  }

  return *this;
}

Vector & Vector::add(const Vector & vec)
{
  int       i;
  double    crd;

  for(i = 0 ; i < dim ; i++){
    crd = this->getvalue(i);
    this->setvalue(i, crd+vec.getvalue(i));
  }

  return *this;
}

double   Vector::dot(Vector & vec)
{
  double sum;
  int i;

  for(sum = ZeroValue, i = 0 ; i < dim ; i++) sum += comp[i]*vec.getvalue(i);

  return sum;
}

Vector  Vector::cross(Vector & vec)
{
  Vector   tmp(D3);
  double   c1, c2, c3;

  c1 = this->getvalue(1)*vec.getvalue(2)-this->getvalue(2)*vec.getvalue(1);
  c2 = this->getvalue(2)*vec.getvalue(0)-this->getvalue(0)*vec.getvalue(2);
  c3 = this->getvalue(0)*vec.getvalue(1)-this->getvalue(1)*vec.getvalue(0);

  tmp.setvalue(0, c1);
  tmp.setvalue(1, c2);
  tmp.setvalue(2, c3);

  return  tmp;
}

double  Vector::triple(Vector & b_vec, Vector & c_vec) // this*(bxC)
{
  Vector   tmp(D3);

  tmp = b_vec.cross(c_vec);

  return  this->dot(tmp);
}

double   Vector::norm()
{
  double sum, value;
  int i;

  for(sum = ZeroValue, i = 0 ; i < dim ; i++) {
    value = this->getvalue(i);
    sum += value*value;
  }

  return sqrt(sum);
}

double & Vector::operator[](int off)
{
  if(off >= 0 && off < this->getdim()) return comp[off];
  else {
    cout << "no reference at [" << off << "] for ";
    cout << this->getdim() << "D vector..." << endl;
    exit(1);
  }
}

Vector  Vector::operator+(Vector & post)
{
  Vector   tmp(*this);

  tmp.add(post);

  return tmp;
}

Vector  Vector::operator-(Vector & post)
{
  Vector   tmp(post*(-UnitValue));

  tmp.add(*this);

  return tmp;
}

Vector  Vector::operator*(double sc)
{
  Vector   tmp(*this);

  tmp.scalar(sc);

  return   tmp;
}

double  Vector::operator*(Vector & rhs)
{
  int      i, cdim;
  double   value;

  value = ZeroValue;

  cdim = (this->getdim() < rhs.getdim()) ? this->getdim():rhs.getdim();

  for(i = 0 ; i < cdim ; i++) value += (*this)[i]*rhs[i];

  return   value;
}

Vector & Vector::operator=(const Vector & src)
{
  int  i;

  if(this == &src) return *this;

  if(src.getdim() != this->getdim()) {
    if(comp) {delete [] comp; comp=NULL;}

    dim = src.getdim();
    comp = new double[dim];
    if(!comp){
      cout << "no memory for " << dim << "D vector..." << endl;
      exit(1);
    }
  }

  for(i = 0 ; i < dim ; i++) comp[i] = src.getvalue(i);

  return *this;
}

Vector & Vector::operator+=(const Vector & src)
{
  int  i;

  if(this->getdim() != src.getdim()) {
    if(comp) delete [] comp;

    dim = src.getdim();
    comp = new double[dim];

    for(i = 0 ; i < dim ; i++) comp[i] = ZeroValue;
  }

  this->add(src);

  return *this;
}

Vector &  Vector::showcontents(const char *front, const char *end)
{
  int      i;

  cout << front;

  cout << "[" << dim << "D]";

  cout << "<";
  cout.precision(SIGDIGITS);
  cout << scientific;
  for (i = 0 ; i < this->getdim() ; i++)
  cout << " " << this->getvalue(i) << " ";
  cout << ">" << end;

  return  *this;
}

Vector &  Vector::showcontents(const char *str)
{
  int      i;

  cout << str << "[" << dim << "D]";

  cout.precision(SIGDIGITS);
  cout << scientific;
  cout << "<";
  for (i = 0 ; i < this->getdim() ; i++)
  cout << " " << this->getvalue(i) << " ";
  cout << ">" << endl;

  return  *this;
}

Vector &  Vector::print(const char *end)
{
  int    i;

  for(i = 0 ; i < this->getdim() ; i++) printf("  %24.16e", this->getvalue(i));
  printf("%s", end);

  return  *this;
}

Vector &  Vector::print(FILE *fp, const char *end)
{
  int    i;

  for(i = 0 ; i < this->getdim() ; i++) fprintf(fp,"  %24.16e", this->getvalue(i));
  fprintf(fp,"%s", end);

  return  *this;
}

Vector &  Vector::print(ofstream & myout, const char *end)
{
  int    i;

  for(i = 0 ; i < this->getdim() ; i++) myout << "  " << this->getvalue(i);
  myout <<  end;

  return  *this;
}

Vector &  Vector::gnuplot(FILE *fp, const char *end)
{
  int    i;

  for(i = 0 ; i < this->getdim() ; i++) fprintf(fp,"  %24.16e", ZeroValue);
  for(i = 0 ; i < this->getdim() ; i++) fprintf(fp,"  %24.16e", this->getvalue(i));
  fprintf(fp,"%s", end);

  return  *this;
}

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                (CoLMat) Matrix of column vectors Class                 |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

CoLMaT::CoLMaT()
:col_dim(0),vec_n(0)
{
  colvec = NULL;
}

CoLMaT::CoLMaT(int dms)
:col_dim(dms),vec_n(dms)
{
  int  i;

  colvec = new Vector[dms];
  for(i = 0 ; i < dms ; i++) colvec[i] = Vector(dms);
}

CoLMaT::CoLMaT(int cdim, int rdim)
:col_dim(cdim),vec_n(rdim)
{
  int  i;

  colvec = new Vector[rdim];
  for(i = 0 ; i < rdim ; i++) colvec[i] = Vector(cdim);
}

CoLMaT::CoLMaT(const CoLMaT & src)
{
  int  i;

  col_dim = src.getcoldim();
  vec_n = src.getrowdim();

  if(vec_n == 0) colvec = NULL;
  else colvec = new Vector[vec_n];

  for(i = 0 ; i < vec_n ; i++) colvec[i] = src.getcolvector(i);
}

CoLMaT::~CoLMaT()
{
  if(colvec) delete [] colvec;
}

Vector &  CoLMaT::getcolvector(int k) const
{
  if(!colvec) {
    fprintf(stderr,"Stop at CoLMaT::getcolvector() [colvec is NULL]\n");
    exit(0);
  }

  return  colvec[k];
}

CoLMaT &  CoLMaT::setvalue(int i, int j, double value)
{
  if(!colvec) return *this;

  if(i < this->getcoldim() && j < this->getrowdim())
  colvec[j].setvalue(i, value);

  return *this;
}

CoLMaT &  CoLMaT::reset()
{
  int      i, j;

  for(j = 0 ; j < this->getrowdim() ; j++)
  for(i = 0 ; i < this->getcoldim() ; i++)
  this->setvalue(i, j, ZeroValue);

  return *this;
}

CoLMaT &  CoLMaT::setOnxn(int ndim)
{
  Vector   vec(ndim);

  vec.setzerovector(ndim);
  this->putvector(ndim-1, vec);
  this->reset();

  return *this;
}

CoLMaT &  CoLMaT::setInxn(int ndim)
{
  int      i;

  this->setOnxn(ndim);

  for(i = 0 ; i < ndim ; i++) colvec[i].setvalue(i, UnitValue);

  return *this;
}

CoLMaT &  CoLMaT::cut(int dim_new)
{
  int    i;

  for(i = 0 ; i < this->getrowdim() ; i++) colvec[i].cut(dim_new);

  col_dim = dim_new;

  return *this;
}

CoLMaT &  CoLMaT::putvector(int kth, Vector & vec)
{
  int      i, rowdim_org=this->getrowdim();
  CoLMaT  *hold=NULL;

  if(kth < 0) return *this;

  this->cut(vec.getdim());

  if(kth < rowdim_org) colvec[kth] = vec;
  else {
    hold = new CoLMaT(*this);

    if(colvec) delete [] colvec;

    colvec = new Vector[kth+1];

    for(i = 0 ; i < rowdim_org ; i++)
    colvec[i] = hold->getcolvector(i);

    colvec[kth] = vec;
  }

  vec_n = kth+1;

  delete hold;

  return *this;
}

double  CoLMaT::getvalue(int i, int j)
{
  double value;

  if(j < this->getrowdim())
  value = this->getcolvector(j).getvalue(i);
  else
  value = ZeroValue;

  return value;
}

CoLMaT &  CoLMaT::scalar(double alpha)
{
  int       i, j;
  double    value;

  for(j = 0 ; j < this->getrowdim() ; j++)
  for(i = 0 ; i < this->getcoldim() ; i++) {
    value = this->getvalue(i, j);
    colvec[j].setvalue(i, alpha*value);
  }

  return *this;
}

CoLMaT &  CoLMaT::rowmatrix(Vector & vec)
{
  int    vec_dim = vec.getdim();
  int    i;

  if(colvec) delete [] colvec;

  colvec = new Vector[vec_dim];

  for(i = 0 ; i < vec_dim ; i++)
  colvec[i].cut(1).setvalue(0, vec.getvalue(i));

  col_dim = 1;
  vec_n = vec_dim;

  return  *this;
}

CoLMaT &  CoLMaT::colmatrix(Vector & vec)
{
  int    vec_dim = vec.getdim();

  if(colvec) delete [] colvec;

  colvec = new Vector[1];

  colvec[0] = vec;

  col_dim = vec_dim;
  vec_n = 1;

  return  *this;
}

CoLMaT &  CoLMaT::rankone(Vector & u, Vector & v)
{
  CoLMaT  umat, vmat;

  umat.colmatrix(u);
  vmat.rowmatrix(v);

  *this = umat*vmat;

  return  *this;
}

CoLMaT &  CoLMaT::showsize(const char *str)
{
  fprintf(stderr,"%s[%d X %d]\n", str, col_dim, vec_n);

  return *this;
}

CoLMaT &  CoLMaT::showcontents(const char *str)
{
  int   i, j;

  fprintf(stderr,"(%d X %d)-size\n", this->getcoldim(), this->getrowdim());
  for(i = 0 ; i < this->getcoldim() ; i++){
    fprintf(stderr,"%s[ ", str);

    for(j = 0 ; j < this->getrowdim() ; j++)
    fprintf(stderr," %16.8e ", this->getvalue(i, j));

    fprintf(stderr," ]\n");
  }

  return *this;
}

CoLMaT & CoLMaT::add(CoLMaT & src)
{
  int  i, scdim = src.getcoldim();
  int  srdim = src.getrowdim();

  if(this->getcoldim() != scdim || this->getrowdim() != srdim) {
    fprintf(stderr,"size mismatch in CoLMaT::add.\n");
    exit(1);
  }

  for(i = 0 ; i < this->getcoldim() ; i++)
  this->getcolvector(i).add(src.getcolvector(i));

  return *this;
}

Vector & CoLMaT::operator[](int off)
{
  if(off >= 0 && off < this->getrowdim()) return this->getcolvector(off);
  else {
    fprintf(stderr,"no reference at [%d]column vector for %d column vectors.\n", off, this->getrowdim());
    exit(1);
  }
}

CoLMaT  CoLMaT::operator+(CoLMaT & src)
{
  CoLMaT   tmp(*this);

  tmp.add(src);

  return tmp;
}

Vector  CoLMaT::operator*(Vector & vec)
{
  int      i;
  int      m = this->getcoldim();
  int      n = this->getrowdim();
  Vector   hold(m), tmp(m);

  for(i = 0 ; i < n ; i++){
    tmp = this->getcolvector(i);
    hold += tmp.scalar(vec.getvalue(i));
  }

  return hold;
}

CoLMaT  CoLMaT::operator*(CoLMaT & src)
{
  int      mcol = this->getcoldim();
  int      nrow = src.getrowdim();
  int      i;
  CoLMaT   dest(mcol, nrow);

  for(i = 0 ; i < nrow ; i++)
  dest[i] = (*this)*src.getcolvector(i);

  return  dest;
}

CoLMaT & CoLMaT::operator=(const CoLMaT & src)
{
  int  i;

  if(this == &src) return *this;

  if(this->getrowdim() != src.getrowdim()) {
    if(colvec) {delete [] colvec; colvec = NULL;}

    col_dim = src.getcoldim();
    vec_n = src.getrowdim();

    colvec = new Vector[vec_n];
    if(!colvec){
      fprintf(stderr,"colvec[%d] fails at CoLMaT::operator=.\n", vec_n);
      exit(1);
    }
  }

  for(i = 0 ; i < vec_n ; i++) colvec[i] = src.getcolvector(i);

  return *this;
}

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|           (OPRT) Operators of projectors and Rotators Class            |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

OPRT::OPRT()
:dim(0)
{
  mat=NULL;
}

OPRT::OPRT(int dm)
:dim(dm)
{
  mat=new CoLMaT(dm);
}

OPRT::OPRT(const OPRT & src)
{
  int i, j;

  dim = src.getdim();

  mat = new CoLMaT(dim);

  for(i = 0 ; i < dim ; i++)
  for(j = 0 ; j < dim ; j++)
  this->setvalue(i, j, src.getvalue(i,j));
}

OPRT::~OPRT()
{
  if(mat) delete mat;
}

OPRT & OPRT::setvalue(int i, int j, double value)
{
  if(!this->getmatrix()) return *this;

  this->getmatrix()->setvalue(i, j, value);

  return *this;
}

OPRT & OPRT::setmatrix(int dm)
{
  if(!this->getmatrix())  {
    mat=new CoLMaT(dm);
    dim = dm;
  }

  return *this;
}

double OPRT::getvalue(int i, int j) const
{
  double   value;

  value = this->getmatrix()->getvalue(i, j);

  return value;
}

OPRT & OPRT::add(const OPRT & src)
{
  int       i, j;
  double    value;

  if(dim != src.getdim()) {
    fprintf(stderr,"matrix size mismatch at OPRT::add.\n");
    exit(1);
  }

  for(i = 0 ; i < dim ; i++)
  for(j = 0 ; j < dim ; j++){
    value = this->getvalue(i, j);
    this->setvalue(i, j, value+src.getvalue(i, j));
  }

  return   *this;
}

Vector  OPRT::acting(Vector & x)
{
  Vector   tmp_vec;
  CoLMaT   *mp;

  //    if(dim != x.getdim()) {
  //        fprintf(stderr,"size mismatch [Pxv]at OPRT::acting.\n");
  //        exit(1);
  //   }

  mp = this->getmatrix();

  tmp_vec = (*mp)*x;

  return  tmp_vec;
}

OPRT & OPRT::OP_onto(Vector & vec)
{
  int      i, j, qdm=vec.getdim();
  double   value, scaling=sqrt(vec.dot(vec));
  Vector   q(qdm);

  if(dim != vec.getdim()) {
    if(mat) delete mat;

    dim = qdm; // caution
    mat = new CoLMaT(qdm);
    if(!mat) {
      fprintf(stderr,"memory fails OPRT::OP_onto.\n");
      exit(1);
    }
  }

  q = vec; q.scalar(UnitValue/scaling);
  for(i = 0 ; i < dim ; i++)
  for(j = 0 ; j < dim ; j++){
    value = q.getvalue(i)*q.getvalue(j);
    this->setvalue(i, j, value);
  }

  return  *this;
}

OPRT & OPRT::COP_onto(Vector & vec)
{
  int      i;
  double   value;

  this->OP_onto(vec); // P onto <vec>

  this->getmatrix()->scalar(-UnitValue); // -P

  for(i = 0 ; i < dim ; i++){
    value = this->getvalue(i, i);
    this->setvalue(i, i, UnitValue+value);
  } // I - P onto <vec>^orth

  return  *this;
}

OPRT & OPRT::OP_onto(Vector & vec1, Vector & vec2)
{
  int        i, j;
  double     value;
  Vector     tmp_vec;
  OPRT  P, Q;

  if(vec1.getdim() != vec2.getdim()) {
    fprintf(stderr,"size mismatch [Pxv]at OPRT::OP_onto.\n");
    exit(1);
  }

  if(dim != vec1.getdim()) {
    if(mat) delete mat;

    dim = vec1.getdim(); // caution
    mat = new CoLMaT(dim);
    if(!mat) {
      fprintf(stderr,"memory fails OPRT::OP_onto.\n");
      exit(1);
    }
  }

  P.OP_onto(vec1);

  Q.COP_onto(vec1);
  tmp_vec = Q.acting(vec2);

  Q.OP_onto(tmp_vec);

  for(i = 0 ; i < dim ; i++)
  for(j = 0 ; j < dim ; j++){
    value = P.getvalue(i, j) + Q.getvalue(i, j);
    this->setvalue(i, j, value);
  }

  return  *this;
}

OPRT & OPRT::OP_onto(int nvec, Vector *vecs)
{
  int        i, j;
  int        cdim, this_dim;
  double     value;
  Vector    *qvec;
  CoLMaT     tmp_mat;
  OPRT  ImP;

  if(nvec == 0 || vecs == NULL) {
    this_dim = this->getdim();
    this->setmatrix(this_dim).getmatrix()->setOnxn(this_dim);
    return   *this;
  }

  cdim = vecs[0].getdim();
  ImP.setmatrix(cdim).getmatrix()->setOnxn(cdim); // Initializing ImP
  this->setmatrix(cdim).getmatrix()->setOnxn(cdim);  // Initializing this matrix

  qvec = new Vector[nvec];
  if(!qvec) {
    fprintf(stderr,"No memory for vectors in OP_onto()...\n");
    exit(1);
  }

  for(i = 0 ; i < nvec ; i++) qvec[i] = vecs[i];

  for(i = 1 ; i < nvec+1 ; i++) {
    value = qvec[i-1].norm();
    if(value < NearZero) {
      fprintf(stderr,"not LI vectors...\n");
      return *this;
    }
    qvec[i-1].scalar(UnitValue/value);
    tmp_mat.rankone(qvec[i-1], qvec[i-1]);
    this->getmatrix()->add(tmp_mat);

    ImP.COP_onto(qvec[i-1]);
    for(j = i ; j < nvec ; j++) qvec[j] = ImP.acting(qvec[j]);
  }

  delete [] qvec;

  return  *this;
}

OPRT &  OPRT::Omega(Vector & omg)
{
  int      dim3 = omg.getdim();
  double   value;
  Vector   tmp(omg);

  if(dim3 != D3) {
    fprintf(stderr,"It works only in 3-dimension.\n");
    exit(0);
  }

  value = omg.norm();
  tmp.scalar(UnitValue/value);

  this->setmatrix(dim3).getmatrix()->setOnxn(dim3);

  this->getmatrix()->getcolvector(1).setvalue(0, -tmp[2]);
  this->getmatrix()->getcolvector(2).setvalue(0,  tmp[1]);
  this->getmatrix()->getcolvector(2).setvalue(1, -tmp[0]);

  this->getmatrix()->getcolvector(0).setvalue(1,  tmp[2]);
  this->getmatrix()->getcolvector(0).setvalue(2, -tmp[1]);
  this->getmatrix()->getcolvector(1).setvalue(2,  tmp[0]);

  return *this;
}

OPRT &  OPRT::showcontents(const char *str)
{
  this->getmatrix()->showcontents(str);

  return  *this;
}

Vector  OPRT::operator*(Vector & post)
{
  return   this->acting(post);
}

OPRT & OPRT::operator=(const OPRT & src)
{
  int  i, j;

  if(this == &src) return *this;

  if(this->getdim() != src.getdim()){
    dim = src.getdim();
    if(mat) delete mat;

    mat = new CoLMaT(dim);
  }

  for(i = 0 ; i < dim ; i++)
  for(j = 0 ; j < dim ; j++)
  this->setvalue(i, j, src.getvalue(i,j));

  return *this;
}

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                     (ROT) Rotators(3D only) Class                      |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

Rot::Rot()
:omega(D3),Pw(D3),ImPw(D3),Omg(D3)
{
}

Rot::Rot(Vector & w)
:omega(w),Pw(D3),ImPw(D3),Omg(D3)
{
  Pw.OP_onto(w);
  ImPw.COP_onto(w);
  Omg.Omega(w);
}

Rot::Rot(const Rot & src)
:omega(src.getRotAxis()),Pw(src.getPw()),ImPw(src.getImPw()),Omg(src.getOmg())
{
}

Rot::~Rot() {}

Rot & Rot::setRotAxis(Vector & vec)
{
  omega = vec;

  return *this;
}

Vector   Rot::acting(Vector & vec, double theta)
{
  Vector  res(D3), tmp(D3);

  res = this->getOPRT(0).acting(vec);

  tmp = this->getOPRT(1).acting(vec);
  tmp.scalar(cos(theta));
  res.add(tmp);

  tmp = this->getOPRT(2).acting(vec);
  tmp.scalar(sin(theta));
  res.add(tmp);

  return res;
}

Rot & Rot::setOPRT(Vector & w)
{
  omega = w;

  Pw.OP_onto(w);
  ImPw.COP_onto(w);
  Omg.Omega(w);

  return *this;
}

OPRT  Rot::getOPRT(int off) const
{
  switch(off) {
    case 0: return Pw;
    case 1: return ImPw;
    case 2: return Omg;
    default: fprintf(stderr,"Wrong offset in Rot:getOPRT...\n");
    exit(0);
  }
}

Rot & Rot::showcontents(const char *str )
{
  fprintf(stderr,"------(Rotation %s)------\n", str);
  omega.showcontents("rot-axis");
  getPw().showcontents("Pw");
  getImPw().showcontents("ImPw");
  getOmg().showcontents("Omg");

  return *this;
}

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                 (Cell) Unit cell of linked list Class                  |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

Cell::Cell()
:index(-1)
{
  nextcell = NULL;
}

Cell::Cell(int dim)
:index(-1),vec(dim)
{
  nextcell = NULL;
}

Cell::Cell(Vector & src)
:index(-1)
{
  vec = src;
  nextcell = NULL;
}

Cell::Cell(Vector & src, char bdc, double bdv, char st)
:index(-1)
{
  vec = src;
  nextcell = NULL;
  cond = bdc;
  value = bdv;
  status = st;
}

Cell::Cell(Vector & src, char bdc, double bdv, char st, Vector & normalvector)
:index(-1)
{
  vec = src;
  nextcell = NULL;
  cond = bdc;
  value = bdv;
  status = st;
  normal = normalvector;
}

Cell::Cell(int offs, Vector & src)
:index(offs)
{
  vec = src;
  nextcell = NULL;
}

Cell::Cell(const Cell & src)
{
  index = src.Index();
  vec = src.Item();
  nextcell = src.Next();
}

Cell::~Cell()
{}

  Cell & Cell::SetItem(Vector & src)
  {
    index = -1;

    vec = src;

    return  *this;
  }

  Cell & Cell::SetItem(int offs, Vector & src)
  {
    index = offs;

    vec = src;

    return  *this;
  }

  Cell & Cell::SetNormal(Vector & src)
  {
    index = -1;

    normal = src;

    return  *this;
  }

  Cell & Cell::SetNormal(int offs, Vector & src)
  {
    index = offs;

    normal = src;

    return  *this;
  }

  Cell & Cell::SetNext(Cell* src)
  {
    nextcell = src;

    return  *this;
  }

  Cell & Cell::SetCond(char bdc)
  {
    cond = bdc;

    return  *this;
  }

  Cell & Cell::SetValue(double bdv)
  {
    value = bdv;

    return  *this;
  }

  Cell & Cell::SetStatus(char st)
  {
    status = st;

    return  *this;
  }

  Cell & Cell::SetAxialLine (AxialLine *xsrc, AxialLine *ysrc) {

    xaxialline = xsrc;
    yaxialline = ysrc;

    return *this;
  }

  Cell & Cell::AppendCell(int nVD)
  {
    Cell *pcl;

    if(!this->Next()) {
      nextcell = new Cell(nVD);
      if(!CellMemoryCheck(nextcell)) exit(0);
    } else {
      pcl = this->Next();
      nextcell = new Cell(nVD);
      if(!CellMemoryCheck(nextcell)) exit(0);
      nextcell->SetNext(pcl);
    }

    return  *this;
  }

  Cell & Cell::AppendCell(Vector & src)
  {
    Cell *pcl;

    if(!this->Next()) {
      nextcell = new Cell;
      if(!CellMemoryCheck(nextcell)) exit(0);
      nextcell->SetItem(src);
    } else {
      pcl = this->Next();
      nextcell = new Cell;
      if(!CellMemoryCheck(nextcell)) exit(0);
      nextcell->SetItem(src);
      nextcell->SetNext(pcl);
    }

    return  *this;
  }

  Cell & Cell::AppendCell(Vector & src, char bdc, double bdv, char st)
  {
    Cell *pcl;

    if(!this->Next()) {
      nextcell = new Cell;
      if(!CellMemoryCheck(nextcell)) exit(0);
      nextcell->SetItem(src);
      nextcell->SetCond(bdc);
      nextcell->SetValue(bdv);
      nextcell->SetStatus(st);
    } else {
      pcl = this->Next();
      nextcell = new Cell;
      if(!CellMemoryCheck(nextcell)) exit(0);
      nextcell->SetItem(src);
      nextcell->SetNext(pcl);
      nextcell->SetCond(bdc);
      nextcell->SetValue(bdv);
      nextcell->SetStatus(st);
    }

    return  *this;
  }

  Cell & Cell::AppendCell(Vector & src, char bdc, double bdv, char st, Vector & normal)
  {
    Cell *pcl;

    if(!this->Next()) {
      nextcell = new Cell;
      if(!CellMemoryCheck(nextcell)) exit(0);
      nextcell->SetItem(src);
      nextcell->SetCond(bdc);
      nextcell->SetValue(bdv);
      nextcell->SetStatus(st);
      nextcell->SetNormal(normal);
    } else {
      pcl = this->Next();
      nextcell = new Cell;
      if(!CellMemoryCheck(nextcell)) exit(0);
      nextcell->SetItem(src);
      nextcell->SetNext(pcl);
      nextcell->SetCond(bdc);
      nextcell->SetValue(bdv);
      nextcell->SetStatus(st);
      nextcell->SetNormal(normal);
    }

    return  *this;
  }

  Cell & Cell::AppendCell(int offs, Vector & src)
  {
    this->AppendCell(src);
    this->Next()->SetIndex(offs);

    return  *this;
  }

  Cell & Cell::AppendCell(int offs, Vector & src, char bdc, double bdv, char st)
  {
    this->AppendCell(src, bdc, bdv, st);
    this->Next()->SetIndex(offs);

    return  *this;
  }

  Cell & Cell::AppendCell(int offs, Vector & src, char bdc, double bdv, char st, Vector & normal)
  {
    this->AppendCell(src, bdc, bdv, st, normal);
    this->Next()->SetIndex(offs);

    return  *this;
  }

  Cell & Cell::showcontents(const char *comnt)
  {
    char wording[1024];

    sprintf(wording,"%s(%d)", comnt, this->Index());
    printf("%s%c\t%s%f\t%s%c\n", "contition = ", this->cond, "value = ", this->value, "status = ", this->status);
    this->Item().showcontents(wording);

    return  *this;
  }

  /****************************************************************************/
  /*+------------------------------------------------------------------------+*/
  /*|               (Thing) Master of Cells(Linked List) Class               |*/
  /*+------------------------------------------------------------------------+*/
  /****************************************************************************/

  Thing::Thing()
  :ncell(0)
  {
    headcell = NULL;
    run = headcell;
  }

  Thing::Thing(const Thing & src)
  {
    Cell *srccell, *destcell;
    Cell *pc;
    Vector  vtmp;

    ncell = src.Howmany();

    if(!src.Head()) {
      ncell = 0;
      run = headcell = NULL;
    } else {
      // desthead
      vtmp = src.Head()->Item();
      headcell = new Cell(vtmp);
      if(!CellMemoryCheck(headcell)) exit(0);
      headcell->SetIndex(src.Head()->Index());

      // pointing headcells
      srccell = src.Head();
      destcell = this->Head();

      // copy  cells until null
      while(srccell->Next()) {
        vtmp = srccell->Next()->Item();
        pc = new Cell(vtmp);
        if(!CellMemoryCheck(pc)) exit(0);
        pc->SetIndex(srccell->Next()->Index());
        destcell->SetNext(pc);

        srccell = srccell->Next();
        destcell = destcell->Next();
      }
    }
  }

  Thing::~Thing()
  {
    this->GoTail();

    while(this->Run()) this->RemoveCell();
  }

  Thing & Thing::MakeCell(int nVD)
  {
    if(!headcell) {
      headcell = new Cell(nVD);
      if(!CellMemoryCheck(headcell)) exit(0);
      run = headcell;

      ncell++;

      return  *this;
    }

    if(!run) return *this;

    run->AppendCell(nVD); // append a cell at [run]->next
    run = run->Next();

    ncell++;

    return *this;
  }

  Thing & Thing::MakeCell(Vector & src)
  {
    if(!headcell) {
      headcell = new Cell(src);
      if(!CellMemoryCheck(headcell)) exit(0);
      run = headcell;

      ncell++;

      return  *this;
    }

    if(!run) return *this;

    run->AppendCell(src); // Append a cell at [run]->next
    run = run->Next();

    ncell++;

    return *this;
  }

  Thing & Thing::MakeCell(Vector & src, char bdc, double bdv, char st)
  {
    if(!headcell) {
      headcell = new Cell(src, bdc, bdv, st);
      if(!CellMemoryCheck(headcell)) exit(0);
      run = headcell;

      ncell++;

      return  *this;
    }

    if(!run) return *this;

    run->AppendCell(src, bdc, bdv, st); // Append a cell at [run]->next
    run = run->Next();

    ncell++;

    return *this;
  }

  Thing & Thing::MakeCell(Vector & src, char bdc, double bdv, char st, Vector & normal)
  {
    if(!headcell) {
      headcell = new Cell(src, bdc, bdv, st, normal);
      if(!CellMemoryCheck(headcell)) exit(0);
      run = headcell;

      ncell++;

      return  *this;
    }

    if(!run) return *this;

    run->AppendCell(src, bdc, bdv, st, normal); // Append a cell at [run]->next
    run = run->Next();

    ncell++;

    return *this;
  }

  Thing & Thing::MakeCell(int offs, Vector & src)
  {
    if(!headcell) {
      headcell = new Cell(src);
      if(!CellMemoryCheck(headcell)) exit(0);
      headcell->SetIndex(offs);
      run = headcell;

      ncell++;

      return  *this;
    }

    if(!run) return *this;

    run->AppendCell(offs, src); // Append a cell at [run]->next
    run = run->Next();

    ncell++;

    return *this;
  }

  Thing & Thing::MakeCell(int offs, Vector & src, char bdc, double bdv, char st)
  {
    if(!headcell) {
      headcell = new Cell(src, bdc, bdv, st);
      if(!CellMemoryCheck(headcell)) exit(0);
      headcell->SetIndex(offs);
      run = headcell;

      ncell++;

      return  *this;
    }

    if(!run) return *this;

    run->AppendCell(offs, src, bdc, bdv, st); // Append a cell at [run]->next
    run = run->Next();

    ncell++;

    return *this;
  }

  Thing & Thing::MakeCell(int offs, Vector & src, char bdc, double bdv, char st, Vector & normal)
  {
    if(!headcell) {
      headcell = new Cell(src, bdc, bdv, st, normal);
      if(!CellMemoryCheck(headcell)) exit(0);
      headcell->SetIndex(offs);
      run = headcell;

      ncell++;

      return  *this;
    }

    if(!run) return *this;

    run->AppendCell(offs, src, bdc, bdv, st, normal); // Append a cell at [run]->next
    run = run->Next();

    ncell++;

    return *this;
  }

  Thing & Thing::RemoveCell()
  {
    Cell *pc=headcell;
    Cell *chold;

    if(!run) return *this;

    // find run-position
    while(pc)
    if(pc->Next() == run) break;
    else  pc = pc->Next();

    // erase run-cell
    chold = run->Next();
    run->SetNext(NULL);
    delete  run;

    // hold the next-cell
    if(pc) pc->SetNext(chold);
    else   pc = headcell = chold;

    run = pc;

    ncell--;

    return *this;
  }

  Thing & Thing::RemoveAllCell()
  {
    this->GoHead();

    while(this->Run()) this->RemoveCell();

    return *this;
  }

  Thing & Thing::TakeoffHeadCell () {

    ncell = 0;
    headcell = NULL;
    run = headcell;

    return *this;

  }

  Thing & Thing::GoPost()
  {
    if(!run) return *this;

    run = run->Next();

    return  *this;
  }

  Thing & Thing::GoPost(int nP)
  {
    int  i;

    for(i = 0 ; i < nP ; i++) this->GoPost();

    return  *this;
  }

  Thing & Thing::GoPre()
  {
    Cell *pc=headcell;

    if(!run) return *this;

    if(run == headcell) return *this;

    while(pc) {
      if(pc->Next() != run) pc = pc->Next();
      else                    break;
    }

    run = pc;

    return  *this;
  }

  Thing & Thing::GoPre(int ntimes)
  {
    int  i;

    for(i = 0 ; i < ntimes ; i++) this->GoPre();

    return  *this;
  }

  Thing & Thing::GoTail()
  {
    if(!this->Head()) {run = NULL; return *this;}

    run = this->Head();
    while(run->Next()) this->GoPost();

    return  *this;
  }

  Thing & Thing::GoHead()
  {
    if(!this->Head()) {run = NULL; return *this;}

    run = this->Head();

    return  *this;
  }

  Thing & Thing::MoveRun(int offs)
  {
    run = this->WhichCell(offs);

    return  *this;
  }

  Thing & Thing::MoveRun(Cell *src)
  {
    run = src;

    return  *this;
  }

  Thing & Thing::SetHead(Cell *src)
  {
    headcell = src;

    return  *this;
  }

  Thing & Thing::HeadChange()
  {
    Cell  *pc;

    if(!this->Run()) return *this; // not chaned

    pc = this->Run(); // save run position

    this->GoTail(); // move run to tail
    if(!this->Run()) return *this;

    this->Run()->SetNext(headcell); // make it closed

    headcell = pc; // move head to the run position

    while(pc->Next() != headcell) pc = pc->Next();
    pc->SetNext(NULL); // make it open

    run = headcell;

    return  *this;
  }

  Thing & Thing::HeadChange(int offs)
  {
    run = this->WhichCell(offs);
    this->HeadChange();

    return  *this;
  }

  Thing & Thing::Swapping()
  {
    Cell *pcl=NULL;

    // do nothing if run == NULL  or run->next == NULL
    if(!this->Run()) return *this;
    if(!this->Run()->Next()) return *this;

    pcl = this->Run(); // save run position

    if(pcl != headcell) {
      this->GoPre(); // move run up
      this->Run()->SetNext(pcl->Next());
      pcl->SetNext(this->Run()->Next()->Next());
      this->Run()->Next()->SetNext(pcl);
    } else {
      this->SetHead(pcl->Next());
      pcl->SetNext(this->Head()->Next());
      this->Head()->SetNext(pcl);
    }

    this->MoveRun(pcl); // recover the run position

    return  *this;
  }

  Thing & Thing::Swapping(Cell * src)
  {
    Cell *pcl=NULL;

    if(!src) return *this;

    pcl = this->Run(); // save run position

    this->MoveRun(src);
    this->Swapping(); // swapping src and its next

    this->MoveRun(pcl); // recover the run position

    return  *this;
  }

  int  Thing::Ordering(char *object, int offs)
  // offs = the offset of [Item]
  {
    Cell *pc;

    if(!this->OnRun()) return 0;
    if(!this->Run()->Next()) return 0;

    pc = this->Run();
    if(!strcmp(object, "Index")) {
      if(pc->Index() > pc->Next()->Index()) return -1;
      else                                    return  1;
    } else if(!strcmp(object, "Item")) {
      if(pc->Item(offs) > pc->Next()->Item(offs)) return -1;
      else                                              return  1;
    } else      return 0;
  }

  Thing & Thing::Sorting(char *cmpHD, char pm, int offs)
  // Exchanging adjacent items if pm is violated
  // [*cmpHD] = [Index] or [Item]
  // offs = the offset of [Item]
  // pm = +(ascending) or -(descending)
  {
    int    ord;

    switch(pm) {
      case '+':
      ord = this->Ordering(cmpHD, offs);
      if(ord == -1) this->Swapping();
      else            this->GoPost();

      break;
      case '-':
      ord = this->Ordering(cmpHD, offs);
      if(ord == 1)  this->Swapping();
      else            this->GoPost();

      break;
      default: break;
    }

    return *this;
  }

  Thing &  Thing::Sorting(const char *cmp)
  // [*cmp]= comparing object of format [X...X{n..n}+(-)]
  // X...X = [Index] or [Item]
  // n...n = the offset of Item which may not be given
  // +(ascending) or -(descending)
  {
    const char  *pcmp;
    char  *pHD, *pNB;
    char   cc, cmpNB[64], cmpHD[64];
    int    offs, count, downto;

    // parsing [cmp]
    pcmp = cmp; pHD = cmpHD; pNB = cmpNB;
    while(*pcmp != '\0') {

      if(IsIn(*pcmp, "+-")) cc = *pcmp;
      else if(IsIn(*pcmp, "0123456789")) *pNB++ = *pcmp;
      else                                   *pHD++ = *pcmp;

      pcmp++;
    }
    *pNB = *pHD = '\0';
    offs = atoi(cmpNB); // offset

    if(!IsIn(cc, "+-")) {
      cout << "[cmp]-string = [...+(-)] ended with + or -" << endl;
      exit(0);
    }

    // Sorting Process
    if(!this->OnRun()) return *this;

    // counting down to the end by initial sorting
    count = 0;
    while(this->Run()->Next()) {this->Sorting(cmpHD, cc, offs); count++;}

    // sorting by swapping successively
    while(count) {
      this->GoPre(count--);

      downto = 0;
      while(downto++ < count) this->Sorting(cmpHD, cc, offs);
    }

    return *this;
  }

  Thing &  Thing::Sorting(const char *cmp, int cell_N)
  // [*cmp]= comparing object of format [X...X{n..n}+(-)]
  // cell_N = the number of cells down from [run] to sort
  // Cell[(run)]->Cell[+1]->Cell[+2]->...->Cell[+cell_N]
  {
    const char  *pcmp;
    char  *pHD, *pNB;
    char   cc, cmpNB[64], cmpHD[64];
    int    offs, count, downto;

    // parsing [cmp]
    pcmp = cmp; pHD = cmpHD; pNB = cmpNB;
    while(*pcmp != '\0') {

      if(IsIn(*pcmp, "+-")) cc = *pcmp;
      else if(IsIn(*pcmp, "0123456789")) *pNB++ = *pcmp;
      else                               *pHD++ = *pcmp;

      pcmp++;
    }
    *pNB = *pHD = '\0';
    offs = atoi(cmpNB); // offset

    if(!IsIn(cc, "+-")) {
      cout << "[cmp]-string = [...+(-)] ended with + or -" << endl;
      exit(0);
    }

    // Sorting Process
    if(!this->OnRun()) return *this;

    // the number of down cells to sort
    count = cell_N;

    // sorting successively
    while(count) {
      downto = 0;
      while(downto++ < count && this->Run()->Next() != NULL )
      this->Sorting(cmpHD, cc, offs);

      count = downto-1; // important!!!

      this->GoPre(count--);
    }

    return *this;
  }

  Thing &  Thing::Import(Thing & src, char pre_post)
  //  Putting  [src] into [*this] at [this-run] position
  //  i.e.  [src]-cells move into [this-thing]
  //  pre_post = '+' (place them rear of [run])
  //             '-' (place them front of [run])
  {
    Cell   *pthis;
    int     n_src, m_remain;

    // src->run is empty of cells
    if(!src.Run()) return *this;

    // this is not empty of cells but this->run is null
    if(this->Run() == NULL && this->Head() != NULL) return *this;

    // putting process

    n_src = src.Howmany();


    if(!this->Head()) {
      pthis = NULL;
      this->SetHead(src.Run());
    } else if(this->Run() == this->Head() && pre_post == '-') {
      pthis = this->Head();
      this->SetHead(src.Run());
    } else {
      if(pre_post == '-') this->GoPre();
      pthis = this->Run()->Next();
      this->Run()->SetNext(src.Run());
    }

    if(src.Run() != src.Head()) {
      src.GoPre();
      src.Run()->SetNext(NULL);
    } else   src.SetHead(NULL);

    this->GoTail();
    this->Run()->SetNext(pthis);

    m_remain = src.CountCells(); // how many cells remain in [src]
    src.SetNcell(m_remain);

    ncell += n_src - m_remain; // how many moving cells

    return  *this;
  }

  Thing &  Thing::Import(Thing & src, int nPost, char pre_post)
  {
    int     ntail;
    Thing   tmpTH;
    Cell   *prun;

    if(!src.Run()) return *this;

    if(nPost >= src.Howmany()) {
      this->Import(src, pre_post);
      return *this;
    }

    prun = src.Run(); // Save [src]-run position

    // cutting after [src]-run + nPost
    src.GoPost(nPost);
    tmpTH.Import(src, '+');
    ntail = tmpTH.Howmany();

    // Recover [src]-run
    src.MoveRun(prun);

    // Import [src]-section
    this->Import(src, pre_post);

    // Address cut [src]
    tmpTH.GoHead();
    src.Import(tmpTH, '+');

    // Move [src]-run to the cut position of [src]
    src.GoPre(ntail-1);

    return *this;
  }

  int   Thing::OnRun()
  {
    if(run)  return  1;
    else       return  0;
  }

  int   Thing::Numbering()
  {
    Cell *pc;
    int   count;

    pc = this->Head(); count = 0;
    while(pc) {
      pc->SetIndex(count++);
      pc = pc->Next();
    }

    return  count;
  }

  int   Thing::Numbering(int num)
  {
    Cell *pc;
    int   count;

    pc = this->Head(); count = num;
    while(pc) {
      pc->SetIndex(count++);
      pc = pc->Next();
    }

    return  count;
  }

  int   Thing::CountCells()
  {
    Cell *pc;
    int   count;

    pc = this->Head(); count = 0;
    while(pc) {count++; pc = pc->Next();}

    return  count;
  }

  Cell *  Thing::WhichCell(int idx)
  {
    Cell *pc;

    pc = this->Head();
    while(pc) {
      if(pc->Index() == idx) return pc;
      pc = pc->Next();
    }

    return  pc;
  }

  Thing & Thing::operator=(Thing & src)
  {
    Cell *srccell, *destcell;
    Cell *pc;
    Vector  vtmp;

    if(this == &src)  return *this;

    // remove all cells of this object
    this->RemoveAllCell();

    if(!src.Head()) {
      ncell = 0;
      run = headcell = NULL;
      return *this;
    }

    // desthead
    vtmp = src.Head()->Item();
    headcell = new Cell(vtmp);
    if(!CellMemoryCheck(headcell)) exit(0);
    headcell->SetIndex(src.Head()->Index());

    // pointing headcells
    srccell = src.Head();
    destcell = this->Head();

    // copy  cells until null
    while(srccell->Next()) {
      vtmp = srccell->Next()->Item();
      pc = new Cell(vtmp);
      if(!CellMemoryCheck(pc)) exit(0);

      pc->SetIndex(srccell->Next()->Index());
      destcell->SetNext(pc);

      srccell = srccell->Next();
      destcell = destcell->Next();
    }

    // how many cells
    ncell = src.Howmany();

    //  run position
    run = src.Run();

    return *this;
  }

  Cell & Thing::operator[](int offs)
  {
    Cell  *pc=NULL;

    pc = this->WhichCell(offs);
    if(!pc) {
      cout << "no cell to the offset [" << offs << "]..." << endl;
      exit(0);
    }

    return *pc;
  }

  Thing &  Thing::showcontents(const char *comnt)
  {
    Cell *pc;

    pc = headcell;

    cout << comnt << "[[" << this->Howmany() << "]]--->" << endl;
    while(pc){
      if(pc == run) pc->showcontents("*cell");
      else          pc->showcontents(" cell");
      pc = pc->Next();
    }

    return *this;
  }

  #endif
