#ifndef CLASS
#define CLASS

#include "Stuff++.hpp"

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                             Section Class                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

Section::Section () {

  name = new char;
  nE = 0;
  SL2D = NULL;
  SL3D = NULL;

}

Section::Section (string nm) {

  name = nm;
  nE = 0;
  SL2D = NULL;
  SL3D = NULL;

}

Section::Section (string nm, int n) {

  name = nm;
  nE = n;
  SL2D = NULL;
  SL3D = NULL;

}

Section::Section (string nm, int n, Line2D *ln2d) {

  name = nm;
  nE = n;
  SL2D = ln2d;
  SL3D = NULL;

}

Section::Section (string nm, int n, Triangle *tr3d) {

  name = nm;
  nE = n;
  SL2D = NULL;
  SL3D = tr3d;

}


Section::~Section () {

  delete [] SL2D;
  delete [] SL3D;

}

Section & Section::setSection (string nm, int n) {

  name = nm;
  nE = n;
  SL2D = NULL;
  SL3D = NULL;

  return *this;

}

Section & Section::setSection (string nm, int n, Line2D *ln2d) {

  name = nm;
  nE = n;
  SL2D = ln2d;
  SL3D = NULL;

  return *this;

}

Section & Section::setSection (string nm, int n, Triangle *tr3d) {

  name = nm;
  nE = n;
  SL2D = NULL;
  SL3D = tr3d;

  return *this;

}

Section & Section::showcontents (const char *comnt) {

  printf("%s\n", comnt);

  printf("%s%s%s\n", "SECTION", " ", this->getName().c_str());

  for (size_t i = 0; i < this->getEltnum(); i++) {

    this->get2DElt(i)->showcontents("");

  }

  return *this;

}

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                             Region Class                               |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

Region::Region () {

  nS = 0;
  sctn = NULL;
  sign = NULL;
  xgrid2d = NULL;
  ygrid2d = NULL;
  xgrid3d = NULL;
  ygrid3d = NULL;
  zgrid3d = NULL;
  nXgrid = 0,
  nYgrid = 0;
  nZgrid = 0;

}

Region::Region (string nm) {

  name = nm;
  nS = 0;
  sctn = NULL;
  sign = NULL;
  xgrid2d = NULL;
  ygrid2d = NULL;
  xgrid3d = NULL;
  ygrid3d = NULL;
  zgrid3d = NULL;
  nXgrid = 0,
  nYgrid = 0;
  nZgrid = 0;

}

Region::Region (string nm, int n) {

  name = nm;
  nS = n;
  sctn = NULL;
  sign = NULL;
  xgrid2d = NULL;
  ygrid2d = NULL;
  xgrid3d = NULL;
  ygrid3d = NULL;
  zgrid3d = NULL;
  nXgrid = 0,
  nYgrid = 0;
  nZgrid = 0;

}

Region::Region (string nm, int n, Vector *src) {

  name = nm;
  nS = n;
  gridinfo = *src;
  sctn = NULL;
  sign = NULL;
  xgrid2d = NULL;
  ygrid2d = NULL;
  xgrid3d = NULL;
  ygrid3d = NULL;
  zgrid3d = NULL;
  nXgrid = 0,
  nYgrid = 0;
  nZgrid = 0;

}

Region::~Region () {

  delete [] sctn;
  delete [] sign;
  delete   &gridinfo;
  delete [] xgrid2d;
  delete [] ygrid2d;
  delete [] xgrid3d;
  delete [] ygrid3d;
  delete [] zgrid3d;

}

Region & Region::setSection (Section *sct) {

  *sctn = sct;
  return *this;

}

Region & Region::setSection (Section *sct, int n) {

  if (!sctn) sctn = new Section*[this->getSectionnum()];

  sctn[n] = sct;
  return *this;

}

Region & Region::setSign (string sgn, int n) {

  if (!sign) sign = new double[this->getSectionnum()];

  if (!sgn.compare("+")) sign[n] =  UnitValue;
  if (!sgn.compare("-")) sign[n] = -UnitValue;

  return *this;

}

Region & Region::make2Dgrid () {

  nXgrid = 0, nYgrid = 0;

  for (double i = gridinfo[2] + gridinfo[6]; i < gridinfo[4] - 0.5 * gridinfo[6]; i += gridinfo[6]) nXgrid++;
  for (double i = gridinfo[1] + gridinfo[5]; i < gridinfo[3] - 0.5 * gridinfo[5]; i += gridinfo[5]) nYgrid++;

  xgrid2d = new Line2D[nXgrid];
  ygrid2d = new Line2D[nYgrid];

  nXgrid = 0, nYgrid = 0;

  for (double i = gridinfo[2] + gridinfo[6]; i < gridinfo[4] - 0.5 * gridinfo[6]; i += gridinfo[6]) {

    xgrid2d[nXgrid].SetLine(round(gridinfo[1] - 0.5 * gridinfo[6], Margin), round(i, Margin), round(gridinfo[3] + 0.5 * gridinfo[6], Margin), round(i, Margin)); nXgrid++;

  }

  for (double i = gridinfo[1] + gridinfo[5]; i < gridinfo[3] - 0.5 * gridinfo[5]; i += gridinfo[5]) {

    ygrid2d[nYgrid].SetLine(round(i, Margin), round(gridinfo[2] - 0.5 * gridinfo[5], Margin), round(i, Margin), round(gridinfo[4] + 0.5 * gridinfo[5], Margin)); nYgrid++;

  }

  return *this;

}

Region & Region::make3Dgrid () {

  // nXgrid = 0, nYgrid = 0, nZgrid = 0;
  //
  // for (double i = gridinfo[2] + gridinfo[6]; i < gridinfo[4] - 0.5 * gridinfo[6]; i += gridinfo[6]) nXgrid++;
  // for (double i = gridinfo[1] + gridinfo[5]; i < gridinfo[3] - 0.5 * gridinfo[5]; i += gridinfo[5]) nYgrid++;
  //
  // xgrid2d = new Line2D[nXgrid];
  // ygrid2d = new Line2D[nYgrid];
  //
  // nXgrid = 0, nYgrid = 0;
  //
  // for (double i = gridinfo[2] + gridinfo[6]; i < gridinfo[4] - 0.5 * gridinfo[6]; i += gridinfo[6]) {
  //
  //   xgrid2d[nXgrid].SetLine(gridinfo[1] - 0.5 * gridinfo[6], i, gridinfo[3] + 0.5 * gridinfo[6], i); nXgrid++;
  //
  // }
  //
  // for (double i = gridinfo[1] + gridinfo[5]; i < gridinfo[3] - 0.5 * gridinfo[5]; i += gridinfo[5]) {
  //
  //   ygrid2d[nYgrid].SetLine(i, gridinfo[2] - 0.5 * gridinfo[5], i, gridinfo[4] - 0.5 * gridinfo[5]); nYgrid++;
  //
  // }

  return *this;

}

Region & Region::setRegion (string nm, int n) {

  name = nm;
  nS = n;
  sctn = NULL;
  sign = NULL;

  return *this;

}

Region & Region::setRegion (string nm, int n, Vector *src) {

  name = nm;
  nS = n;
  gridinfo = *src;
  sctn = NULL;
  sign = NULL;

  return *this;

}

Region & Region::showcontents (const char *comnt) {

  printf("%s\n", comnt);

  printf("%s%s%s\n", "REGION", " ", this->getName().c_str());

  printf("%s\n", "Grid information ");
  gridinfo.showcontents("");

  // for (size_t i = 0; i < this->getSectionnum(); i++) printf("%f\t%s\n", this->getSign(i), this->getSection(i)->getName());
  //
  // if (this->getZgridnum() == 0) {
  //
  //   for (size_t i = 0; i < this->getXgridnum(); i++) {printf("%lu%s", i, " x-grid-line = "); get2Dxgrid(i)->showcontents("");}
  //   for (size_t i = 0; i < this->getYgridnum(); i++) {printf("%lu%s", i, " y-grid-line = "); get2Dygrid(i)->showcontents("");}
  //
  // } else {
  //
  //   for (size_t i = 0; i < this->getXgridnum(); i++) {printf("%lu%s", i, " x-grid-line = "); get3Dxgrid(i)->showcontents("");}
  //   for (size_t i = 0; i < this->getYgridnum(); i++) {printf("%lu%s", i, " y-grid-line = "); get3Dygrid(i)->showcontents("");}
  //   for (size_t i = 0; i < this->getZgridnum(); i++) {printf("%lu%s", i, " z-grid-line = "); get3Dzgrid(i)->showcontents("");}
  //
  // }

  return *this;

}

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                            Axialline Class                             |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

AxialLine::AxialLine ():index(-1) {

  start_pt      = NULL;
  end_pt        = NULL;
  nextaxialline = NULL;

}

AxialLine::AxialLine (Cell *spt, Cell *ept):index(-1) {

  start_pt      = spt;
  end_pt        = ept;
  nextaxialline = NULL;

  this->CalcCenterPt ();

}

AxialLine::AxialLine (int offs, Cell *spt, Cell *ept):index(offs) {

  start_pt      = spt;
  end_pt        = ept;
  nextaxialline = NULL;

  this->CalcCenterPt ();

}

AxialLine::~AxialLine () {}

AxialLine & AxialLine::SetIndex (int offs) {

  index = offs;

  return *this;

}

AxialLine & AxialLine::SetNext (AxialLine *src) {

  nextaxialline = src;

  return *this;
}

AxialLine & AxialLine::setpoints (Cell *spt, Cell *ept) {

  start_pt = spt;
  end_pt   = ept;

  this->CalcCenterPt ();

  return *this;
}

AxialLine & AxialLine::AppendAxialLine (Cell *spt, Cell *ept) {

  AxialLine *pcl = NULL;

  if (!this->Next()) {

    nextaxialline = new AxialLine;
    if (!AxialLineMemoryCheck(nextaxialline)) exit(0);
    nextaxialline->setpoints(spt, ept);

  } else {

    pcl = this->Next();
    nextaxialline = new AxialLine;
    if (!AxialLineMemoryCheck(nextaxialline)) exit(0);
    nextaxialline->setpoints(spt, ept);
    nextaxialline->SetNext(pcl);

  }

  return *this;
}

AxialLine & AxialLine::AppendAxialLine (int offs, Cell *spt, Cell *ept) {

  this->AppendAxialLine(spt, ept);
  this->Next()->SetIndex(offs);

  return *this;
}

AxialLine & AxialLine::CalcCenterPt () {

  if (!start_pt) return *this;
  if (!end_pt)   return *this;

  Vector src(start_pt->Item().getdim() - 1);

  for (size_t i = 0; i < src.getdim(); i++) src.setvalue(i, (start_pt->Item()[i] + end_pt->Item()[i]) / 2.0);

  center_pt.setvector(src);

  return *this;
}

AxialLine & AxialLine::showcontents (const char *comnt) {

  char wording[1024];

  sprintf(wording,"%s(%d) %s", comnt, this->Index(), "Starting point ");
  this->Start()->showcontents(wording);

  sprintf(wording,"%s(%d) %s", comnt, this->Index(), "Ending point ");
  this->End()->showcontents(wording);

  sprintf(wording,"%s(%d) %s", comnt, this->Index(), "Center point ");
  this->Center()->showcontents(wording);

  return *this;
}

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                         Axialline Linked list                          |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

Thing_AxialLine::Thing_AxialLine ():nline(0) {

  headline = NULL;
  run = headline;

}

Thing_AxialLine::Thing_AxialLine (const Thing_AxialLine & src) {

  AxialLine *srcline, *destline;
  AxialLine *pc;
  Cell *spt, *ept;

  nline = src.Howmany();

  if (!src.Head()) {

    nline = 0;
    run = headline = NULL;

  } else {
    // desthead
    spt = src.Head()->Start();
    ept = src.Head()->End();

    headline = new AxialLine(spt, ept);
    if (!AxialLineMemoryCheck(headline)) exit(0);
    headline->SetIndex(src.Head()->Index());

    // pointing headlines;
    srcline = src.Head();
    destline = this->Head();

    // copy lines until null
    while (srcline->Next()) {

      spt = srcline->Next()->Start();
      ept = srcline->Next()->End();
      pc = new AxialLine(spt, ept);
      if (!AxialLineMemoryCheck(pc)) exit(0);
      pc->SetIndex(srcline->Next()->Index());
      destline->SetNext(pc);

      srcline  = srcline->Next();
      destline = destline->Next();

    }

  }

}

Thing_AxialLine::~Thing_AxialLine () {

  this->GoTail();

  while (this->Run()) this->RemoveLine();

}

Thing_AxialLine & Thing_AxialLine::MakeLine (Cell *spt, Cell *ept) {

  if (!headline) {
    headline = new AxialLine(spt, ept);
    if (!AxialLineMemoryCheck(headline)) exit(0);
    run = headline;

    nline++;

    return *this;
  }

  if (!this->Run()) return *this;

  run->AppendAxialLine(spt, ept);
  run = run->Next();

  nline++;

  return *this;
}

Thing_AxialLine & Thing_AxialLine::MakeLine (int offs, Cell *spt, Cell *ept) {

  if (!headline) {
    headline = new AxialLine(spt, ept);
    if (!AxialLineMemoryCheck(headline)) exit(0);
    headline->SetIndex(offs);
    run = headline;

    nline++;

    return *this;
  }

  if (!this->Run()) return *this;

  run->AppendAxialLine(offs, spt, ept);
  run = run->Next();

  nline++;

  return *this;
}

Thing_AxialLine & Thing_AxialLine::RemoveLine () {

  AxialLine *pc = headline;
  AxialLine *chold;
  AxialLine *nullpointer = NULL;

  if (!run) return *this;

  // find run-position
  while (pc) {
    if (pc->Next() == run) break;
    else                   pc = pc->Next();
  }

  // erase run-line
  chold = run->Next();
  run->SetNext(nullpointer);
  delete run;

  // hold the next-line
  if (pc) pc->SetNext(chold);
  else    pc = headline = chold;

  run = pc;

  nline--;

  return *this;
}

Thing_AxialLine & Thing_AxialLine::RemoveAllLine () {

  this->GoHead();

  while (this->Run()) this->RemoveLine();

  return *this;

}

Thing_AxialLine & Thing_AxialLine::TakeoffHeadLine () {

  nline = 0;
  headline = NULL;
  run = headline;

  return *this;

}


Thing_AxialLine & Thing_AxialLine::GoPost () {

  if (!run) return *this;

  run = run->Next();

  return *this;
}

Thing_AxialLine & Thing_AxialLine::GoPost (int nP) {

  for (size_t i = 0; i < nP; i++) this->GoPost();

  return *this;
}

Thing_AxialLine & Thing_AxialLine::GoPre () {

  AxialLine *pc = headline;

  if (!run) return *this;

  if (run == headline) return *this;

  while (pc) {
    if (pc->Next() != run) pc = pc->Next();
    else                   break;
  }

  run = pc;

  return *this;
}

Thing_AxialLine & Thing_AxialLine::GoPre (int ntimes) {

  for (size_t i = 0; i < ntimes; i++) this->GoPre();

  return *this;
}

Thing_AxialLine & Thing_AxialLine::GoHead () {

  if (!this->Head()) {run = NULL; return *this;}

  run = this->Head();

  return *this;
}

Thing_AxialLine & Thing_AxialLine::GoTail () {

  if (!this->Head()) {run = NULL; return *this;}

  run = this->Head();
  while (run->Next()) this->GoPost();

  return *this;
}

Thing_AxialLine & Thing_AxialLine::MoveRun (int offs) {

  run = this->WhichLine(offs);

  return *this;
}

Thing_AxialLine & Thing_AxialLine::MoveRun (AxialLine *src) {

  run = src;

  return *this;
}

Thing_AxialLine & Thing_AxialLine::SetHead (AxialLine *src) {

  headline = src;

  return *this;
}

int Thing_AxialLine::CountLines() {
  AxialLine *pc;
  int count;

  pc = this->Head(); count = 0;
  while (pc) {count++; pc = pc->Next();}

  return  count;
}

int Thing_AxialLine::Numbering() {
  AxialLine *pc;
  int        count;

  pc = this->Head(); count = 0;
  while (pc) {
    pc->SetIndex(count++);
    pc = pc->Next();
  }

  return count;
}

int Thing_AxialLine::Numbering(int num) {
  AxialLine *pc;
  int        count;

  pc = this->Head(); count = num;
  while (pc) {
    pc->SetIndex(count++);
    pc = pc->Next();
  }

  return count;
}

AxialLine * Thing_AxialLine::WhichLine (int idx) {

  AxialLine *pc;

  pc = this->Head();
  while (pc) {
    if (pc->Index() == idx) return pc;
    pc = pc->Next();
  }

  return pc;
}

Thing_AxialLine & Thing_AxialLine::Import (Thing_AxialLine & src, char pre_post) {
  //  Putting  [src] into [*this] at [this-run] position
  //  i.e.  [src]-cells move into [this-thing]
  //  pre_post = '+' (place them rear of [run])
  //             '-' (place them front of [run])

  AxialLine *pthis;
  int n_src, m_remain;

  // src->run is empty of cells
  if(!src.Run()) return *this;

  // this is not empty of cells but this->run is null
  if(this->Run() == NULL && this->Head() != NULL) return *this;

  // putting process

  n_src = src.Howmany();

  if (!this->Head()) {
    pthis = NULL;
    this->SetHead(src.Run());
  } else if (this->Run() == this->Head() && pre_post == '-') {
    pthis = this->Head();
    this->SetHead(src.Run());
  } else {
    if (pre_post == '-') this->GoPre();
    pthis = this->Run()->Next();
    this->Run()->SetNext(src.Run());
  }

  if (src.Run() != src.Head()) {
    src.GoPre();
    src.Run()->SetNext(NULL);
  } else   src.SetHead(NULL);

  this->GoTail();
  this->Run()->SetNext(pthis);

  m_remain = src.CountLines(); // how many cells remain in [src]
  src.SetNline(m_remain);

  nline += n_src - m_remain; // how many moving cells

  return *this;
}

Thing_AxialLine & Thing_AxialLine::Import (Thing_AxialLine & src, int nPost, char pre_post) {

  int             ntail;
  Thing_AxialLine tmpTH;
  AxialLine      *prun;

  if (!src.Run()) return *this;

  if (nPost >= src.Howmany()) {
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

Thing_AxialLine & Thing_AxialLine::showcontents (const char *comnt) {

  AxialLine *pc = headline;

  printf("%s%s%d%s\n", comnt, "[[", this->Howmany(), "]]--->");

  while (pc) {
    if (pc == run) pc->showcontents("*AxialLine");
    else           pc->showcontents(" AxialLine");
    pc = pc->Next();
  }

  return *this;
}

Thing_AxialLine & Thing_AxialLine::operator=(Thing_AxialLine & src) {
  AxialLine *srcline, *destline;
  AxialLine *pc;
  Cell      *spt, *ept;

  if (this == &src)  return *this;

  // remove all cells of this object
  this->RemoveAllLine();

  if(!src.Head()) {
    nline = 0;
    run = headline = NULL;
    return *this;
  }

  // desthead
  spt = src.Head()->Start();
  ept = src.Head()->End();
  headline = new AxialLine(spt, ept);
  if (!AxialLineMemoryCheck(headline)) exit(0);
  headline->SetIndex(src.Head()->Index());

  // pointing headlines
  srcline = src.Head();
  destline = this->Head();

  // copy  cells until null
  while (srcline->Next()) {
    spt = srcline->Next()->Start();
    ept = srcline->Next()->End();
    pc = new AxialLine(spt, ept);
    if (!AxialLineMemoryCheck(pc)) exit(0);

    pc->SetIndex(srcline->Next()->Index());
    destline->SetNext(pc);

    srcline = srcline->Next();
    destline = destline->Next();
  }

  // how many cells
  nline = src.Howmany();

  //  run position
  run = src.Run();

  return *this;
}

AxialLine & Thing_AxialLine::operator[](int offs)
{
  AxialLine *pc=NULL;

  pc = this->WhichLine(offs);
  if (!pc) {
    printf("%s%d%s\n", "no cell to the offset [", offs, "]...");
    exit(0);
  }

  return *pc;
}

#endif
