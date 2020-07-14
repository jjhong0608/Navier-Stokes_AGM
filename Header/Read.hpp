#ifndef READ
#define READ

#include "Class++.hpp"

int IsIn (char cc, const char *tokens) {
  // return 1  <------ cc happens to {tokens}
  // return 0  <------ otherwise
  const char *tmp = tokens;

  while (*tmp != '\0')
  if (cc == *tmp++) return 1;

  return 0;
}

int copyupto (char *src, char *dest, const char *end) {
  char *cc = src;
  char *pstr = dest;

  do {
    if (!IsIn(*cc, end)) *pstr++ = *cc++;
    else                 break;
  } while(*cc != '\0');

  *pstr = '\0';

  if(*cc == '\0') return  0;
  else            return  1;
}

int copyupto (const char *src, char *dest) {
  const char *cc = src;
  char *pstr = dest;

  do {
    *pstr++ = *cc++;
  } while(*cc != '\0');

  *pstr = '\0';

  if(*cc == '\0') return  0;
  else            return  1;
}

int duplicate (char cc, const char *tokens) {
  //   return n  <------ cc happens n times to {tokens}
  const char *tmp;
  int count = 0;

  tmp = tokens;

  while (*tmp != '\0')
  if (cc == *tmp++) count++;

  return count;
}

int all_IsIn (char *src, const char *tokens) {
  //   return 1  <------ {src} belong to {tokens}
  //   return 0  <------ otherwise
  char *tmp = src;

  if (*tmp == '\0' || *tmp == '\n') return 0;

  do {
    if(!IsIn(*tmp, tokens)) return 0;
    else                    tmp++;
  } while (*tmp != '\0' && *tmp != '\n');

  return 1;
}

char readingblanks(ifstream & mf, const char *blks) {
  // keep reading characters from the file [mf]
  // while it is in [blks] (the user-defined blank characters)
  // return  the first non-blank charater not in [blks]
  char cc;

  while (mf.get(cc)) if(!IsIn(cc, blks)) break;

  return cc;
}

char readingupto (ifstream & mf, const char *end) {
  // keep reading characters from the file [mf]
  // while it is NOT in [end] (the user-defined ending characters)
  // return  the first ending charater encountered in [blks]
  char cc;

  while (mf.get(cc)) if(IsIn(cc,end)) break;

  return cc;

}

int getword(char *str, ifstream & mf, const char *preblank, const char *postblank) {
  //   return value 1 --->  there is '\n' at the end of the word
  //   return value 0 --->  there is a blank except for '\n'
  char cc, *pstr = str;

  cc = readingblanks(mf, preblank);

  do {
    if (!IsIn(cc, postblank)) *pstr++ = cc;
    else                       break;
  } while (mf.get(cc));

  *pstr = '\0';

  if (cc == '\n') return  1;
  else            return  0;
}

// int getword(char *str, ifstream & mf) {
//   //   return value 1 --->  there is '\n' at the end of the word
//   //   return value 0 --->  there is a blank except for '\n'
//   char cc, *pstr = str;
//
//   cc = readingblanks(mf, " \t");
//   if (cc == '\n') {*pstr = '\0'; return  1;}
//
//   do {
//     if (!IsIn(cc, " \t\n")) *pstr++ = cc;
//     else                     break;
//   } while (mf.get(cc));
//
//   *pstr = '\0';
//
//   if(cc == '\n') return  1;
//   else           return  0;
// }

int getword (string *str, ifstream & mf) {

  mf >> *str;

  return 0;

}

int getword(char *str, ifstream & mf, const char *tokens) {
  char cc, *pstr = str;

  cc = readingblanks(mf, " \t");

  do {
    if (!IsIn(cc, " \t\n\r")) *pstr++ = cc;
    else                       break;
  } while (mf.get(cc));

  *pstr = '\0';

  // every entry in str appears in tokens->1
  if (*str == '\n') return 0;

  return all_IsIn(str, tokens);
}

int getwordupto(char *str, ifstream & mf, const char *end) {
  char cc, *pstr=str;

  cc = readingblanks(mf, " \t");

  do {
    if (!IsIn(cc, end)) {
      *pstr++ = cc;
    } else {
      break;
    }
  } while(mf.get(cc));

  *pstr = '\0';

  if(cc == '\n') return  1;
  else           return  0;
}

int rmtailblanks(char *src, const char *blks) {
  char  *p;
  int    count=0;

  p = src;
  while(*p != '\0') {
    if(!IsIn(*p, blks)) {
      count++; p++;
    } else {
      break;
    }
  }

  *p = '\0';

  return count;
}

int ReadDimension (int *nD) {

  *nD = D2;

  printf("%d%s\n\n", *nD, "D Axial line generator...");

  return 0;

}

int ReadGeoInfoName (const string Geometry_filename, ifstream & GeoInfo) {

  GeoInfo.open(Geometry_filename);

  if (!GeoInfo.is_open()) {

    printf("%s%s%s\n", "FILE [", Geometry_filename.c_str(), "] exists?");
    exit(1);

  }

  printf("%s%s%s\n", "FILE [", Geometry_filename.c_str(), "] is opened");

  return 0;

}

int ReadOutputName (const string AGL_output_file, FILE* &output) {

  output = fopen(AGL_output_file.c_str(), "w");

  if (!output) {

    printf("%s%s%s\n", "FILE [", AGL_output_file.c_str(), "] exists?");
    exit(1);

  }

  printf("%s%s%s\n", "FILE [", AGL_output_file.c_str(), "] is opened");

  return 0;

}

int ReadInputFile (const string Geometry_filename, const string AGL_output_file, ifstream & GeoInfo, FILE* &output, int *nD) {

  ReadDimension(nD);
  ReadGeoInfoName(Geometry_filename, GeoInfo);
  ReadOutputName(AGL_output_file, output);

  return 1;

}

int ReadGeoInfoNumber (ifstream & GeoInfo, int *nS, int *nR, int *nE, int *nG) {

  string words;

  while (!GeoInfo.eof()) {

    getword(&words, GeoInfo);

    if (!words.compare("SECTION"))   (*nS)++;
    if (!words.compare("REGION"))    (*nR)++;
    if (!words.compare("LINE"))      (*nE)++;
    if (!words.compare("SURFACE"))   (*nE)++;
    if (!words.compare("+"))         (*nG)++;
    if (!words.compare("-"))         (*nG)++;

  }

  printf("%s%d\t%s%d\t%s%d\t%s%d\n", "#Section = ", *nS, "#Region = ", *nR, "#Element = ", *nE, "#Sign = ", *nG);

  GeoInfo.clear();
  GeoInfo.seekg(0,ios::beg);

  return 0;

}

int GeoInfoToVec (ifstream & GeoInfo, int *nSE, int *nTE, int nDV, Vector &vtmp, char *bdc, double *bdv) {

  string words;

  getword(&words, GeoInfo);

  if (!words.compare("ENDSECTION")) return 1;

  for (size_t i = 0; i < nDV; i++) {getword(&words, GeoInfo); vtmp[i] = atof(words.c_str());}

  getword(&words, GeoInfo);

  *bdc = words[0];

  if (*bdc != 'I') getword(&words, GeoInfo); *bdv = atof(words.c_str());

  (*nSE)++;
  (*nTE)++;

  return 0;

}

int RegionInfoToVec (ifstream & GeoInfo, int *nRE, int *nTE, int nDV, string *sign, string *section_name, Vector &vtmp) {

  string words;

  if ((*nRE) == 0) for (size_t i = 0; i < nDV; i++) {getword(&words, GeoInfo); vtmp[i] = atof(words.c_str());}

  getword(&words, GeoInfo);

  if (!words.compare("ENDREGION")) return 1;

  *sign = words; getword(&words, GeoInfo); *section_name = words;

  (*nRE)++;
  (*nTE)++;

  return 0;

}

int ReadGeoInfo (ifstream & GeoInfo, Section* & section, Region* & region, int nS, int nR, int nE, int nD, int nG, int tt) {

  string words, name, nMovex, nMovey;
  int  nSE = 0, nTE = -1, nRE = 0, k = 0, nt = 0;
  char bdc;
  double bdv;
  section = new Section[nS];
  region  = new Region[nR];

  if (nD == D2) {

    Line2D *SL2D = new Line2D[nE];
    int nDV = 3 * D2; // [x1 y1] [x2 y2] [n1 n2]
    // int nDV = 2 * D2; // [x1 y1] [x2 y2]
    Vector vtmp(nDV);

    getword(&words, GeoInfo);

    while(words.compare("REGION")) {

      if (!words.compare("SECTION")) {

        printf("%s ", "SECTION");
        getword(&name, GeoInfo);
        printf("%s %s\n", name.c_str(), "Reading..."); nt++;

        getword(&nMovex, GeoInfo);
        getword(&nMovey, GeoInfo);

        nSE = 0;
        while (!GeoInfoToVec(GeoInfo, &nSE, &nTE, nDV, vtmp, &bdc, &bdv)) {

          SL2D[nTE].SetLine(vtmp[0] + tt * atof(nMovex.c_str()), vtmp[1] + tt * atof(nMovey.c_str()), vtmp[2] + tt * atof(nMovex.c_str()), vtmp[3] + tt * atof(nMovey.c_str()));
          // SL2D[nTE].SetNormal (vtmp[4], vtmp[5]);
          SL2D[nTE].SetBD(bdc, bdv);

        }

        section[k].setSection(name, nSE, &SL2D[nTE - nSE + 1]);
        k++;

      }

      readingupto(GeoInfo, "\n");
      getword(&words, GeoInfo);

    }

    nDV = 3 * D2 + 1; // [conductivity x_min y_min x_max y_max delta_x delta_y]

    vtmp.setzerovector(nDV);
    string *sign = new string[nG];
    string *section_name = new string[nG];

    k = 0;
    nTE = 0;
    int j = 0;

    while (!GeoInfo.eof()) {

      if (!words.compare("REGION")) {

        printf("%s ", "REGION");
        getword(&name, GeoInfo);
        printf("%s %s\n", name.c_str(), "Reading...");

        region[k].setName(name);

        nRE = 0;
        while (!RegionInfoToVec(GeoInfo, &nRE, &nTE, nDV, &sign[nTE], &section_name[nTE], vtmp));

        region[k].setSectionnum(nRE);
        region[k].setGridinfo(&vtmp);

        for (size_t n = 0; n < nRE; n++) {

          for (size_t i = 0; i < nS; i++) {

            if (!section_name[j].compare(section[i].getName())) {

              region[k].setSection(&section[i], n);
              region[k].setSign(sign[j], n);
              j++;
              break;

            }

          }

        }

        k++;

      }

      getword(&words, GeoInfo);

    }

  }

  // if (nD == D3) Triangle *SL3D = new Triangle[nE];

  for (size_t i = 0; i < nR; i++) region[i].make2Dgrid();

  GeoInfo.close();

  return 0;

}

int MasterMeetSlave(int nMS, Line2D* MS2D, int nSL, Line2D* SL2D, Thing & crspt, double sign) {

  int      i, nFind, nPTS = 0;
  Cell    *pc;

  // cout << endl; // by junhong, jo

  // finding intersecting points of 2D master lines with 2D slave lines
  for(nPTS = 0, i = 0 ; i < nMS ; i++) {
    //          sprintf(tmpstr,"M-Line[%d]",i);
    //          MS2D[i].showcontents(tmpstr);

    pc = crspt.Run(); // save initial [run] position
    nFind = MS2D[i].ScanCross(nSL, SL2D, crspt, sign); // MS3D[i] meets TRI3D

    nPTS += nFind;
    // cout << "L[" << i << "]-----(+" << nFind << "<" << nPTS << ">" << ")" << endl; // by junhong, jo

    if(nFind) {

      if(!pc) pc = crspt.Head();
      else    pc = pc->Next();

      while(pc) {
        pc->SetIndex(i);
        pc = pc->Next();
      }

      //             sprintf(tmpstr,"crspt[%d]",i);
      //	     crspt.showcontents(tmpstr);
    }

    //	  cout << "+---------------------------------------+" << endl << endl;
  }

  return nPTS;
}

int MasterMeetSlave(int nMS, Line3D* & MS3D, int nSL, Triangle* & TRI3D, Thing & crspt) {

  int      i, nFind, nPTS;
  Cell    *pc;

  cout << endl;

  // finding intersecting points of 3D master lines with 3D slave Triangles
  for(nPTS = 0, i = 0 ; i < nMS ; i++) {
    //          sprintf(tmpstr,"M-Line[%d]",i);
    //          MS3D[i].showcontents(tmpstr);

    pc = crspt.Run(); // save initial [run] position

    nFind = MS3D[i].ScanCross(nSL, TRI3D, crspt); // MS3D[i] meets TRI3D

    nPTS += nFind;
    cout << "L[" << i << "]-----(+" << nFind << "<" << nPTS << ">" << ")" << endl;

    if(nFind) {

      if(!pc) pc = crspt.Head();
      else    pc = pc->Next();

      while(pc) {
        pc->SetIndex(i);
        pc = pc->Next();
      }

      //             sprintf(tmpstr,"crspt[%d]",i);
      //	     crspt.showcontents(tmpstr);
    }

    // Sorting in ascending order(+) of Item[3]
    if(nFind) crspt.GoPre(nFind-1).Sorting("Item3+").GoTail();

    //	  cout << "+---------------------------------------+" << endl << endl;
  }

  return nPTS;
}

int DelDupPt (Thing & src) {

  Cell  *cel1;
  Cell  *cel2;
  Vector vec1;
  Vector vec2;

  src.GoHead();
  while (src.Run()) {

    cel1 = src.Run(); vec1 = src.Run()->Item(); src.GoPost(); if (!src.Run()) break;
    cel2 = src.Run(); vec2 = src.Run()->Item();

    if (IsEqualDouble(vec1[0], vec2[0]) && IsEqualDouble(vec1[1], vec2[1])) {

      if (cel1->Status() != cel2->Status()) cel1->SetStatus('N');

      src.RemoveCell(); //printf("%s%f\n", "(vec1 - vec2).norm() = ", (vec1 - vec2).norm());

    }

  }

  return 1;

}

Cell * CheckSamePoint (Thing & all, Cell * pt) {

  all.GoHead();
  Vector vec1;
  Vector vec2;

  vec2 = pt->Item();

  while (all.Run()) {

    vec1 = all.Run()->Item();

    if (IsEqualDouble(vec1[0], vec2[0]) && IsEqualDouble(vec1[1], vec2[1])) {

      return all.Run();

    }

    all.GoPost();

  }

  return pt;
}

int CopyTempToPerm (Thing & perm, Thing & temp, Thing_AxialLine & permline, Thing_AxialLine & templine, char XY) {

  Cell *spt, *ept;

  if (temp.Howmany() < 1) return 0;

  temp.GoHead();
  if (XY == 'X') temp.Sorting("Item0+", temp.Howmany());
  if (XY == 'Y') temp.Sorting("Item1+", temp.Howmany());
  DelDupPt(temp); temp.GoHead();

  for (size_t i = 0; i < temp.Howmany() - 1; i++) {

    spt = temp.Run(); temp.GoPost();
    ept = temp.Run();

    spt = CheckSamePoint(perm, spt);
    ept = CheckSamePoint(perm, ept);

    templine.MakeLine(spt, ept);

  }

  temp.GoHead();

  perm.GoTail();
  perm.Import(temp, temp.Howmany(), '+');

  // temp.RemoveAllCell();
  temp.TakeoffHeadCell();

  templine.GoHead();

  permline.GoTail();
  permline.Import(templine, templine.Howmany() ,'+');

  // templine.RemoveAllLine();
  templine.TakeoffHeadLine();

  return 1;

}

int Delete2DAxialLine (Region *reg, Thing_AxialLine & aline) {

  int     idx = 0, nA = aline.Howmany(), count = 0, removed = 0;
  double  seglen = 0.0;
  Vector *segend = new Vector[8];
  Line2D *check_seg = new Line2D[8];
  Thing   cross_pt;
  int kkk = 0;

  seglen
  = (reg->getGridinfo()->getvalue(3) - reg->getGridinfo()->getvalue(1))
  + (reg->getGridinfo()->getvalue(4) - reg->getGridinfo()->getvalue(2));


  for (size_t i = 0; i < reg->getSectionnum(); i++) {

    if (reg->getSign(i) > ZeroValue) {

      count = 0; nA = aline.Howmany();
      printf("\n%s%s%s%zu%s%d%s\n", "Removing ", "Axial lines outside of  ", "Section[", i + 1, "/", reg->getSectionnum(), "]");

      aline.GoHead(); idx = 0;

      while (aline.Run()) {

        segend[0].setvector(*(aline.Run()->Center())); segend[0].setvalue(0, aline.Run()->Center()->getvalue(0) + seglen); segend[0].setvalue(1, aline.Run()->Center()->getvalue(1) + seglen);
        segend[1].setvector(*(aline.Run()->Center())); segend[1].setvalue(0, aline.Run()->Center()->getvalue(0) + seglen); segend[1].setvalue(1, aline.Run()->Center()->getvalue(1) - seglen);
        segend[2].setvector(*(aline.Run()->Center())); segend[2].setvalue(0, aline.Run()->Center()->getvalue(0) - seglen); segend[2].setvalue(1, aline.Run()->Center()->getvalue(1) + seglen);
        segend[3].setvector(*(aline.Run()->Center())); segend[3].setvalue(0, aline.Run()->Center()->getvalue(0) - seglen); segend[3].setvalue(1, aline.Run()->Center()->getvalue(1) - seglen);

        segend[4].setvector(*(aline.Run()->Center())); segend[4].setvalue(0, aline.Run()->Center()->getvalue(0) + 10.0 * NearZero); segend[4].setvalue(1, aline.Run()->Center()->getvalue(1) + 10.0 * NearZero);
        segend[5].setvector(*(aline.Run()->Center())); segend[5].setvalue(0, aline.Run()->Center()->getvalue(0) - 10.0 * NearZero); segend[5].setvalue(1, aline.Run()->Center()->getvalue(1) - 10.0 * NearZero);
        segend[6].setvector(*(aline.Run()->Center())); segend[6].setvalue(0, aline.Run()->Center()->getvalue(0) - 10.0 * NearZero); segend[6].setvalue(1, aline.Run()->Center()->getvalue(1) + 10.0 * NearZero);
        segend[7].setvector(*(aline.Run()->Center())); segend[7].setvalue(0, aline.Run()->Center()->getvalue(0) + 10.0 * NearZero); segend[7].setvalue(1, aline.Run()->Center()->getvalue(1) - 10.0 * NearZero);

        check_seg[0].SetLine(*(aline.Run()->Center()), segend[0]);
        check_seg[1].SetLine(*(aline.Run()->Center()), segend[1]);
        check_seg[2].SetLine(*(aline.Run()->Center()), segend[2]);
        check_seg[3].SetLine(*(aline.Run()->Center()), segend[3]);

        check_seg[4].SetLine(segend[4], segend[5]);
        check_seg[5].SetLine(segend[6], segend[7]);

        removed = 0;

        if (MasterMeetSlave(1, &check_seg[0], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, -UnitValue) > 0) {

          if (cross_pt.Howmany() > 1) {cross_pt.GoHead(); cross_pt.Sorting("Item2+", cross_pt.Howmany()); DelDupPt(cross_pt);}

          if (cross_pt.Howmany() > 0) if (cross_pt.Head()->Status() == 'O') {aline.RemoveLine(); idx--; removed++;}

          cross_pt.RemoveAllCell();

        } else if (MasterMeetSlave(1, &check_seg[1], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, -UnitValue) > 0) {

          if (cross_pt.Howmany() > 1) {cross_pt.GoHead(); cross_pt.Sorting("Item2+", cross_pt.Howmany()); DelDupPt(cross_pt);}

          if (cross_pt.Howmany() > 0) if (cross_pt.Head()->Status() == 'O') {aline.RemoveLine(); idx--; removed++;}

          cross_pt.RemoveAllCell();

        } else if (MasterMeetSlave(1, &check_seg[2], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, -UnitValue) > 0) {

          if (cross_pt.Howmany() > 1) {cross_pt.GoHead(); cross_pt.Sorting("Item2+", cross_pt.Howmany()); DelDupPt(cross_pt);}

          if (cross_pt.Howmany() > 0) if (cross_pt.Head()->Status() == 'O') {aline.RemoveLine(); idx--; removed++;}

          cross_pt.RemoveAllCell();

        } else if (MasterMeetSlave(1, &check_seg[3], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, -UnitValue) > 0) {

          if (cross_pt.Howmany() > 1) {cross_pt.GoHead(); cross_pt.Sorting("Item2+", cross_pt.Howmany()); DelDupPt(cross_pt);}

          if (cross_pt.Howmany() > 0) if (cross_pt.Head()->Status() == 'O') {aline.RemoveLine(); idx--; removed++;}

          cross_pt.RemoveAllCell();

        } else {

          aline.RemoveLine(); idx--;

          cross_pt.RemoveAllCell();
        }

        if (removed == 0) {

          if (MasterMeetSlave(1, &check_seg[4], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, -UnitValue) > 0 && cross_pt.Howmany() > 0) {

            aline.RemoveLine(); idx--;

            cross_pt.RemoveAllCell();

          } else if (MasterMeetSlave(1, &check_seg[5], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, -UnitValue) > 0 && cross_pt.Howmany() > 0) {

            aline.RemoveLine(); idx--;

            cross_pt.RemoveAllCell();

          }

        }

        idx++;

        if (idx > 1) aline.GoPost(); if (idx > 1) kkk++;

        printf("%s%d%s%d", "\r", count++, "/", nA);

      }

      printf("\n%s\n", "All Axial lines are removed");

    } else if (reg->getSign(i) < ZeroValue) {

      count = 0; nA = aline.Howmany();
      printf("\n%s%s%s%zu%s%d%s\n", "Removing ", "Axial lines inside of  ", "Section[", i + 1, "/", reg->getSectionnum(), "]");

      aline.GoHead(); idx = 0;
      while (aline.Run()) {

        segend[0].setvector(*(aline.Run()->Center())); segend[0].setvalue(0, aline.Run()->Center()->getvalue(0) + seglen); segend[0].setvalue(1, aline.Run()->Center()->getvalue(1) + seglen);
        segend[1].setvector(*(aline.Run()->Center())); segend[1].setvalue(0, aline.Run()->Center()->getvalue(0) + seglen); segend[1].setvalue(1, aline.Run()->Center()->getvalue(1) - seglen);
        segend[2].setvector(*(aline.Run()->Center())); segend[2].setvalue(0, aline.Run()->Center()->getvalue(0) - seglen); segend[2].setvalue(1, aline.Run()->Center()->getvalue(1) + seglen);
        segend[3].setvector(*(aline.Run()->Center())); segend[3].setvalue(0, aline.Run()->Center()->getvalue(0) - seglen); segend[3].setvalue(1, aline.Run()->Center()->getvalue(1) - seglen);

        segend[4].setvector(*(aline.Run()->Center())); segend[4].setvalue(0, aline.Run()->Center()->getvalue(0) + 10.0 * NearZero); segend[4].setvalue(1, aline.Run()->Center()->getvalue(1) + 10.0 * NearZero);
        segend[5].setvector(*(aline.Run()->Center())); segend[5].setvalue(0, aline.Run()->Center()->getvalue(0) - 10.0 * NearZero); segend[5].setvalue(1, aline.Run()->Center()->getvalue(1) - 10.0 * NearZero);
        segend[6].setvector(*(aline.Run()->Center())); segend[6].setvalue(0, aline.Run()->Center()->getvalue(0) - 10.0 * NearZero); segend[6].setvalue(1, aline.Run()->Center()->getvalue(1) + 10.0 * NearZero);
        segend[7].setvector(*(aline.Run()->Center())); segend[7].setvalue(0, aline.Run()->Center()->getvalue(0) + 10.0 * NearZero); segend[7].setvalue(1, aline.Run()->Center()->getvalue(1) - 10.0 * NearZero);

        check_seg[0].SetLine(*(aline.Run()->Center()), segend[0]);
        check_seg[1].SetLine(*(aline.Run()->Center()), segend[1]);
        check_seg[2].SetLine(*(aline.Run()->Center()), segend[2]);
        check_seg[3].SetLine(*(aline.Run()->Center()), segend[3]);

        check_seg[4].SetLine(segend[4], segend[5]);
        check_seg[5].SetLine(segend[6], segend[7]);

        removed = 0;

        if (MasterMeetSlave(1, &check_seg[0], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, UnitValue) > 0) {

          if (cross_pt.Howmany() > 1) {cross_pt.GoHead(); cross_pt.Sorting("Item2+", cross_pt.Howmany()); DelDupPt(cross_pt);}

          if (cross_pt.Howmany() > 0) if (cross_pt.Head()->Status() == 'O') {aline.RemoveLine(); idx--; removed++;}

          cross_pt.RemoveAllCell();

        } else if (MasterMeetSlave(1, &check_seg[1], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, UnitValue) > 0) {

          if (cross_pt.Howmany() > 1) {cross_pt.GoHead(); cross_pt.Sorting("Item2+", cross_pt.Howmany()); DelDupPt(cross_pt);}

          if (cross_pt.Howmany() > 0) if (cross_pt.Head()->Status() == 'O') {aline.RemoveLine(); idx--; removed++;}

          cross_pt.RemoveAllCell();

        } else if (MasterMeetSlave(1, &check_seg[2], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, UnitValue) > 0) {

          if (cross_pt.Howmany() > 1) {cross_pt.GoHead(); cross_pt.Sorting("Item2+", cross_pt.Howmany()); DelDupPt(cross_pt);}

          if (cross_pt.Howmany() > 0) if (cross_pt.Head()->Status() == 'O') {aline.RemoveLine(); idx--; removed++;}

          cross_pt.RemoveAllCell();

        } else if (MasterMeetSlave(1, &check_seg[3], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, UnitValue) > 0) {

          if (cross_pt.Howmany() > 1) {cross_pt.GoHead(); cross_pt.Sorting("Item2+", cross_pt.Howmany()); DelDupPt(cross_pt);}

          if (cross_pt.Howmany() > 0) if (cross_pt.Head()->Status() == 'O') {aline.RemoveLine(); idx--; removed++;}

          cross_pt.RemoveAllCell();

        }

        if (removed == 0) {

          if (MasterMeetSlave(1, &check_seg[4], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, UnitValue) > 0) {

            aline.RemoveLine(); idx--;

            cross_pt.RemoveAllCell();

          } else if (MasterMeetSlave(1, &check_seg[5], reg->getSection(i)->getEltnum(), reg->getSection(i)->get2DElt(0), cross_pt, UnitValue) > 0) {

            aline.RemoveLine(); idx--;

            cross_pt.RemoveAllCell();

          }

        }

        idx++;

        if (idx > 1) aline.GoPost();

        printf("%s%d%s%d", "\r", count++, "/", nA);

      }

      printf("\n%s\n", "All Axial lines are removed");

    }

  }

  aline.Numbering();

  return 1;
}

int Make2DCrossPoints (Region *reg, Thing_AxialLine &xaline, Thing_AxialLine &yaline, Thing &CrossPoints, int PTSnum) {

  xaline.GoHead();
  Line2D xLINE, yLINE;
  int count = 0;

  printf("\n%s\n", "Making Cross points:");

  while (xaline.Run()) {

    xLINE.SetLine(xaline.Run()->Start()->Item()[0], xaline.Run()->Start()->Item()[1], xaline.Run()->End()->Item()[0], xaline.Run()->End()->Item()[1]);

    yaline.GoHead();

    while (yaline.Run()) {

      yLINE.SetLine(yaline.Run()->Start()->Item()[0], yaline.Run()->Start()->Item()[1], yaline.Run()->End()->Item()[0], yaline.Run()->End()->Item()[1]);

      if (!(IsEqualDouble(xaline.Run()->Start()->Item()[1], yaline.Run()->Start()->Item()[1])))
      if (!(IsEqualDouble(xaline.Run()->Start()->Item()[1], yaline.Run()->End()->Item()[1])))
      if (!(IsEqualDouble(yaline.Run()->Start()->Item()[0], xaline.Run()->Start()->Item()[0])))
      if (!(IsEqualDouble(yaline.Run()->Start()->Item()[0], xaline.Run()->End()->Item()[0])))
      if (MasterMeetSlave(1, &xLINE, 1, &yLINE, CrossPoints, 1.0) > 0) CrossPoints.Run()->SetAxialLine(xaline.Run(), yaline.Run());
      yaline.GoPost();
      printf("%s%d%s%d", "\r", ++count, "/", xaline.Howmany() * yaline.Howmany());
    }

    xaline.GoPost();
  }

  CrossPoints.Numbering(PTSnum);

  printf("\n%s\n", "Cross points are made");

  return 1;
}

int CopyCrossPoints (Thing &CrossPoints, Thing &CrossPoints_temp, int PTSnum) {

  Vector src(2);

  CrossPoints_temp.GoHead();

  while (CrossPoints_temp.Run()) {

    src.setvalue(0, CrossPoints_temp.Run()->Item()[0]); src.setvalue(1, CrossPoints_temp.Run()->Item()[1]);
    CrossPoints.MakeCell(src); CrossPoints.Run()->SetAxialLine(CrossPoints_temp.Run()->Xaxialline(), CrossPoints_temp.Run()->Yaxialline());

    CrossPoints_temp.GoPost();
  }

  CrossPoints.Numbering(PTSnum);

  return 1;
}

int Write2DCrossPoints (FILE* fp, Thing &CrossPoints, Region *reg, int *PTSnum) {

  int count = 0;

  fprintf(fp, "%s %s\n\n", "REGION", reg->getName().c_str());
  fprintf(fp, "%s %23.16e\n\n", "Material Property =", reg->getGridinfo()->getvalue(0));
  fprintf(fp, "%s %d\n\n", "# The number of cross points =", CrossPoints.Howmany());

  CrossPoints.GoHead();
  printf("\n%s\n", "Writing Cross points:");
  while (CrossPoints.Run()) {
    fprintf(fp, "%d %23.16e %23.16e\n", *PTSnum, CrossPoints.Run()->Item()[0], CrossPoints.Run()->Item()[1]); count++;
    printf("%s%d%s%d", "\r", count, "/", CrossPoints.Howmany());
    CrossPoints.Run()->SetIndex(*PTSnum);
    CrossPoints.GoPost(); (*PTSnum)++;
  }
  fprintf(fp, "\n");
  printf("\n%s\n", "Cross points are writed");

  return 1;
}

int Count2DBoundaryPoints (FILE* fp, Thing_AxialLine &xaline, Thing_AxialLine &yaline, Region *reg, int *PTSnum) {

  int idx = 0;

  xaline.GoHead();
  while (xaline.Run()) {

    if (xaline.Run()->Start()->Index() < 1 && xaline.Run()->Start()->Index() > -10) {xaline.Run()->Start()->SetIndex(-100);
      idx++;
    }

    if (xaline.Run()->End()->Index() < 1 && xaline.Run()->End()->Index() > -10) {xaline.Run()->End()->SetIndex(-100);
      idx++;
    }
    xaline.GoPost();
  }

  yaline.GoHead();
  while (yaline.Run()) {

    if (yaline.Run()->Start()->Index() < 1 && yaline.Run()->Start()->Index() > -10) {yaline.Run()->Start()->SetIndex(-100);
      idx++;
    }

    if (yaline.Run()->End()->Index() < 1 && yaline.Run()->End()->Index() > -10) {yaline.Run()->End()->SetIndex(-100);
      idx++;
    }
    yaline.GoPost();
  }

  fprintf(fp, "%s %d\n\n", "# The number of boundary points =", idx);

  return idx;
}

int Count2DBoundaryPoints (FILE* fp, Thing_AxialLine &xaline, Thing_AxialLine &yaline) {

  int idx = 0;

  xaline.GoHead();
  while (xaline.Run()) {

    if (xaline.Run()->Start()->Index() > 0) {xaline.Run()->Start()->SetIndex(-100);
      idx++;
    }

    if (xaline.Run()->End()->Index() > 0) {xaline.Run()->End()->SetIndex(-100);
      idx++;
    }
    xaline.GoPost();
  }

  yaline.GoHead();
  while (yaline.Run()) {

    if (yaline.Run()->Start()->Index() > 0) {yaline.Run()->Start()->SetIndex(-100);
      idx++;
    }

    if (yaline.Run()->End()->Index() > 0) {yaline.Run()->End()->SetIndex(-100);
      idx++;
    }
    yaline.GoPost();
  }

  fprintf(fp, "%s %d\n\n", "# The number of boundary points =", idx);

  return idx;
}

int Write2DBoundaryPoints (FILE* fp, Thing_AxialLine &xaline, Thing_AxialLine &yaline, Region *reg, int *PTSnum, int boundaryNum) {

  int count = 0;

  xaline.GoHead();
  printf("\n%s\n", "Writing Boundary points:");
  while (xaline.Run()) {

    if (xaline.Run()->Start()->Index() < 1) {

      xaline.Run()->Start()->SetIndex(*PTSnum);

      if (xaline.Run()->Start()->Cond() == 'I') {

        fprintf(fp, "%d %23.16e %23.16e %c\n", *PTSnum,
        xaline.Run()->Start()->Item()[0], xaline.Run()->Start()->Item()[1], xaline.Run()->Start()->Cond());

      } else {

        fprintf(fp, "%d %23.16e %23.16e %c %23.16e %23.16e %23.16e\n", *PTSnum,
        xaline.Run()->Start()->Item()[0], xaline.Run()->Start()->Item()[1], xaline.Run()->Start()->Cond(), xaline.Run()->Start()->Value(), xaline.Run()->Start()->Normal()[0], xaline.Run()->Start()->Normal()[1]);

      }
      count++;
      printf("%s%d%s%d", "\r", count, "/", boundaryNum);
      (*PTSnum)++;

    }

    if (xaline.Run()->End()->Index() < 1) {

      xaline.Run()->End()->SetIndex(*PTSnum);

      if (xaline.Run()->End()->Cond() == 'I') {

        fprintf(fp, "%d %23.16e %23.16e %c\n", *PTSnum,
        xaline.Run()->End()->Item()[0], xaline.Run()->End()->Item()[1], xaline.Run()->End()->Cond());

      } else {

        fprintf(fp, "%d %23.16e %23.16e %c %23.16e %23.16e %23.16e\n", *PTSnum,
        xaline.Run()->End()->Item()[0], xaline.Run()->End()->Item()[1], xaline.Run()->End()->Cond(), xaline.Run()->End()->Value(), xaline.Run()->End()->Normal()[0], xaline.Run()->End()->Normal()[1]);

      }
      count++;
      printf("%s%d%s%d", "\r", count, "/", boundaryNum);
      (*PTSnum)++;

    }
    xaline.GoPost();
  }

  yaline.GoHead();
  while (yaline.Run()) {

    if (yaline.Run()->Start()->Index() < 1) {

      yaline.Run()->Start()->SetIndex(*PTSnum);

      if (yaline.Run()->Start()->Cond() == 'I') {

        fprintf(fp, "%d %23.16e %23.16e %c\n", *PTSnum,
        yaline.Run()->Start()->Item()[0], yaline.Run()->Start()->Item()[1], yaline.Run()->Start()->Cond());

      } else {

        fprintf(fp, "%d %23.16e %23.16e %c %23.16e %23.16e %23.16e\n", *PTSnum,
        yaline.Run()->Start()->Item()[0], yaline.Run()->Start()->Item()[1], yaline.Run()->Start()->Cond(), yaline.Run()->Start()->Value(), yaline.Run()->Start()->Normal()[0], yaline.Run()->Start()->Normal()[1]);

      }
      count++;
      printf("%s%d%s%d", "\r", count, "/", boundaryNum);
      (*PTSnum)++;

    }

    if (yaline.Run()->End()->Index() < 1) {

      yaline.Run()->End()->SetIndex(*PTSnum);

      if (yaline.Run()->End()->Cond() == 'I') {

        fprintf(fp, "%d %23.16e %23.16e %c\n", *PTSnum,
        yaline.Run()->End()->Item()[0], yaline.Run()->End()->Item()[1], yaline.Run()->End()->Cond());

      } else {

        fprintf(fp, "%d %23.16e %23.16e %c %23.16e %23.16e %23.16e\n", *PTSnum,
        yaline.Run()->End()->Item()[0], yaline.Run()->End()->Item()[1], yaline.Run()->End()->Cond(), yaline.Run()->End()->Value(), yaline.Run()->End()->Normal()[0], yaline.Run()->End()->Normal()[1]);

      }
      count++;
      printf("%s%d%s%d", "\r", count, "/", boundaryNum);
      (*PTSnum)++;

    }
    yaline.GoPost();
  }
  fprintf(fp, "\n");
  printf("\n%s\n", "Boundary points are writed");

  return 1;
}

int Write2DXaxialline (FILE* fp, Thing &CrossPoints, Thing_AxialLine &xaline, Region *reg, int PTSnum) {

  int count = 0;

  fprintf(fp, "%s %d\n\n", "# The number of cross points along x-axis =", CrossPoints.Howmany() + 2 * xaline.Howmany());

  xaline.GoHead(); CrossPoints.GoHead(); printf("\n%s\n", "Writing x-axial lines:");
  while (xaline.Run()) {

    fprintf(fp, "%d ", xaline.Run()->Start()->Index());
    count++;
    printf("%s%d%s%d", "\r", count, "/", CrossPoints.Howmany() + 2 * xaline.Howmany());

    if (CrossPoints.Run()) {

      while (CrossPoints.Run()->Xaxialline() == xaline.Run()) {

        fprintf(fp, "%d ", CrossPoints.Run()->Index());
        count++;
        printf("%s%d%s%d", "\r", count, "/", CrossPoints.Howmany() + 2 * xaline.Howmany());

        CrossPoints.GoPost();

        if (!CrossPoints.Run()) break;

      }

    }

    fprintf(fp, "%d\n", xaline.Run()->End()->Index());
    count++;
    printf("%s%d%s%d", "\r", count, "/", CrossPoints.Howmany() + 2 * xaline.Howmany());

    xaline.GoPost();
  }

  fprintf(fp, "\n");
  printf("\n%s\n", "x-axial lines are writed");

  return 1;
}

int Write2DYaxialline (FILE* fp, Thing &CrossPoints, Thing_AxialLine &yaline, Region *reg, int PTSnum) {

  int count = 0;

  fprintf(fp, "%s %d\n\n", "# The number of cross points along y-axis =", CrossPoints.Howmany() + 2 * yaline.Howmany());

  yaline.GoHead(); printf("\n%s\n", "Writing y-axial lines:");
  while (yaline.Run()) {

    CrossPoints.GoHead();

    fprintf(fp, "%d ", yaline.Run()->Start()->Index());
    count++;
    printf("%s%d%s%d", "\r", count, "/", CrossPoints.Howmany() + 2 * yaline.Howmany());

    while (CrossPoints.Run()) {

      if (CrossPoints.Run()->Yaxialline() == yaline.Run()) {

        fprintf(fp, "%d ", CrossPoints.Run()->Index());
        count++;
        printf("%s%d%s%d", "\r", count, "/", CrossPoints.Howmany() + 2 * yaline.Howmany());

      }

      CrossPoints.GoPost();
    }

    fprintf(fp, "%d\n", yaline.Run()->End()->Index());
    count++;
    printf("%s%d%s%d", "\r", count, "/", CrossPoints.Howmany() + 2 * yaline.Howmany());

    yaline.GoPost();
  }

  fprintf(fp, "\n%s\n\n", "ENDREGION");

  printf("\n%s\n", "y-axial lines are writed");

  if (count != CrossPoints.Howmany() + 2 * yaline.Howmany()) {printf("%s\n", "Write2DYaxialline Error Occur"); exit(1);}

  return 1;
}


int DeleteCrossPoints (Thing & CrossPoints, Thing & BoundaryPoints) {

  int idx = 0;;

  CrossPoints.GoHead();

  while (CrossPoints.Run()) {

    BoundaryPoints.GoHead();

    while (BoundaryPoints.Run()) {

      if (IsEqualDouble(CrossPoints.Run()->Item()[0], BoundaryPoints.Run()->Item()[0]) && IsEqualDouble(CrossPoints.Run()->Item()[1], BoundaryPoints.Run()->Item()[1])) {

        CrossPoints.RemoveCell(); idx--; break;

      }

      BoundaryPoints.GoPost();

    }

    idx++;

    if (idx > 0) CrossPoints.GoPost();

  }

  CrossPoints.Numbering();

  return 1;

}

int CellMemoryCheck(Cell * src) {

  if(!src) {
    printf("%s\n", "No more memory for [Cell]...");
    return  0;
  }

  return 1;
}

int AxialLineMemoryCheck (AxialLine *src) {

  if (!src) {
    printf("%s\n", "No more memory for [AxialLine]...");
    return 0;
  }

  return 1;
}

#endif
