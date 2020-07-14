#ifndef AXIALDATA_H
#define AXIALDATA_H

#include "ControlData.hpp"

/* 생성자 */
AxialData::AxialData () {
  // 내부점의 개수
  in_pts_num = 0;
  // 경계점의 개수
  bd_pts_num = 0;
  // 모든점의 개수
  pts_num = 0;
  // 축선생성기의 결과파일의 x-axial line위의 점의 개수
  nx = 0;
  // 축선생성기의 결과파일의 y-axial line위의 점이 개수
  ny = 0;
  // x-axial line위의 점의 개수
  xaxial_num = 0;
  // y-axial line위의 점의 개수
  yaxial_num = 0;
  // y-axial line의 개수
  xxaxial_num = 0;
  // y-axial line의 개수
  yyaxial_num = 0;
}

/* 소멸자 */
AxialData::~AxialData () {
  // 점을 포함하는 축선의 index를 저장하는 변수
  for (size_t i = 0; i < pts_num; i++)     delete [] axial_index[i]; delete [] axial_index;
  // 점의 왼쪽, 오른쪽, 위쪽, 아래쪽의 점의 index를 저장하는 변수
  for (size_t i = 0; i < pts_num; i++)     delete [] EWNS_index[i];  delete [] EWNS_index;
  // 점의 좌표를 저장하는 변수
  for (size_t i = 0; i < pts_num; i++)     delete [] pts[i];         delete [] pts;
  // x-축선의 정보를 저장하는 변수
  for (size_t i = 0; i < xxaxial_num; i++) delete [] xaxial[i];      delete [] xaxial;
  // y-축선의 정보를 저장하는 변수
  for (size_t i = 0; i < yyaxial_num; i++) delete [] yaxial[i];      delete [] yaxial;
  // x-축선위의 점들의 index를 저장하는 변수
  delete [] xaxial_index;
  // y-축선위의 점들의 index를 저장하는 변수
  delete [] yaxial_index;
  // x-축선의 첫번째 점의 xaxial_index에서의 위치를 저장하는 변수
  delete [] xxaxial_index;
  // y-축선의 첫번째 점의 yaxial_index에서의 위치를 저장하는 변수
  delete [] yyaxial_index;
  // 모든점의 index에서 내부점의 index를 저장하는 변수
  delete [] ptsTOin_pts;
  // 내부점의 index에서 모든점의 index를 저장하는 변수
  delete [] in_ptsTOpts;
  // 각 점의 경계조건의 값을 저장하는 변수
  delete [] b_u;
  // 각 점의 경계조건을 저장하는 변수
  delete [] bc_u;
  // 각 점의 conductivity를 저장하는 변수
  delete [] mp_u;
}

/* 점의 좌표 */
double AxialData::Pts (int i, char xy) {
  // 에러메시지를 출력하기 위한 buf
  char buf[256];

  // index가 모든점의 개수보다 크다면 에러메시지를 출력
  if (i >= pts_num) {
    sprintf (buf, "AxialData::Pts, i (= %d) is greater than pts_num (= %d)", i, this->Pts_Num ());
    PrintError (buf);
  }
  // x좌표를 return
  if (xy == 'x' || xy == 'X') return pts[i][0];
  // y좌표를 return
  if (xy == 'y' || xy == 'Y') return pts[i][1];

  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("AxialData::Pts");
  exit (1);
}

/* 축선의 시작점과 끝점의 정보 */
double AxialData::XYaxial (char xy, int i, int j) {
  // 에러메시지를 출력하기 위한 buf
  char buf[256];

  // x-축선의 경우
  if (xy == 'x' || xy == 'X') {
    // index가 x-축선의 개수보다 크다면 에러메시지를 출력
    if (i >= xxaxial_num) {
      sprintf (buf, "AxialData::XYaxial, i (= %d) is greater than xxaxial_num (= %d)", i, this->XXYYaxial_Num ('x'));
      PrintError (buf);
    }
    // 잘못된 index인 경우 에러메시지를 출력
    if (j > 3) {
      sprintf (buf, "AxialData::XYaxial, j (= %d) must be less than 3", j);
      PrintError (buf);
    }
    // x-축선의 점의 정보를 return
    return xaxial[i][j];
  }
  // y-축선의 경우
  if (xy == 'y' || xy == 'Y') {
    // index가 y-축선의 개수보다 크다면 에러메시지를 출력
    if (i >= yyaxial_num) {
      sprintf (buf, "AxialData::XYaxial, i (= %d) is greater than yyaxial_num (= %d)", i, this->XXYYaxial_Num ('y'));
      PrintError (buf);
    }
    // 잘못된 index인 경우 에러메시지를 출력
    if (j > 3) {
      sprintf (buf, "AxialData::XYaxial, j (= %d) must be less than 3", j);
      PrintError (buf);
    }
    // y-축선의 점의 정보를 return
    return yaxial[i][j];
  }
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("AxialData::XYaxial");
  exit (1);
}

/* 경계조건의 값 */
double AxialData::Boundaryvalue (int i) {
  // 에러메시지를 출력하기 위한 buf
  char buf[256];
  // index가 모든 점의 개수보다 크다면 에러메시지를 출력
  if (i >= pts_num) {
    sprintf (buf, "AxialData::Boundaryvalue, i (= %d) is greater than pts_num (= %d)", i, this->Pts_Num ());
    PrintError (buf);
  }
  // 경계조건의 값을 return
  return b_u[i];
}

/* 각 점에서의 conductivity */
double AxialData::MaterialProperty (int i) {
  // 에러메시지를 출력하기 위한 buf
  char buf[256];
  // index가 모든 점의 개수보다 크다면 에러메시지를 출력
  if (i >= pts_num) {
    sprintf (buf, "AxialData::MaterialProperty, i (= %d) is greater than pts_num (= %d)", i, this->Pts_Num ());
    PrintError (buf);
  }
  // 각 점에서의 conductivity를 return
  return mp_u[i];
}

/* 각 점에서의 normal vector */
double AxialData::Normal (int i, char xy) {
  // 에러메시지를 출력하기 위한 buf
  char buf[256];
  // index가 모든 점의 개수보다 크다면 에러메시지를 출력
  if (i >= pts_num) {
    sprintf (buf, "AxialData::Normal, i (= %d) is greater than pts_num (= %d)", i, this->Pts_Num ());
    PrintError (buf);
  }
  // x-성분을 return하는 경우
  if (xy == 'X' || xy == 'x') return this->normal[i][0];
  // y-성분을 return하는 경우
  if (xy == 'Y' || xy == 'y') return this->normal[i][1];
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("AxialData::Normal");
  exit (1);
}

/* 각 점에서의 경계조건 */
char AxialData::Boundarycondition (int i) {
  // 에러메시지를 출력하기 위한 buf
  char buf[256];
  // index가 모든 점의 개수보다 크다면 에러메시지를 출력
  if (i >= pts_num) {
    sprintf (buf, "AxialData::Boundarycondition, i (= %d) is greater than pts_num (= %d)", i, this->Pts_Num ());
    PrintError (buf);
  }
  // 경계조건을 return
  return bc_u[i];
}

/* x-축선과 y-축선위의 점의 개수 */
int AxialData::XYaxial_Num (char xy) {
  // x-축선위의 점의 개수를 return
  if (xy == 'x' || xy == 'X') return xaxial_num;
  // y-축선위의 점의 개수를 return
  if (xy == 'y' || xy == 'Y') return yaxial_num;
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("AxialData::XYaxial_num");
  exit (1);
}

/* x-축선과 y-축선의 개수 */
int AxialData::XXYYaxial_Num (char xy) {
  // x-축선의 개수를 return
  if (xy == 'x' || xy == 'X') return xxaxial_num;
  // y-축선이 개수를 return
  if (xy == 'y' || xy == 'Y') return yyaxial_num;
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("AxialData::XXYYaxial_num");
  exit (1);
}

/* 각 점의 왼쪽, 오른쪽, 위쪽, 아래쪽의 index */
int AxialData::EWNS_Index (int i, char EWNS) {
  // 에러메시지를 출력하기 위한 buf
  char buf[256];
  // index가 모든 점의 개수보다 크다면 에러메시지를 출력
  if (i >= pts_num) {
    sprintf (buf, "AxialData::EWNS_Index, i (= %d) is greater than pts_num (= %d)", i, this->Pts_Num ());
    PrintError (buf);
  }
  // 오른쪽점의 index
  if (EWNS == 'e' || EWNS == 'E') return EWNS_index[i][0];
  // 왼쪽점의 index
  if (EWNS == 'w' || EWNS == 'W') return EWNS_index[i][1];
  // 위쪽점의 index
  if (EWNS == 'n' || EWNS == 'N') return EWNS_index[i][2];
  // 아래쪽점의 index
  if (EWNS == 's' || EWNS == 'S') return EWNS_index[i][3];

  //위치의 참조가 잘못되었을 경우 에러메시지를 출력하고 종료
  PrintError ("AxialData::EWNS_Index");
  exit (1);
}

/* x-축선과 y-축선위의 점의 index */
int AxialData::XYaxial_Index (char xy, int i) {
  // 에러메시지를 출력하기 위한 buf
  char buf[256];
  // x-축선의 경우
  if (xy == 'x' || xy == 'X') {
    // index가 x-축선위의 점의 개수보다 크다면 에러메시지를 출력
    if (i >= xaxial_num) {
      sprintf (buf, "AxialData::XYaxial_Index, i (= %d) is greater than xaxial_num (= %d)", i, this->XYaxial_Num ('x'));
      PrintError (buf);
    }
    // x-축선위의 점의 index를 return
    return xaxial_index[i];
  }
  // y-축선의 경우
  if (xy == 'y' || xy == 'Y') {
    // index가 ¥-축선위의 점의 개수보다 크다면 에러메시지를 출력
    if (i >= yaxial_num) {
      sprintf (buf, "AxialData::XYaxial_Index, i (= %d) is greater than yaxial_num (= %d)", i, this->XYaxial_Num ('y'));
      PrintError (buf);
    }
    // ¥-축선위의 점의 index를 return
    return yaxial_index[i];
  }
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("AxialData::XYaxial_Index");
  exit (1);
}

/* x-축선의 첫번째 점의 xaxial_index의 위치와 y-축선의 첫번재 점의 yaxial_index의 위치 */
int AxialData::XXYYaxial_Index (char xy, int i) {
  // 에러메시지를 출력하기 위한 buf
  char buf[256];
  // x-축선의 경우
  if (xy == 'x' || xy == 'X') {
    // index가 x-축선위의 개수보다 크다면 에러메시지를 출력
    if (i > xxaxial_num) {
      sprintf (buf, "AxialData::XXYYaxial_Index, i (= %d) is greater than xxaxial_num (= %d)", i, this->XXYYaxial_Num ('x'));
      PrintError (buf);
    }
    // x-축선의 첫번째 점의 xaxial_index의 위치를 return
    return xxaxial_index[i];
  }
  // y-축선의 경우
  if (xy == 'y' || xy == 'Y') {
    // index가 y-축선위의 개수보다 크다면 에러메시지를 출력
    if (i > yyaxial_num) {
      sprintf (buf, "AxialData::XXYYaxial_Index, i (= %d) is greater than yyaxial_num (= %d)", i, this->XXYYaxial_Num ('y'));
      PrintError (buf);
    }
    // y-축선의 첫번째 점의 xaxial_index의 위치를 return
    return yyaxial_index[i];
  }
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("XXYYaxial_Index");
  exit (1);
}

/* 모든점의 index와 내부점의 index, 모든점의 index와 phi의 값을 계산하는 점의 index의 참조 */
int AxialData::PtsTOPts (char IP, int i) {
  // 에러메시지를 출력하기 위한 buf
  char buf[256];
  // 내부점의 index에서 모든점의 index를 참조
  if (IP == 'I' || IP == 'i') {
    // index가 내부점의 개수보다 크다면 에러메시지를 출력
    if (i >= in_pts_num) {
      sprintf (buf, "AxialData::PtsTOPts, i (= %d) is greater than in_pts_num (= %d)", i, this->In_Pts_Num ());
      PrintError (buf);
    }
    // 내부점의 index에서 모든점의 index를 return
    return in_ptsTOpts[i];
  }
  // 모든점의 index에서 내부점이 index를 참조
  if (IP == 'P' || IP == 'p') {
    // index가 모든점의 개수보다 크다면 에러메시지를 출력
    if (i >= pts_num) {
      sprintf (buf, "AxialData::PtsTOPts, i (= %d) is greater than pts_num (= %d)", i, this->Pts_Num ());
      PrintError (buf);
    }
    // 모든점의 index에서 내부점의 index를 return
    return ptsTOin_pts[i];
  }
  // 모든점의 index에서 phi의 값을 계산하는 점의 index를 참조
  if (IP == 'T' || IP == 't') {
    // index가 모든점의 개수보다 크다면 에러메시지를 출력
    if (i >= pts_num) {
      sprintf (buf, "AxialData::PtsTOPts, i (= %d) is greater than pts_num (= %d)", i, this->Pts_Num ());
      PrintError (buf);
    }
    // 모든점의 index에서 phi의 값을 계산하는 점의 index를 return
    return ptsTOphi_pts[i];
  }
  // phi의 값을 계산하는 점의 index에서 모든점의 index를 참조
  if (IP == 'H' || IP == 'h') {
    // index가 phi의 값을 계산하는 점의 개수보다 크다면 에러메시지를 출력
    if (i >= phi_pts_num) {
      sprintf (buf, "AxialData::PtsTOPts, i (= %d) is greater than phi_pts_num (= %d)", i, this->Phi_Pts_Num ());
    }
    // phi의 값을 계산하는 점의 index에서 모든점의 index를 return
    return phi_ptsTOpts[i];
  }
  //참조가 잘못되었을 경우 에러메시지를 출력하고 종료
  PrintError ("AxialData::PtsTOPts");
  exit (1);
}

/* 각 점을 포함하는 x-축선과 y축선의 index */
int AxialData::Axial_Index (int i, char xy) {
  // 각 점을 포함하는 x-축선의 index를 return
  if (xy == 'X' || xy== 'x') return axial_index[i][0];
  // 각 점을 포함하는 y-축선의 index를 return
  if (xy == 'Y' || xy== 'y') return axial_index[i][1];
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("AxialData::Axial_Index");
  exit (1);
}

/* 모든점의 index와 내부점의 index, 모든점의 index와 phi의 값을 계산하는 점의 index의 저장 */
AxialData & AxialData::SetPtsTOpts (char IP, int i, int value) {
  // 에러메시지를 출력하기 위한 buf
  char buf[256];
  // 내부점의 index에서 모든점의 index를 저장
  if (IP == 'I' || IP == 'i') {
    // index가 내부점의 개수보다 크다면 에러메시지를 출력
    if (i >= in_pts_num) {
      sprintf (buf, "AxialData::SetPtsTOpts, i (= %d) is greater than in_pts_num (= %d)", i, this->In_Pts_Num ());
      PrintError (buf);
    }
    // 내부점의 index에서 모든점의 index를 저장
    in_ptsTOpts[i] = value; return *this;
  }
  // 모든점의 index에서 내부점이 index를 저장
  if (IP == 'P' || IP == 'p') {
    // index가 모든점의 개수보다 크다면 에러메시지를 출력
    if (i >= pts_num) {
      sprintf (buf, "AxialData::SetPtsTOpts, i (= %d) is greater than pts_num (= %d)", i, this->Pts_Num ());
      PrintError (buf);
    }
    // 모든점의 index에서 내부점의 index를 저장
    ptsTOin_pts[i] = value;; return *this;
  }
  // 모든점의 index에서 phi의 값을 계산하는 점의 index를 저장
  if (IP == 'T' || IP == 't') {
    // index가 모든점의 개수보다 크다면 에러메시지를 출력
    if (i >= pts_num) {
      sprintf (buf, "AxialData::SetPtsTOpts, i (= %d) is greater than pts_num (= %d)", i, this->Pts_Num ());
      PrintError (buf);
    }
    // 모든점의 index에서 phi의 값을 계산하는 점의 index를 저장
    ptsTOphi_pts[i] = value;; return *this;
  }
  // phi의 값을 계산하는 점의 index에서 모든점의 index를 저장
  if (IP == 'H' || IP == 'h') {
    // index가 phi의 값을 계산하는 점의 개수보다 크다면 에러메시지를 출력
    if (i >= phi_pts_num) {
      sprintf (buf, "AxialData::SetPtsTOpts, i (= %d) is greater than phi_pts_num (= %d)", i, this->Phi_Pts_Num ());
    }
    // phi의 값을 계산하는 점의 index에서 모든점의 index를 저장
    phi_ptsTOpts[i] = value; return *this;
  }
  //참조가 잘못되었을 경우 에러메시지를 출력하고 종료
  PrintError ("AxialData::SetPtsTOpts");
  exit (1);
}

/* 모든점의 index와 내부점의 index, 모든점의 index와 phi의 값을 계산하는 점의 index를 저장하는 변수의 allocation */
AxialData & AxialData::AllocatePhipts (int phinum) {
  // phi의 값을 계산하는 점의 개수를 저장
  phi_pts_num = phinum;
  // 모든점의 index에서 내부점이 index를 저장하는 변수의 선언
  ptsTOin_pts = new int[this->Pts_Num ()];
  // 내부점의 index에서 모든점의 index를 저장하는 변수의 선언
  in_ptsTOpts = new int[this->In_Pts_Num ()];
  // 모든점의 index에서 phi의 값을 계산하는 점의 index를 저장하는 변수의 선언
  ptsTOphi_pts = new int[this->Pts_Num ()];
  // phi의 값을 계산하는 점의 index에서 모든점의 index를 저장하는 변수의 선언
  phi_ptsTOpts = new int[this->Phi_Pts_Num ()];

  // 각 변수의 값을 -1로 초기화
  for (size_t i = 0; i < this->Pts_Num (); i++)     {SetPtsTOpts ('P', i, -1); SetPtsTOpts ('T', i, -1);}
  for (size_t i = 0; i < this->In_Pts_Num (); i++)  {SetPtsTOpts ('I', i, -1);}
  for (size_t i = 0; i < this->Phi_Pts_Num (); i++) {SetPtsTOpts ('H', i, -1);}

  return *this;
}

/* 각 점의 왼쪽, 오른쪽, 위쪽, 아래쪽의 index를 저장하는 변수의 초기화 */
AxialData & AxialData::SortEWNS () {
  for (size_t i = 0; i < this->Pts_Num (); i++) {
    // 왼쪽의 index를 저장하는 변수의 초기화
    if (this->EWNS_Index (i, 'E') == i) EWNS_index[i][0] = -1;
    // 오른쪽의 index를 저장하는 변수의 초기화
    if (this->EWNS_Index (i, 'W') == i) EWNS_index[i][1] = -1;
    // 위쪽의 index를 저장하는 변수의 초기화
    if (this->EWNS_Index (i, 'N') == i) EWNS_index[i][2] = -1;
    // 아래쪽의 index를 저장하는 변수의 초기화
    if (this->EWNS_Index (i, 'S') == i) EWNS_index[i][3] = -1;
  }

  return *this;
}

/* AxialData를 내보내기 */
AxialData & AxialData::ExportAxialData (ControlData *cdat) {

  // 점의 좌표정보를 내보내기
  FILE *pts_output = fopen ("pts.dat", "w");
  if (pts_output == NULL) {
    printf ("====== Some error has occured ======\n");
    printf ("     Failed to open \"pts.dat\"     \n");
    printf ("====================================\n");
    exit (1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf (pts_output, "%9lu\t%23.16e\t%23.16e\n", i, pts[i][0], pts[i][1]);
  fclose (pts_output);

  // 점의 경계조건의 값을 내보내기
  FILE *b_u_output = fopen ("b_u.dat", "w");
  if (b_u_output == NULL) {
    printf ("====== Some error has occured ======\n");
    printf ("     Failed to open \"b_u.dat\"     \n");
    printf ("====================================\n");
    exit (1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf (b_u_output, "%9lu\t%23.16e\n", i, b_u[i]);
  fclose (b_u_output);

  // 점의 경계조건을 내보내기
  FILE *bc_u_output = fopen ("bc_u.dat", "w");
  if (bc_u_output == NULL) {
    printf ("======= Some error has occured ======\n");
    printf ("     Failed to open \"bc_u.dat\"     \n");
    printf ("=====================================\n");
    exit (1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf (bc_u_output, "%9lu\t%c\n", i, bc_u[i]);
  fclose (bc_u_output);

  // 점의 왼쪽, 오른쪽, 위쪽, 아래쪽의 index를 내보내기
  FILE *EWNS_index_output = fopen ("EWNS_index.dat", "w");
  if (EWNS_index_output == NULL) {
    printf ("========== Some error has occured =========\n");
    printf ("     Failed to open \"EWNS_index.dat\"     \n");
    printf ("===========================================\n");
    exit (1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf (EWNS_index_output, "%9lu\t%9d\t%9d\t%9d\t%9d\n", i, EWNS_index[i][0], EWNS_index[i][1], EWNS_index[i][2], EWNS_index[i][3]);
  fclose (EWNS_index_output);

  // x-축선의 점들의 index를 내보내기
  FILE *xaxial_index_output = fopen ("xaxial_index.dat", "w");
  if (xaxial_index_output == NULL) {
    printf ("=========== Some error has occured ==========\n");
    printf ("     Failed to open \"xaxial_index.dat\"     \n");
    printf ("=============================================\n");
    exit (1);
  }
  for (size_t i = 0; i < xaxial_num; i++)
  fprintf (xaxial_index_output, "%9lu\t%9d\n", i, xaxial_index[i]);
  fclose (xaxial_index_output);

  // y-축선의 점들의 index를 내보내기
  FILE *yaxial_index_output = fopen ("yaxial_index.dat", "w");
  if (yaxial_index_output == NULL) {
    printf ("=========== Some error has occured ==========\n");
    printf ("     Failed to open \"yaxial_index.dat\"     \n");
    printf ("=============================================\n");
    exit (1);
  }
  for (size_t i = 0; i < yaxial_num; i++)
  fprintf (yaxial_index_output, "%9lu\t%9d\n", i, yaxial_index[i]);
  fclose (yaxial_index_output);

  // x-축선의 첫번째 점의 xaxial_index의 위치를 내보내기
  FILE *xxaxial_index_output = fopen ("xxaxial_index.dat", "w");
  if (xxaxial_index_output == NULL) {
    printf ("=========== Some error has occured ===========\n");
    printf ("     Failed to open \"xxaxial_index.dat\"     \n");
    printf ("==============================================\n");
    exit (1);
  }
  for (size_t i = 0; i < xxaxial_num; i++)
  fprintf (xxaxial_index_output, "%9lu\t%9d\n", i, xxaxial_index[i]);
  fclose (xxaxial_index_output);

  // y-축선의 첫번째 점의 yaxial_index의 위치를 내보내기
  FILE *yyaxial_index_output = fopen ("yyaxial_index.dat", "w");
  if (yyaxial_index_output == NULL) {
    printf ("=========== Some error has occured ===========\n");
    printf ("     Failed to open \"xxaxial_index.dat\"     \n");
    printf ("==============================================\n");
    exit (1);
  }
  for (size_t i = 0; i < yyaxial_num; i++)
  fprintf (yyaxial_index_output, "%9lu\t%9d\n", i, yyaxial_index[i]);
  fclose (yyaxial_index_output);

  // x-축선의 시작점과 끝점의 정보를 내보내기
  FILE *xaxial_output = fopen (cdat->Xaxial ().c_str (), "w");
  if (xaxial_output == NULL) {
    printf ("======== Some error has occured =======\n");
    printf ("     Failed to open \"xaxial.dat\"     \n");
    printf ("=======================================\n");
    exit (1);
  }
  for (size_t i = 0; i < xxaxial_num; i++)
  fprintf (xaxial_output, "%9lu\t%23.16e\t%23.16e\t%23.16e\n", i, xaxial[i][0], xaxial[i][1], xaxial[i][2]);
  fclose (xaxial_output);

  // y-축선의 시작점과 끝점의 정보를 내보내기
  FILE *yaxial_output = fopen (cdat->Yaxial ().c_str (), "w");
  if (yaxial_output == NULL) {
    printf ("======== Some error has occured =======\n");
    printf ("     Failed to open \"yaxial.dat\"     \n");
    printf ("=======================================\n");
    exit (1);
  }
  for (size_t i = 0; i < yyaxial_num; i++)
  fprintf (yaxial_output, "%9lu\t%23.16e\t%23.16e\t%23.16e\n", i, yaxial[i][0], yaxial[i][1], yaxial[i][2]);
  fclose (yaxial_output);

  // 모든점의 index에서 내부점의 index의 참조를 내보내기
  FILE *ptsTOin_pts_output = fopen ("ptsTOin_pts.dat", "w");
  if (ptsTOin_pts_output == NULL) {
    printf ("========== Some error has occured ==========\n");
    printf ("     Failed to open \"ptsTOin_pts.dat\"     \n");
    printf ("============================================\n");
    exit (1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf (ptsTOin_pts_output, "%9lu\t%9d\n", i, ptsTOin_pts[i]);
  fclose (ptsTOin_pts_output);

  // 모든점의 index에서 phi를 계산하는 점들의 index의 참조를 내보내기
  FILE *ptsTOphi_pts_output = fopen ("ptsTOphi_pts.dat", "w");
  if (ptsTOphi_pts_output == NULL) {
    printf ("========== Some error has occured ==========\n");
    printf ("     Failed to open \"ptsTOin_pts.dat\"     \n");
    printf ("============================================\n");
    exit (1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf (ptsTOphi_pts_output, "%9lu\t%9d\n", i, ptsTOphi_pts[i]);
  fclose (ptsTOphi_pts_output);

  // 내부점의 index에서 모든점의 index의 참조를 내보내기
  FILE *in_ptsTOpts_output = fopen ("in_ptsTOpts.dat", "w");
  if (in_ptsTOpts_output == NULL) {
    printf ("========== Some error has occured ==========\n");
    printf ("     Failed to open \"in_ptsTOpts.dat\"     \n");
    printf ("============================================\n");
    exit (1);
  }
  for (size_t i = 0; i < in_pts_num; i++)
  fprintf (in_ptsTOpts_output, "%9lu\t%9d\n", i, in_ptsTOpts[i]);
  fclose (in_ptsTOpts_output);

  // phi를 계산하는 점등릐 index에서 모든점의 index의 참조를 내보내기
  FILE *phi_ptsTOpts_output = fopen ("phi_ptsTOpts.dat", "w");
  if (phi_ptsTOpts_output == NULL) {
    printf ("========== Some error has occured ==========\n");
    printf ("     Failed to open \"phi_ptsTOpts.dat\"     \n");
    printf ("============================================\n");
    exit (1);
  }
  for (size_t i = 0; i < phi_pts_num; i++)
  fprintf (phi_ptsTOpts_output, "%9lu\t%9d\n", i, phi_ptsTOpts[i]);
  fclose (phi_ptsTOpts_output);

  return *this;
}

/* 축선생성기로부터 만들어진 점과 축선의 정보들을 읽어들이는 모듈 */
AxialData & AxialData::LoadAxialData (string AxialFile_input) {
  // 축선생성기로부터 만들어진 점과 축선의 정보들을 저장한 파일을 열기
  ifstream AxialFile (AxialFile_input);

  // 파일로부터 한 줄을 읽어서 저장하는 변수
  string inputString;

  // 필요없는 문자열을 임시로 저장하는 변수
  string tempstring[3];

  // 내부점과 경계점, x-축선위의 점과, y-축선위의 점을 읽어들이기 위한 표시
  int ind_int = 0;

  // 각각의 region을 읽어들이기 위한 표시
  int ind_rgn = 0;

  // x-축선과 y-축선의 개수를 세기위한 변수
  int Line_num_x = 0, Line_num_y = 0;

  // x-축선과 y-축선의 점들의 개수의 참조를 위한 변수
  int ix = 0, jy = 0, tmp = 0;

  // x-축선과 y-축선위의 점의 index를 저장하기 위한 변수
  int pt_present = 0, pt_previous = 0;

  // 이전의 모든점의 개수와 이전의 내부점의 개수
  int Opts_num = 0, Oin_pts_num = 0;

  // x-축선과 y-축선위의 점의 개수를 세기위한 변수
  int Nx = 0, Ny = 0;

  // conductivity를 읽어들이기 위한 변수
  double mp;

  // 에러메시지를 출력하기 위한 buf
  char buf[256];

  // 축선생성기로부터 만들어진 파일을 열지못했을 경우에 에러메시지를 출력하고 종료
  if (AxialFile.is_open () == false) {
    sprintf (buf, "No Axial Data file: %s\n\tPlease check Axial data file name", AxialFile_input.c_str ());
    PrintError (buf);
  }

  // 축선생성기로부터 만들어진 파일을 성공적으로 열었다는 표시
  printf ("Axial file: \"%s\" open\n", AxialFile_input.c_str ());

  // 축선생성기로부터 만들어진 파일을 한줄씩 읽어서 만들어야 할 변수의 크기를 측정한다
  while (!AxialFile.eof ()) {
    // 축선생성기로부터 만들어진 파일에서 한줄을 읽어서 inputString에 저장
    getline (AxialFile, inputString);

    // inputString의 크기가 0이라면 다시 한번 더 읽어서 inputString에 저장
    if (inputString.size () == 0) getline (AxialFile, inputString);

    // 읽어들인 줄에서 "REGION"이라는 단어를 찾으면 진입
    if (inputString.find ("REGION ") != string::npos) {
      // 읽어들인 region의 개수를 더한다.
      ind_rgn += 1;

      // 파일로부터 한줄을 읽어서 inputString에 저장
      getline (AxialFile, inputString);

      // "Material",  "Property",  "="를 각각 tempstring에 저장하고 conductivity를 mp에 저장한다.
      AxialFile >> tempstring[0] >> tempstring[1] >> tempstring[2] >> mp;
    }

    // 읽어들인 줄에서 "ENDREGION"이라는 단어를 찾으면 ind_int를 0으로 초기화
    if (inputString.find ("ENDREGION") != string::npos) ind_int = 0;

    // 읽어들인 줄에서 "="을 찾으면 진입
    if (inputString.find ("=") != string::npos) {
      // ind_int에 1을 더한다
      ind_int += 1;

      // ind_int의 값에 따라서 switch문에 진입 (1: 내부점, 2: 경계점, 3: x-축선위 점, 4: y-축선의 점)
      switch (ind_int) {

        // 내부점의 개수를 읽어들인다
        case 1:

        // "="의 뒤에 있는 내부점의 개수를 in_pts_num에 더한다
        in_pts_num += stoi (inputString.substr (inputString.find ("=") + 2, inputString.size ()));

        // 지금까지의 내부점의 개수와 경계점의 개수를 더해 Oin_pts_num에 저장
        Oin_pts_num = in_pts_num + bd_pts_num;
        break;

        // 경계점의 개수를 읽어들인다
        case 2:

        // "="의 뒤에 있는 경계점의 개수를 bd_pts_num에 더한다
        bd_pts_num += stoi (inputString.substr (inputString.find ("=") + 2, inputString.size ()));

        // 지금까지의 내부점의 개수와 경계점의 개수를 더해, pts_num에 저장
        pts_num = in_pts_num + bd_pts_num;

        // Opts_num에 pts_num의 값을 저장
        Opts_num = pts_num;
        break;

        // x-축선위의 점의 개수를 읽어들인다
        case 3:

        // "="의 뒤에 있는 x-축선위의 점의 개수를 nx에 더한다
        nx += stoi (inputString.substr (inputString.find ("=") + 2, inputString.size ()));

        // 첫번째 region인 경우, Line_num_x를 0으로 초기화
        if (ind_rgn == 1 ) Line_num_x = 0;
        break;

        // y-축선위의 점의 개수를 읽어들인다
        case 4:

        // "="의 뒤에 있는 y-축선위의 점의 개수를 ny에 더한다
        ny += stoi (inputString.substr (inputString.find ("=") + 2, inputString.size ()));

        // 첫번째 region인 경우, Line_num_y를 0으로 초기화
        if (ind_rgn == 1 ) Line_num_y = 0;
        break;
      }
    }

    // 읽어들인 줄에서 "="을 찾기 못했을 경우 (x-축선과 y-축선의 개수를 센다)
    else {
      // ind_int의 값에 따라서 switch문에 진입 (3: x-축선위 점, 4: y-축선의 점)
      switch (ind_int) {
        // x-축선의 개수를 센다
        case 3:
        Line_num_x += 1;
        break;

        // y-축선의 개수를 센다
        case 4:
        Line_num_y += 1;
        break;
      }
    }
  }

  // 파일에서 커서의 위치를 다시 처음으로 되돌린다
  AxialFile.clear ();
  AxialFile.seekg (0, ios::beg);

  // 파일로부터 읽어들인 변수의 크기를 이용해서 변수들을 선언
  // 각 점의 좌표를 저장하는 변수
  pts  = new double*[this->pts_num]; for (size_t i = 0; i < this->pts_num; i++) pts[i] = new double[2];

  // 각 점의 경계조건의 값을 저장하는 변수
  b_u  = new double[this->pts_num];

  // 각 점의 경계조건을 저장하는 변수 (C: 내부점, I: Interface, D: Dirichlet, N: Neumann, F: Infinity, S: Singularity)
  bc_u = new char[this->pts_num];

  // 각 점의 conductivity를 저장하는 변수
  mp_u = new double[this->pts_num];

  // 각 점의 normal vector를 저장하는 변수
  normal = new double*[this->pts_num]; for (size_t i = 0; i < this->pts_num; i++) normal[i] = new double[2];

  // 각 점의 왼쪽, 오른쪽, 위쪽, 아래쪽의 index를 저장하는 변수
  EWNS_index = new int*[this->pts_num]; for (size_t i = 0; i < this->pts_num; i++) EWNS_index[i] = new int[4];

  // 각 점을 포함하는 x-축선과 y-축선의 index를 저장하는 변수
  axial_index = new int*[this->pts_num]; for (size_t i = 0; i < this->pts_num; i++) axial_index[i] = new int[2];

  // x-축선의 시작점과 끝점의 정보를 저장하는 변수 (시작점과 끝점의 y-좌표, 시작점의 x-좌표, 끝점의 x-좌표 순서로 저장)
  xaxial = new double*[Line_num_x]; for (size_t i = 0; i < Line_num_x; i++) xaxial[i] = new double[3];

  // y-축선의 시작점과 끝점의 정보를 저장하는 변수 (시작점과 끝점의 x-좌표, 시작점의 y-좌표, 끝점의 y-좌표 순서로 저장)
  yaxial = new double*[Line_num_y]; for (size_t i = 0; i < Line_num_y; i++) yaxial[i] = new double[3];

  // x-축선위의 점들의 index를 저장하는 변수
  xaxial_index = new int[this->nx];

  // y-축선위의 점들의 index를 저장하는 변수
  yaxial_index = new int[this->ny];

  // 각 x-축선의 첫번째 점의 xaxial_index의 위치를 저장하는 변수 (마지막은 저장된 원소의 총 개수)
  xxaxial_index = new int[Line_num_x + 1];

  // 각 y-축선의 첫번째 점의 yaxial_index의 위치를 저장하는 변수 (마지막은 저장된 원소의 총 개수)
  yyaxial_index = new int[Line_num_y + 1];

  // 각 점의 경계조건을 'C'로 초기화
  for (size_t i = 0; i < this->pts_num; i++) bc_u[i] = 'C';

  // 각 점의 normal vector를 0으로 초기화
  for (size_t i = 0; i < this->pts_num; i++) for (size_t j = 0; j < 2; j++) normal[i][j] = ZeroValue;

  // 각 점의 왼쪽, 오른쪽, 위쪽, 아래쪽의 index를 -1로 초기화
  for (size_t i = 0; i < this->pts_num; i++) for (size_t j = 0; j < 4; j++) EWNS_index[i][j] = -1;

  // 각 점을 포함하는 x-축선과 y-축선의 index를 -1로 초기화
  for (size_t i = 0; i < this->pts_num; i++) for (size_t j = 0; j < 2; j++) axial_index[i][j] = -1;

  // 내부점의 개수를 0으로 초기화
  in_pts_num = 0;

  // 경계점의 개수를 0으로 초기화
  bd_pts_num = 0;

  // 모든점의 개수를 0으로 초기화
  pts_num = 0;

  // x-축선위의 점의 개수를 0으로 초기화
  xaxial_num = nx;

  // y-축선위의 점의 개수를 0으로 초기화
  yaxial_num = ny;

  // x-축선위의 점의 개수를 0으로 초기화
  nx = 0;

  // y-축선위의 점의 개수를 0으로 초기화
  ny = 0;

  // x-축선의 개수에 위에서 읽어들인 Line_num_x를 저장
  xxaxial_num = Line_num_x;

  // y-축선의 개수에 위에서 읽어들인 Line_num_y를 저장
  yyaxial_num = Line_num_y;

  // 내부점과 경계점, x-축선위의 점과, y-축선위의 점을 읽어들이기 위한 표시를 0으로 초기화
  ind_int = 0;

  // 각각의 region을 읽어들이기 위한 표시를 0으로 초기화
  ind_rgn = 0;

  // x-축선의 개수와 y-축선의 개수를 0으로 초기화
  Line_num_x = 0, Line_num_y = 0;

  // x-축선과 y-축선의 점들의 개수의 참조를 위한 변수를 0으로 초기화
  ix = 0, jy = 0, tmp = 0;

  // x-축선과 y-축선위의 점의 index를 저장하기 위한 변수를 0으로 초기화
  pt_present = 0, pt_previous = 0;

  // 이전의 모든점의 개수와 이전의 내부점의 개수를 0으로 초기화
  Opts_num = 0, Oin_pts_num = 0;


  // 축선생성기로부터 만들어진 파일을 한줄씩 읽어서 데이터를 읽어들인다
  while (!AxialFile.eof ()) {
    // 축선생성기로부터 만들어진 파일에서 한줄을 읽어서 inputString에 저장
    getline (AxialFile, inputString);

    // inputString의 크기가 0이라면 다시 한번 더 읽어서 inputString에 저장
    while (inputString.size () == 0 && !AxialFile.eof ()) getline (AxialFile, inputString);

    // 읽어들인 줄에서 "REGION"이라는 단어를 찾으면 진입
    if (inputString.find ("REGION ") != string::npos) {
      // 읽어들인 region의 개수를 더한다.
      ind_rgn += 1;

      // ㅇ번째 region을 읽어들인다는 표시
      printf ("REGION %d reading...\n", ind_rgn);

      // 파일로부터 한줄을 읽어서 inputString에 저장
      getline (AxialFile, inputString);

      // "Material",  "Property",  "="를 각각 tempstring에 저장하고 conductivity를 mp에 저장한다.
      AxialFile >> tempstring[0] >> tempstring[1] >> tempstring[2] >> mp;
    }

    // 읽어들인 줄에서 "ENDREGION"이라는 단어를 찾으면 ind_int를 0으로 초기화
    if (inputString.find ("ENDREGION") != string::npos) ind_int = 0;

    // 읽어들인 줄에서 "="을 찾으면 진입
    if (inputString.find ("=") != string::npos) {
      // ind_int에 1을 더한다
      ind_int += 1;

      // ind_int의 값에 따라서 switch문에 진입 (1: 내부점, 2: 경계점, 3: x-축선위 점, 4: y-축선의 점)
      switch (ind_int) {
        // 내부점의 정보를 읽어들인다
        case 1:

        // "="의 뒤에 있는 내부점의 개수를 in_pts_num에 더한다
        in_pts_num += stoi (inputString.substr (inputString.find ("=") + 2, inputString.size ()));

        // 내부점의 index와 x-좌표, y-좌표, conductivity를 읽어서 저장
        for (size_t i = Opts_num; i < in_pts_num + bd_pts_num; i++) {AxialFile >> tmp >> pts[i][0] >> pts[i][1]; mp_u[i] = mp;}

        // 지금까지의 내부점의 개수와 경계점의 개수를 더해 Oin_pts_num에 저장
        Oin_pts_num = in_pts_num + bd_pts_num;
        break;

        // 경계점의 정보를 읽어들인다
        case 2:

        // "="의 뒤에 있는 경계점의 개수를 bd_pts_num에 더한다
        bd_pts_num += stoi (inputString.substr (inputString.find ("=") + 2, inputString.size ()));

        // 지금까지의 내부점의 개수와 경계점의 개수를 더해, pts_num에 저장
        pts_num = in_pts_num + bd_pts_num;

        // 이 region의 경계점을 제외한 모든점에서부터 이 region의 경계점을 포함한 모든점의 개수까지 반복
        for (size_t i = Oin_pts_num; i < pts_num; i++) {
          // 경계점의 index와 x-좌표, y-좌표, 경계조건을 읽어서 저장
          AxialFile >> tmp >> pts[i][0] >> pts[i][1] >> bc_u[i];

          // 경게조건이 interface가 아니라면 경계조건의 값을 읽어서 저장
          if (bc_u[i] != 'I') AxialFile >> b_u[i] >> normal[i][0] >> normal[i][1];

          /*          TEMP          */
          // if (IsEqualDouble (pts[i][1], -3.0E0)) bc_u[i] = 'D';
          // if (IsEqualDouble (pts[i][0],  3.0E0)) bc_u[i] = 'D';
          // if (IsEqualDouble (pts[i][0], -3.0E0)) bc_u[i] = 'N', normal[i][0] = pts[i][0] / abs (pts[i][0]), normal[i][1] = ZeroValue;
          // if (IsEqualDouble (pts[i][1],  3.0E0)) bc_u[i] = 'N', normal[i][1] = pts[i][1] / abs (pts[i][1]), normal[i][0] = ZeroValue;
          // if (IsEqualDouble (pts[i][0], -3.0E0)) bc_u[i] = 'N', normal[i][0] = -1.0E0, normal[i][1] = ZeroValue;
          // if (IsEqualDouble (pts[i][1],  3.0E0)) bc_u[i] = 'N', normal[i][1] = 1.0E0, normal[i][0] = ZeroValue;
          /*          TEMP          */

          // 각 경계점의 conductivity를 저장
          mp_u[i] = mp;
        }

        // Opts_num에 pts_num의 값을 저장
        Opts_num = pts_num;
        break;

        // x-축선의 정보를 읽어들인다
        case 3:
        // "="의 뒤에 있는 x-축선위의 점의 개수를 nx에 더한다
        nx += stoi (inputString.substr (inputString.find ("=") + 2, inputString.size ()));

        // 첫번째 region인 경우, ix 0으로 초기화
        if (ind_rgn == 1 ) ix = 0;

        // 현재 region에서의 x-축선위의 점의 개수보다 작은 동안 반복
        while (Nx < nx) {
          // x-축선의 첫번째 줄의 첫번째 점의 index를 읽어들인다
          AxialFile >> pt_previous;

          // xaxial_index에 점의 index를 저장
          this->xaxial_index[Nx] = pt_previous;

          // xxaxial_index에 첫번째 점의 xaxial_index에서의 위치를 추가
          this->xxaxial_index[ix] = Nx;

          // 읽어들인 점의 개수를 추가
          Nx += 1;

          // 읽어들인 첫번째 점이 내부점이라면 에러메시지를 출력
          if (this->bc_u[pt_previous] == 'C') {
            sprintf (buf, "%s\n\t%s", "The first point of the x-axial line is the cross point", "Please check the Axial Data file");
            PrintError (buf);
          }

          // x-축선의 시작점의 y-좌표를 저장
          this->xaxial[ix][0] = pts[pt_previous][1];

          // x-축선의 시작점의 x-좌표를 저장
          this->xaxial[ix][1] = pts[pt_previous][0];

          // x-축선의 시작점의 왼쪽의 점을 자기자신으로 함
          this->EWNS_index[pt_previous][1] = pt_previous;

          // x-축선의 시작점을 포함하는 x-축선의 index를 저장
          this->axial_index[pt_previous][0] = ix;

          // x-축선의 두번재 점의 index를 읽어서 현재점의 index에 저장
          AxialFile >> pt_present;

          // xaxial_index에 점의 index를 저장
          this->xaxial_index[Nx] = pt_present;

          // 읽어들인 점의 개수를 추가
          Nx += 1;

          // 읽어들인 점이 경계점인 동안 반복
          while (this->bc_u[pt_present] == 'C') {
            // 이전에 읽어들인 점의 오른쪽점의 index에 현재 읽어들인 점의 index를 저장
            this->EWNS_index[pt_previous][0] = pt_present;

            // 현재 읽어들인 점의 왼쪽점의 index에 이전에 읽어들인 점의 index를 저장
            this->EWNS_index[pt_present][1] = pt_previous;

            // 현재점을 포함하는 x-축선의 index를 저장
            this->axial_index[pt_present][0] = ix;

            // 이전에 읽어들인 점에 현재점의 index를 저장
            pt_previous = pt_present;

            // x-축선의 다음점의 index를 읽어서 현재점의 index에 저장
            AxialFile >> pt_present;

            // xaxial_index에 점의 index를 저장
            this->xaxial_index[Nx] = pt_present;

            // 읽어들인 점의 개수를 추가
            Nx += 1;
          }

          // x-축선의 끝점의 오른쪽점을 자기 자신으로 저장
          this->EWNS_index[pt_present][0] = pt_present;

          // x-축선의 끝점의 x-좌표를 저장
          this->xaxial[ix][2] = pts[pt_present][0];

          // x-축선의 끝점을 포함하는 x-축선의 index를 저장
          this->axial_index[pt_present][0] = ix;

          // 이전에 읽어들인 점의 오른쪽점의 index에 x-축선의 끝점의 index를 저장
          this->EWNS_index[pt_previous][0] = pt_present;

          // x-축선의 끝점의 왼쪽점의 index에 이전에 읽어들인 점의 index를 저장
          this->EWNS_index[pt_present][1] = pt_previous;

          // x-축선의 index에 1을 더한다
          ix += 1;
        }
        break;

        // y-축선의 정보를 읽어들인다
        case 4:

        // "="의 뒤에 있는 y-축선위의 점의 개수를 ny에 더한다
        ny += stoi (inputString.substr (inputString.find ("=") + 2, inputString.size ()));

        // 첫번째 region인 경우, jy 0으로 초기화
        if (ind_rgn == 1 ) jy = 0;

        // 현재 region에서의 y-축선위의 점의 개수보다 작은 동안 반복
        while (Ny < ny) {
          // y-축선의 첫번째 줄의 첫번째 점의 index를 읽어들인다
          AxialFile >> pt_previous;

          // yaxial_index에 점의 index를 저장
          this->yaxial_index[Ny] = pt_previous;

          // yyaxial_index에 첫번째 점의 yaxial_index에서의 위치를 추가
          this->yyaxial_index[jy] = Ny;

          // 읽어들인 점의 개수를 추가
          Ny += 1;

          // 읽어들인 첫번째 점이 내부점이라면 에러메시지를 출력
          if (this->bc_u[pt_previous] == 'C') {
            sprintf (buf, "%s\n%s", "The first point of the y-axial line is the cross point", "Please check the Axial Data file");
            PrintError (buf);
          }

          // y-축선의 시작점의 아래쪽의 점을 자기자신으로 저장
          this->EWNS_index[pt_previous][3] = pt_previous;

          // y-축선의 시작점의 x-좌표를 저장
          this->yaxial[jy][0] = pts[pt_previous][0];

          // y-축선의 시작점의 y-좌표를 저장
          this->yaxial[jy][1] = pts[pt_previous][1];

          // y-축선의 시작점을 포함하는 y-축선의 index를 저장
          this->axial_index[pt_previous][1] = jy;

          // y-축선의 두번재 점의 index를 읽어서 현재점의 index에 저장
          AxialFile >> pt_present;

          // yaxial_index에 점의 index를 저장
          this->yaxial_index[Ny] = pt_present;

          // 읽어들인 점의 개수를 추가
          Ny += 1;

          // 읽어들인 점이 경계점인 동안 반복
          while (this->bc_u[pt_present] == 'C') {
            // 이전에 읽어들인 점의 위쪽점의 index에 현재 읽어들인 점의 index를 저장
            this->EWNS_index[pt_previous][2] = pt_present;

            // 현재 읽어들인 점의 아래쪽점의 index에 이전에 읽어들인 점의 index를 저장
            this->EWNS_index[pt_present][3] = pt_previous;

            // 현재점을 포함하는 y-축선의 index를 저장
            this->axial_index[pt_present][1] = jy;

            // 이전에 읽어들인 점에 현재점의 index를 저장
            pt_previous = pt_present;

            // y-축선의 다음점의 index를 읽어서 현재점의 index에 저장
            AxialFile >> pt_present;

            // yaxial_index에 점의 index를 저장
            this->yaxial_index[Ny] = pt_present;

            // 읽어들인 점의 개수를 추가
            Ny += 1;
          }

          // y-축선의 끝점의 위쪽의 점을 자기자신으로 저장
          this->EWNS_index[pt_present][2] = pt_present;

          // y-축선의 끝점의 y-좌표를 저장
          this->yaxial[jy][2] = pts[pt_present][1];

          // y-축선의 끝점을 포함하는 y-축선의 index를 저장
          this->axial_index[pt_present][1] = jy;

          // 이전에 읽어들인 점의 위쪽점의 index에 y-축선의 끝점의 index를 저장
          this->EWNS_index[pt_previous][2] = pt_present;

          // y-축선의 끝점의 아래쪽점의 index에 이전에 읽어들인 점의 index를 저장
          this->EWNS_index[pt_present][3] = pt_previous;

          // y-축선의 index에 1을 더한다
          jy += 1;
        }

        break;
      }
    }
  }
  /* **** TEMP **** */
  // for (size_t i = 0; i < this->pts_num; i++)
  // for (size_t j = 0; j < 2; j++)
  // if (IsEqualDouble (fabs (pts[i][j]), 3.0E0)) this->normal[i][j] = pts[i][j] / fabs (pts[i][j]), this->normal[i][1 - j] = ZeroValue;

  for (size_t i = 0; i < this->pts_num; i++) {
    double x = pts[i][0], y = pts[i][1];
    // if (y > -1.2E0 && y < 1.2E0 && this->bc_u[i] == 'N') this->bc_u[i] = 'D';

    // if (fabs (x) < 1.6E0 && fabs (y) < 1.6E0) this->normal[i][0] = -this->normal[i][0], this->normal[i][1] = -this->normal[i][1];
    // if (IsEqualDouble (x,  1.5E0) && IsEqualDouble (fabs (y), 5.0E-1)) this->normal [i][0] = -1.0E0, this->normal[i][1] =  ZeroValue;
    // if (IsEqualDouble (x, -1.5E0) && IsEqualDouble (fabs (y), 5.0E-1)) this->normal [i][0] =  1.0E0, this->normal[i][1] =  ZeroValue;
    // if (IsEqualDouble (y,  1.5E0) && IsEqualDouble (fabs (x), 5.0E-1)) this->normal [i][1] = -1.0E0, this->normal[i][0] =  ZeroValue;
    // if (IsEqualDouble (y, -1.5E0) && IsEqualDouble (fabs (x), 5.0E-1)) this->normal [i][1] =  1.0E0, this->normal[i][0] =  ZeroValue;
    //
    // if (IsEqualDouble (x,  1.5E0) && IsEqualDouble (fabs (y), 5.0E-1)) this->normal [i][0] = -1.0E0, this->normal[i][1] =  ZeroValue;
    // if (IsEqualDouble (x, -1.5E0) && IsEqualDouble (fabs (y), 5.0E-1)) this->normal [i][0] =  1.0E0, this->normal[i][1] =  ZeroValue;
    // if (IsEqualDouble (y,  1.5E0) && IsEqualDouble (fabs (x), 5.0E-1)) this->normal [i][1] = -1.0E0, this->normal[i][0] =  ZeroValue;
    // if (IsEqualDouble (y, -1.5E0) && IsEqualDouble (fabs (x), 5.0E-1)) this->normal [i][1] =  1.0E0, this->normal[i][0] =  ZeroValue;
    //
    // if (IsEqualDouble (x,  3.0E0) && IsEqualDouble (fabs (y), 1.0E0)) this->normal [i][0] =  1.0E0, this->normal[i][1] =  ZeroValue;
    // if (IsEqualDouble (x, -3.0E0) && IsEqualDouble (fabs (y), 1.0E0)) this->normal [i][0] = -1.0E0, this->normal[i][1] =  ZeroValue;
    // if (IsEqualDouble (y,  3.0E0) && IsEqualDouble (fabs (x), 1.0E0)) this->normal [i][1] =  1.0E0, this->normal[i][0] =  ZeroValue;
    // if (IsEqualDouble (y, -3.0E0) && IsEqualDouble (fabs (x), 1.0E0)) this->normal [i][1] = -1.0E0, this->normal[i][0] =  ZeroValue;
    //
    // if (IsEqualDouble (x,  3.0E0) && IsEqualDouble (fabs (y), 1.0E0)) this->normal [i][0] =  1.0E0, this->normal[i][1] =  ZeroValue;
    // if (IsEqualDouble (x, -3.0E0) && IsEqualDouble (fabs (y), 1.0E0)) this->normal [i][0] = -1.0E0, this->normal[i][1] =  ZeroValue;
    // if (IsEqualDouble (y,  3.0E0) && IsEqualDouble (fabs (x), 1.0E0)) this->normal [i][1] =  1.0E0, this->normal[i][0] =  ZeroValue;
    // if (IsEqualDouble (y, -3.0E0) && IsEqualDouble (fabs (x), 1.0E0)) this->normal [i][1] = -1.0E0, this->normal[i][0] =  ZeroValue;


    this->normal[i][0] = x / sqrt (x * x + y * y), normal[i][1] = y / sqrt (x * x + y * y);
    //
    if (x * x + y * y < 1.1E0) this->normal[i][0] = -this->normal[i][0], this->normal[i][1] = -this->normal[i][1];
    if (Coordinate (x, y) == Coordinate ( 1.0E0, ZeroValue)) this->normal [i][0] = -1.0E0, this->normal[i][1] =  ZeroValue;
    if (Coordinate (x, y) == Coordinate (-1.0E0, ZeroValue)) this->normal [i][0] =  1.0E0, this->normal[i][1] =  ZeroValue;
    if (Coordinate (x, y) == Coordinate (ZeroValue,  1.0E0)) this->normal [i][0] =  ZeroValue, this->normal[i][1] = -1.0E0;
    if (Coordinate (x, y) == Coordinate (ZeroValue, -1.0E0)) this->normal [i][0] =  ZeroValue, this->normal[i][1] =  1.0E0;
    if (Coordinate (x, y) == Coordinate ( 2.0E0, ZeroValue)) this->normal [i][0] =  1.0E0, this->normal[i][1] =  ZeroValue;
    if (Coordinate (x, y) == Coordinate (-2.0E0, ZeroValue)) this->normal [i][0] = -1.0E0, this->normal[i][1] =  ZeroValue;
    if (Coordinate (x, y) == Coordinate (ZeroValue,  2.0E0)) this->normal [i][0] =  ZeroValue, this->normal[i][1] =  1.0E0;
    if (Coordinate (x, y) == Coordinate (ZeroValue, -2.0E0)) this->normal [i][0] =  ZeroValue, this->normal[i][1] = -1.0E0;
  }


  // for (size_t i = 0; i < this->pts_num; i++) if (IsEqualDouble (pts[i][1], 3.0E0)) this->bc_u[i] = 'D';
  // for (size_t i = 0; i < this->pts_num; i++) if (IsEqualDouble (fabs (pts[i][0]), 3.0E0) && IsEqualDouble (pts[i][1], -1.0E0)) this->bc_u[i] = 'D';

  /* **** TEMP **** */


  // 내부점의 개수를 0으로 초기화
  in_pts_num = 0;

  // 모든점의 개수만큼 반복
  for (size_t i = 0; i < this->pts_num; i++) {
    // 경계조건이 interface인 점들도 내부점으로 취급
    if (bc_u[i] == 'C' || bc_u[i] == 'I' || bc_u[i] == 'N' || bc_u[i] == 'D') {
      in_pts_num += 1;
    }
  }

  // x-축선의 첫번째점의 xaxial_index의 위치를 저장하는 변수의 마지막의 원소의 개수를 저장한다
  xxaxial_index[xxaxial_num] = xaxial_num ;

  // y-축선의 첫번째점의 yaxial_index의 위치를 저장하는 변수의 마지막의 원소의 개수를 저장한다
  yyaxial_index[yyaxial_num] = yaxial_num ;

  // 열었던 Axial Data파일을 닫는다
  AxialFile.close ();

  return *this;
}

/* 경계조건과 경계조건의 값을 지징해서 주기위한 모듈 */
AxialData & AxialData::AssignBoundaryValue () {

  // double Radius_Singularity = 0.2;

  // for (size_t i = 0; i < pts_num; i++) {
  //   // 경계점이 Dirichlet boundary condition인 경우
  //   if (this->bc_u[i] == 'D') this->b_u[i] = b_u_ftn (this->pts[i][0], this->pts[i][1]);
  //   // 경계점이 Neumann boundary condition인 경우
  //   if (this->bc_u[i] == 'N') this->b_u[i] = dudx_ftn (this->pts[i][0], this->pts[i][1]) * this->Normal(i, 'x') + dudy_ftn (this->pts[i][0], this->pts[i][1]) * this->Normal(i, 'y');
  // }
  // if (IsEqualDouble (pts[i][0] * pts[i][0] + pts[i][1] * pts[i][1], Radius_Singularity * Radius_Singularity)) bc_u[i] = 'T';

  // if (IsEqualDouble (pts[i][0], 0.0) && fabs (pts[i][1]) < Radius_Singularity - NearZero && bc_u[i] == 'C') bc_u[i] = 'T';
  // if (IsEqualDouble (pts[i][1], 0.0) && fabs (pts[i][0]) < Radius_Singularity - NearZero && bc_u[i] == 'C') bc_u[i] = 'T';
  //
  // if (pts[i][0] * pts[i][0] + pts[i][1] * pts[i][1] < Radius_Singularity * Radius_Singularity - NearZero && bc_u[i] == 'C' && pts[i][0] < - NearZero && pts[i][1] >   NearZero) bc_u[i] = 'T';
  // if (pts[i][0] * pts[i][0] + pts[i][1] * pts[i][1] < Radius_Singularity * Radius_Singularity - NearZero && bc_u[i] == 'C' && pts[i][0] >   NearZero && pts[i][1] >   NearZero) bc_u[i] = 'T';
  // if (pts[i][0] * pts[i][0] + pts[i][1] * pts[i][1] < Radius_Singularity * Radius_Singularity - NearZero && bc_u[i] == 'C' && pts[i][0] < - NearZero && pts[i][1] < - NearZero) bc_u[i] = 'T';
  //
  // if (IsEqualDouble (pts[i][0], 0.0) && IsEqualDouble (pts[i][1], 0.0)) bc_u[i] = 'S';
  // }

  return *this;
}

#endif
