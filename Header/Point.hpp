#ifndef POINT_H
#define POINT_H

#include "Coordinate.hpp"

/* 생성자 */
Point::Point () {
  // coord
  this->coord = (Coordinate (2));
  // minimun coord
  this->min_coord = Coordinate (2);
  // maximum coord
  this->max_coord = Coordinate (2);
  // normal vector
  this->normal = Coordinate (2);
  // mark
  this->mark = "p";
  // solution
  this->u = NULL;
  // derivative
  this->ux = NULL; this->uy = NULL;
  // integration
  this->intg = NULL;
  // 점의 index
  this->index = -1;
  // 점의 행렬에서의 index
  this->mtrxindex = -1;
  // phi
  this->phi = NULL;
  // right-hand side f
  this->f = -1;
  // 오른쪽과 오른쪽의 위, 오른쪽의 아래의 점의 주소
  this->E = NULL; this->EN = NULL; this->ES = NULL;
  // 왼쪽과 왼쪽의 위, 왼쪽의 아래의 점의 주소
  this->W = NULL; this->WN = NULL; this->WS = NULL;
  // 위쪽과 위쪽의 오른쪽, 위쪽의 왼쪽의 점의 주소
  this->N = NULL; this->NE = NULL; this->NW = NULL;
  // 아래쪽과 아래쪽의 오른쪽, 아래쪽의 왼쪽의 점의 주소
  this->S = NULL; this->SE = NULL; this->SW = NULL;
  // 두번째 오른쪽과 두번째 오른쪽의 위, 두번째 오른쪽의 아래의 점의 주소
  this->e = NULL; this->en = NULL; this->es = NULL;
  // 두번째 왼쪽과 두번째 왼쪽의 위, 두번째 왼쪽의 아래의 점의 주소
  this->w = NULL; this->wn = NULL; this->ws = NULL;
  // 두번째 위쪽과 두번째 위쪽의 오른쪽, 두번째 위쪽의 왼쪽의 점의 주소
  this->n = NULL; this->ne = NULL; this->nw = NULL;
  // 두번째 아래쪽과 두번째 아래쪽의 오른쪽, 두번째 아래쪽의 왼쪽의 점의 주소
  this->s = NULL; this->se = NULL; this->sw = NULL;
}

/* 소멸자 */
Point::~Point () {

}

/* 점을 포함하는 국소축선의 시작점과 끝점의 좌표 */
double Point::MinMaxCoordinate (char xy, char mp) {
  // x-축선의 경우
  if (xy == 'X' || xy == 'x') {
    // 시작점의 x-좌표를 return
    if (mp == 'M' || mp == 'm') return this->MinCoord ().Value ('x');
    // 끝점의 x-좌표를 return
    if (mp == 'P' || mp == 'p') return this->MaxCoord ().Value ('x');
  }
  // y-축선의 경우
  if (xy == 'Y' || xy == 'y') {
    // 시작점의 y-좌표를 return
    if (mp == 'M' || mp == 'm') return this->MinCoord ().Value ('y');
    // 끝점의 y-좌표를 return
    if (mp == 'P' || mp == 'p') return this->MaxCoord ().Value ('y');
  }
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("Point::MinMaxCoordinate");
  exit (1);
}

/* 미분의 값 */
Point * Point::Diff (char xy) {
  // x-성분으로의 미분의 값을 return
  if (xy == 'X' || xy == 'x') return this->ux;
  // y-성분으로의 미분의 값을 return
  if (xy == 'Y' || xy == 'y') return this->uy;
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("Point::Diff");
  exit (1);
}

/* conductivity */
double Point::MaterialProperty () {
  // conductivity를 return
  return this->mp_u;
}

/* velocity point */
Point * Point::Velocity (char uv) {
  // u-velocity를 return 하는 경우
  if (uv == 'u' || uv == 'U') return this->uVelocity;
  // v-velocity를 return 하는 경우
  if (uv == 'v' || uv == 'V') return this->vVelocity;
  // 잘못된 참조를 한 경우 에러메시지를 출력하고 종료
  PrintError ("Point::Velocity");
  exit (1);
}

/* 오른쪽, 왼쪽, 위쪽, 아래쪽의 점의 주소 */
Point * Point::EWNS (char EWNS, char ewns) {
  char errorMassage[256];
  // 오른쪽에 있는 점의 참조
  if (EWNS == 'E' || EWNS == 'e') {
    // 오른쪽점의 주소를 return
    if (ewns == 'E' || ewns == 'e') return this->E;
    // 잘못된 잠초, 에러미세지를 출력하고 종료
    if (ewns == 'W' || ewns == 'w') PrintError ("Point::EWNS, EW");
    // 오른쪽의 위쪽점의 주소를 return
    if (ewns == 'N' || ewns == 'n') return this->EN;
    // 오른쪽의 아래쪽의 주소를 return
    if (ewns == 'S' || ewns == 's') return this->ES;
  }
  // 왼쪽에 있는 점의 참조
  if (EWNS == 'W' || EWNS == 'w') {
    // 잘못된 잠초, 에러미세지를 출력하고 종료
    if (ewns == 'E' || ewns == 'e') PrintError ("Point::EWNS, WE");
    // 왼쪽점의 주소를 return
    if (ewns == 'W' || ewns == 'w') return this->W;
    // 왼쪽의 위쪽점의 주소를 return
    if (ewns == 'N' || ewns == 'n') return this->WN;
    // 왼쪽의 아래쪽점의 주소를 return
    if (ewns == 'S' || ewns == 's') return this->WS;
  }
  // 위쪽에 있는 점의 참조
  if (EWNS == 'N' || EWNS == 'n') {
    // 위쪽의 오른쪽점의 주소를 return
    if (ewns == 'E' || ewns == 'e') return this->NE;
    // 위쪽의 왼쪽점의 주소를 return
    if (ewns == 'W' || ewns == 'w') return this->NW;
    // 위쪽점의 주소를 return
    if (ewns == 'N' || ewns == 'n') return this->N;
    // 잘못된 잠초, 에러미세지를 출력하고 종료
    if (ewns == 'S' || ewns == 's') PrintError ("Point::EWNS, NS");
  }
  // 아래쪽에 있는 점의 참조
  if (EWNS == 'S' || EWNS == 's') {
    // 아래쪽의 오른쪽점의 주소를 return
    if (ewns == 'E' || ewns == 'e') return this->SE;
    // 아래쪽의 왼쪽점의 주소를 return
    if (ewns == 'W' || ewns == 'w') return this->SW;
    // 잘못된 잠초, 에러미세지를 출력하고 종료
    if (ewns == 'N' || ewns == 'n') PrintError ("Point::EWNS, SN");
    // 아래쪽점의 주소를 return
    if (ewns == 'S' || ewns == 's') return this->S;
  }
  // 오른쪽, 왼쪽, 위쪽, 아래쪽의 참조가 아닌 경우 에러메시지를 출력하고 종료
  sprintf (errorMassage, "Point::EWNS, EWNS = %c, ewns = %c", EWNS, ewns);
  PrintError (errorMassage);
  exit (1);
}

/* 오른쪽, 왼쪽, 위쪽, 아래쪽의 두번째점의 주소 */
Point * Point::EWNS2nd (char EWNS, char ewns) {
  // 두번째 오른쪽에 있는 점의 참조
  if (EWNS == 'E' || EWNS == 'e') {
    // 두번째 오른쪽점의 주소를 return
    if (ewns == 'E' || ewns == 'e') return this->e;
    // 잘못된 잠초, 에러미세지를 출력하고 종료
    if (ewns == 'W' || ewns == 'w') PrintError ("Point::EWNS2nd, EW");
    // 두번째 오른쪽의 위쪽점의 주소를 return
    if (ewns == 'N' || ewns == 'n') return this->en;
    // 두번째 오른쪽의 아래쪽의 주소를 return
    if (ewns == 'S' || ewns == 's') return this->es;
  }
  // 두번째 왼쪽에 있는 점의 참조
  if (EWNS == 'W' || EWNS == 'w') {
    // 잘못된 잠초, 에러미세지를 출력하고 종료
    if (ewns == 'E' || ewns == 'e') PrintError ("Point::EWNS2nd, WE");
    // 두번째 왼쪽점의 주소를 return
    if (ewns == 'W' || ewns == 'w') return this->w;
    // 두번째 왼쪽의 위쪽점의 주소를 return
    if (ewns == 'N' || ewns == 'n') return this->wn;
    // 두번째 왼쪽의 아래쪽점의 주소를 return
    if (ewns == 'S' || ewns == 's') return this->ws;
  }
  // 두번째 위쪽에 있는 점의 참조
  if (EWNS == 'N' || EWNS == 'n') {
    // 두번째 위쪽의 오른쪽점의 주소를 return
    if (ewns == 'E' || ewns == 'e') return this->ne;
    // 두번째 위쪽의 왼쪽점의 주소를 return
    if (ewns == 'W' || ewns == 'w') return this->nw;
    // 두번째 위쪽점의 주소를 return
    if (ewns == 'N' || ewns == 'n') return this->n;
    // 잘못된 잠초, 에러미세지를 출력하고 종료
    if (ewns == 'S' || ewns == 's') PrintError ("Point::EWNS2nd, NS");
  }
  // 두번째 아래쪽에 있는 점의 참조
  if (EWNS == 'S' || EWNS == 's') {
    // 두번째 아래쪽의 오른쪽점의 주소를 return
    if (ewns == 'E' || ewns == 'e') return this->se;
    // 두번째 아래쪽의 왼쪽점의 주소를 return
    if (ewns == 'W' || ewns == 'w') return this->sw;
    // 잘못된 잠초, 에러미세지를 출력하고 종료
    if (ewns == 'N' || ewns == 'n') PrintError ("Point::EWNS2nd, SN");
    // 두번째 아래쪽점의 주소를 return
    if (ewns == 'S' || ewns == 's') return this->s;
  }
  // 오른쪽, 왼쪽, 위쪽, 아래쪽의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("Point::EWNS2nd");
  exit (1);
}

/* 점의 좌표를 저장 */
Point & Point::SetCoordinate (char xy, double value) {
  // x좌표를 저장
  if (xy == 'x' || xy == 'X') {this->Coord ().SetCoordinate ('x', value); return *this;}
  // y좌표를 저장
  if (xy == 'y' || xy == 'Y') {this->Coord ().SetCoordinate ('y', value); return *this;}
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("Point::SetCoordinate");
  return *this;
}

Point & Point::SetCoordinate (double x, double y) {
  this->Coord ().SetCoordinate (x, y);
  return *this;
}

/* 점을 포함하는 국소축선의 시작점과 끝점의 좌표를 저장 */
Point & Point::SetMinMaxCoordinate (char xy, char mp,  double value) {
  // x-축선의 경우
  if (xy == 'X' || xy == 'x') {
    // 시작점의 x-좌표를 저장
    if (mp == 'M' || mp == 'm') {this->MinCoord ().SetCoordinate ('x', value); return *this;}
    // 끝점의 x-좌표를 저장
    if (mp == 'P' || mp == 'p') {this->MaxCoord ().SetCoordinate ('x', value); return *this;}
  }
  // y-축선의 경우
  if (xy == 'Y' || xy == 'y') {
    // 시작점의 y-좌표를 저장
    if (mp == 'M' || mp == 'm') {this->MinCoord ().SetCoordinate ('y', value); return *this;}
    // 끝점의 y-좌표를 저장
    if (mp == 'P' || mp == 'p') {this->MaxCoord ().SetCoordinate ('y', value); return *this;}
  }
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("Point::SetMinMaxCoordinate");
  exit (1);
}

/* 점을 포함하는 축선의 index */
int Point::Axis (char xy) {
  // x-축선의 경우
  if (xy == 'X' || xy == 'x') return this->axis[0];
  // y-축선의 경우
  if (xy == 'Y' || xy == 'y') return this->axis[1];
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("Point::Axis");
  exit (1);
}

/* 점의 index를 저장 */
Point & Point::SetIndex (int src) {
  // index 저장
  this->index = src;
  // 행렬에서의 index를 -1로 초기화
  this->mtrxindex = -1;
  // 오른쪽과 오른쪽의 위, 오른쪽의 아래의 점의 주소
  this->E = NULL; this->EN = NULL; this->ES = NULL;
  // 왼쪽과 왼쪽의 위, 왼쪽의 아래의 점의 주소
  this->W = NULL; this->WN = NULL; this->WS = NULL;
  // 위쪽과 위쪽의 오른쪽, 위쪽의 왼쪽의 점의 주소
  this->N = NULL; this->NE = NULL; this->NW = NULL;
  // 아래쪽과 아래쪽의 오른쪽, 아래쪽의 왼쪽의 점의 주소
  this->S = NULL; this->SE = NULL; this->SW = NULL;
  // 두번째 오른쪽과 두번째 오른쪽의 위, 두번째 오른쪽의 아래의 점의 주소
  this->e = NULL; this->en = NULL; this->es = NULL;
  // 두번째 왼쪽과 두번째 왼쪽의 위, 두번째 왼쪽의 아래의 점의 주소
  this->w = NULL; this->wn = NULL; this->ws = NULL;
  // 두번째 위쪽과 두번째 위쪽의 오른쪽, 두번째 위쪽의 왼쪽의 점의 주소
  this->n = NULL; this->ne = NULL; this->nw = NULL;
  // 두번째 아래쪽과 두번째 아래쪽의 오른쪽, 두번째 아래쪽의 왼쪽의 점의 주소
  this->s = NULL; this->se = NULL; this->sw = NULL;
  return *this;
}

/* 행렬에서의 index를 저장 */
Point & Point::SetMtrx_Index (int src) {
  // 행렬에서의 index 저장
  this->mtrxindex = src;
  return *this;
}

/* 경계조건을 저장 */
Point & Point::SetCondition (char src) {
  // 경계조건 저장
  this->condition = src;
  return *this;
}

/* 경계조건의 값을 저장 */
Point & Point::SetBoundaryvalue (double src) {
  // 경계조건의 값 저장
  this->b_u = src;
  return *this;
}

/* 점을 포함하는 국소축선의 시작점과 끝점의 좌표를 저장하는 모듈 */
Point & Point::SetMinMax () {
  // 오른쪽점이 존재하는 경우, 그 점의 좌표를 x-축의 끝점의 x-좌표로 저장
  if (this->E != NULL) this->SetMinMaxCoordinate ('x', 'p', this->E->Coord ().Value ('x'));
  // 오른쪽점이 존재하지 않는 경우, 오른쪽의 위쪽점과 오른쪽의 아래쪽점을 이용
  if (this->E == NULL) {
    // 오른쪽의 위쪽점이 존재하는 경우, 그 점의 좌표를 x-축의 끝점의 x-좌표로 저장
    if (this->EN != NULL) this->SetMinMaxCoordinate ('x', 'p', this->EN->Coord ().Value ('x'));
    // 오른쪽의 아래쪽점이 존재하는 경우, 그 점의 좌표를 x-축의 끝점의 x-좌표로 저장
    if (this->ES != NULL) this->SetMinMaxCoordinate ('x', 'p', this->ES->Coord ().Value ('x'));
    // 두 점이 모두 존재하지 않는 경우, 각 점의 index와 에러메시지를 출력하고 종료
    if (this->EN == NULL && this->ES == NULL) {
      printf ("%s%d%s%d%s%d%s%d\n", "pts = ", this->index, " E = ", this->E->Index (), " EN = ", this->EN->Index (), " ES = ", this->ES->Index ());
      PrintError ("Point::SetMinMax");
    }
  }
  // 왼쪽점이 존재하는 경우, 그 점의 좌표를 x-축의 시작점의 x-좌표로 저장
  if (this->W != NULL) this->SetMinMaxCoordinate ('x', 'm', this->W->Coord ().Value ('x'));
  // 왼쪽점이 존재하지 않는 경우, 왼쪽의 위쪽점과 왼쪽의 아래쪽점을 이용
  if (this->W == NULL) {
    // 왼쪽의 위쪽점이 존재하는 경우, 그 점의 좌표를 x-축의 시작점의 x-좌표로 저장
    if (this->WN != NULL) this->SetMinMaxCoordinate ('x', 'm', this->WN->Coord ().Value ('x'));
    // 왼쪽의 아래쪽점이 존재하는 경우, 그 점의 좌표를 x-축의 시작점의 x-좌표로 저장
    if (this->WS != NULL) this->SetMinMaxCoordinate ('x', 'm', this->WS->Coord ().Value ('x'));
    // 두 점이 모두 존재하지 않는 경우, 각 점의 index와 에러메시지를 출력하고 종료
    if (this->WN == NULL && this->WS == NULL) {
      printf ("%s%d%s%d%s%d%s%d\n", "pts = ", this->index, " W = ", this->W->Index (), " WN = ", this->WN->Index (), " WS = ", this->WS->Index ());
      PrintError ("Point::SetMinMax");
    }
  }
  // 위쪽점이 존재하는 경우, 그 점의 좌표를 y-축의 끝점의 y-좌표로 저장
  if (this->N != NULL) this->SetMinMaxCoordinate ('y', 'p', this->N->Coord ().Value ('y'));
  // 위쪽점이 존재하지 않는 경우, 위쪽의 오른쪽점과 위쪽의 왼쪽점을 이용
  if (this->N == NULL) {
    // 위쪽의 오른쪽점이 존재하는 경우, 그 점의 좌표를 y-축의 끝점의 y-좌표로 저장
    if (this->NE != NULL) this->SetMinMaxCoordinate ('y', 'p', this->NE->Coord ().Value ('y'));
    // 위쪽의 왼쪽점이 존재하는 경우, 그 점의 좌표를 y-축의 끝점의 y-좌표로 저장
    if (this->NW != NULL) this->SetMinMaxCoordinate ('y', 'p', this->NW->Coord ().Value ('y'));
    // 두 점이 모두 존재하지 않는 경우, 각 점의 index와 에러메시지를 출력하고 종료
    if (this->NE == NULL && this->NW == NULL) {
      printf ("%s%d%s%d%s%d%s%d\n", "pts = ", this->index, " N = ", this->N->Index (), " NE = ", this->NE->Index (), " NW = ", this->NW->Index ());
      PrintError ("Point::SetMinMax");
    }
  }
  // 아래쪽점이 존재하는 경우, 그 점의 좌표를 y-축의 시작점의 y-좌표로 저장
  if (this->S != NULL) this->SetMinMaxCoordinate ('y', 'm', this->S->Coord ().Value ('y'));
  // 아래쪽점이 존재하지 않는 경우, 아래쪽의 오른쪽점과 아래쪽의 왼쪽점을 이용
  if (this->S == NULL) {
    // 아래쪽의 오른쪽점이 존재하는 경우, 그 점의 좌표를 y-축의 시작점의 y-좌표로 저장
    if (this->SE != NULL) this->SetMinMaxCoordinate ('y', 'm', this->SE->Coord ().Value ('y'));
    // 아래쪽의 왼쪽점이 존재하는 경우, 그 점의 좌표를 y-축의 시작점의 y-좌표로 저장
    if (this->SW != NULL) this->SetMinMaxCoordinate ('y', 'm', this->SW->Coord ().Value ('y'));
    // 두 점이 모두 존재하지 않는 경우, 각 점의 index와 에러메시지를 출력하고 종료
    if (this->SE == NULL && this->SW == NULL) {
      printf ("%s%d%s%d%s%d%s%d\n", "pts = ", this->index, " S = ", this->S->Index (), " SE = ", this->SE->Index (), " SW = ", this->SW->Index ());
      PrintError ("Point::SetMinMax");
    }
  }
  // 점이 Interface위에 있는 경우 Interface위의 점을 포함하는 국소축선의 시작점과 끝점의 좌표를 저장하는 모듈을 실행
  if (this->Condition () == 'I' || this->Condition () == 'M') this->SetInterfaceMinMax ();
  return *this;
}

/* Interface 위의 점을 포함하는 국소축선의 시작점과 끝점의 좌표를 저장하는 모듈 */
Point & Point::SetInterfaceMinMax () {
  // 오른쪽점이 자기자신인 경우
  if (this->EWNS ('E', 'E') == this) {
    // 위쪽점과 위쪽점의 오른쪽점이 존재하는 경우, 위쪽점의 오른쪽점의 x-좌표를 국소 x-축선의 끝점의 x-좌표로 저장
    if (this->EWNS ('N', 'N') != NULL) if (this->EWNS ('N', 'N')->EWNS ('E', 'E') != NULL) SetMinMaxCoordinate ('x', 'p', this->EWNS ('N', 'N')->EWNS ('E', 'E')->Coord ().Value ('x'));
    // 아래쪽점과 아래쪽점의 오른쪽점이 존재하는 경우, 아래쪽점의 오른쪽점의 x-좌표를 국소 x-축선의 끝점의 x-좌표로 저장
    if (this->EWNS ('S', 'S') != NULL) if (this->EWNS ('S', 'S')->EWNS ('E', 'E') != NULL) SetMinMaxCoordinate ('x', 'p', this->EWNS ('S', 'S')->EWNS ('E', 'E')->Coord ().Value ('x'));
    // 위쪽의 오른쪽점이 존재하는 경우, 위쪽의 오른쪽점의 x-좌표를 국소 x-축선의 끝점의 x-좌표로 저장
    if (this->EWNS ('N', 'E') != NULL) SetMinMaxCoordinate ('x', 'p', this->EWNS ('N', 'E')->Coord ().Value ('x'));
    // 아래쪽의 오른쪽점이 존재하는 경우, 아래쪽의 오른쪽점의 x-좌표를 국소 x-축선의 끝점의 x-좌표로 저장
    if (this->EWNS ('S', 'E') != NULL) SetMinMaxCoordinate ('x', 'p', this->EWNS ('S', 'E')->Coord ().Value ('x'));
  }
  // 왼쪽점이 자기자신인 경우
  if (this->EWNS ('W', 'W') == this) {
    // 위쪽점과 위쪽점의 왼쪽점이 존재하는 경우, 위쪽점의 왼쪽점의 x-좌표를 국소 x-축선의 시작점의 x-좌표로 저장
    if (this->EWNS ('N', 'N') != NULL) if (this->EWNS ('N', 'N')->EWNS ('W', 'W') != NULL) SetMinMaxCoordinate ('x', 'm', this->EWNS ('N', 'N')->EWNS ('W', 'W')->Coord ().Value ('x'));
    // 아래쪽점과 아래쪽점의 왼쪽점이 존재하는 경우, 아래쪽점의 왼쪽점의 x-좌표를 국소 x-축선의 시작점의 x-좌표로 저장
    if (this->EWNS ('S', 'S') != NULL) if (this->EWNS ('S', 'S')->EWNS ('W', 'W') != NULL) SetMinMaxCoordinate ('x', 'm', this->EWNS ('S', 'S')->EWNS ('W', 'W')->Coord ().Value ('x'));
    // 위쪽의 왼쪽점이 존재하는 경우, 위쪽의 왼쪽점의 x-좌표를 국소 x-축선의 시작점의 x-좌표로 저장
    if (this->EWNS ('N', 'W') != NULL) SetMinMaxCoordinate ('x', 'm', this->EWNS ('N', 'W')->Coord ().Value ('x'));
    // 아래쪽의 왼쪽점이 존재하는 경우, 아래쪽의 왼쪽점의 x-좌표를 국소 x-축선의 시작점의 x-좌표로 저장
    if (this->EWNS ('S', 'W') != NULL) SetMinMaxCoordinate ('x', 'm', this->EWNS ('S', 'W')->Coord ().Value ('x'));
  }
  // 위쪽점이 자기자신인 경우
  if (this->EWNS ('N', 'N') == this) {
    // 오른쪽점과 오른쪽점의 위쪽점이 존재하는 경우, 오른쪽점의 위쪽점의 y-좌표를 국소 y-축선의 끝점의 y-좌표로 저장
    if (this->EWNS ('E', 'E') != NULL) if (this->EWNS ('E', 'E')->EWNS ('N', 'N') != NULL) SetMinMaxCoordinate ('y', 'p', this->EWNS ('E', 'E')->EWNS ('N', 'N')->Coord ().Value ('y'));
    // 왼쪽점과 왼쪽점의 위쪽점이 존재하는 경우, 왼쪽점의 위쪽점의 y-좌표를 국소 y-축선의 끝점의 y-좌표로 저장
    if (this->EWNS ('W', 'W') != NULL) if (this->EWNS ('W', 'W')->EWNS ('N', 'N') != NULL) SetMinMaxCoordinate ('y', 'p', this->EWNS ('W', 'W')->EWNS ('N', 'N')->Coord ().Value ('y'));
    // 오른쪽의 위쪽점이 존재하는 경우, 오른쪽의 위쪽점의 y-좌표를 국소 y-축선의 끝점의 y-좌표로 저장
    if (this->EWNS ('E', 'N') != NULL) SetMinMaxCoordinate ('y', 'p', this->EWNS ('E', 'N')->Coord ().Value ('y'));
    // 왼쪽의 위쪽점이 존재하는 경우, 왼쪽의 위쪽점의 y-좌표를 국소 y-축선의 끝점의 y-좌표로 저장
    if (this->EWNS ('W', 'N') != NULL) SetMinMaxCoordinate ('y', 'p', this->EWNS ('W', 'N')->Coord ().Value ('y'));
  }
  // 아래쪽점이 자기자신인 경우
  if (this->EWNS ('S', 'S') == this) {
    // 오른쪽점과 오른쪽점의 아래쪽점이 존재하는 경우, 오른쪽점의 아래쪽점의 y-좌표를 국소 y-축선의 시작점의 y-좌표로 저장
    if (this->EWNS ('E', 'E') != NULL) if (this->EWNS ('E', 'E')->EWNS ('S', 'S') != NULL) SetMinMaxCoordinate ('y', 'm', this->EWNS ('E', 'E')->EWNS ('S', 'S')->Coord ().Value ('y'));
    // 왼쪽점과 왼쪽점의 아래쪽점이 존재하는 경우, 왼쪽점의 아래쪽점의 y-좌표를 국소 y-축선의 시작점의 y-좌표로 저장
    if (this->EWNS ('W', 'W') != NULL) if (this->EWNS ('W', 'W')->EWNS ('S', 'S') != NULL) SetMinMaxCoordinate ('y', 'm', this->EWNS ('W', 'W')->EWNS ('S', 'S')->Coord ().Value ('y'));
    // 오른쪽의 아래쪽점이 존재하는 경우, 오른쪽의 아래쪽점의 y-좌표를 국소 y-축선의 시작점의 y-좌표로 저장
    if (this->EWNS ('E', 'S') != NULL) SetMinMaxCoordinate ('y', 'm', this->EWNS ('E', 'S')->Coord ().Value ('y'));
    // 왼쪽의 아래쪽점이 존재하는 경우, 왼쪽의 아래쪽점의 y-좌표를 국소 y-축선의 시작점의 y-좌표로 저장
    if (this->EWNS ('W', 'S') != NULL) SetMinMaxCoordinate ('y', 'm', this->EWNS ('W', 'S')->Coord ().Value ('y'));
  }
  return *this;
}

/* solution을 저장 */
Point & Point::SetValue (double values) {
  // solution 저장
  this->value = values; return *this;
}

/* 직전의 solution을 저장 */
Point & Point::SetPre (Point *src) {
  // solution 저장
  this->pre = src; return *this;
}

/* 전전의 solution을 저장 */
Point & Point::SetOld (Point *src) {
  // solution 저장
  this->old = src; return *this;
}

/* 중간의 solution을 저장 */
Point & Point::SetHat (Point *src) {
  // solution 저장
  this->hat = src; return *this;
}

/* 미분을 저장 */
Point & Point::SetDiff (char xy, Point * src) {
  // x-성분으로의 미분의 값을 저장
  if (xy == 'X' || xy == 'x') {this->ux = src; return *this;}
  // y-성분으로의 미분의 값을 저장
  if (xy == 'Y' || xy == 'y') {this->uy = src; return *this;}
  // x 혹은 y의 참조가 아닌 경우 에러메시지를 출력하고 종료
  PrintError ("Point::SetDiff");
  exit (1);
}

/* 미분하기 전의 점을 저장 */
Point & Point::SetIntg (Point * src) {
  // 미분하기 전의 점을 저장
  this->intg = src; return *this;
}

/* phi를 저장 */
Point & Point::SetPhi (Point *src) {
  // phi를 저장
  this->phi = src; return *this;
}

/* right-hand side f를 저장 */
Point & Point::SetF (double value) {
  // right-hand side f를 저장
  this->f = value; return *this;
}

/* velocity point를 저장 */
Point & Point::SetVelocity (char uv, Point *source) {
  // u-velocity를 저장
  if (uv == 'U' || uv == 'u') {this->uVelocity = source; return *this;}
  // v-velocity를 저장
  if (uv == 'V' || uv == 'v') {this->vVelocity = source; return *this;}
  // 잘못된 참조의 경우 에러메시지를 출력하고 종료
  PrintError ("Point::SetVelocity");
  exit (1);
}

/* pressure point를 저장 */
Point & Point::SetPressure (Point * source) {
  // 압력의 점의 주소를 저장
  this->pressure = source;
  return *this;
}

/* mark를 저장 */
Point & Point::SetMark (string src) {
  // mark를 저장
  this->mark = src; return *this;
}

/* conductivity를 저장 */
Point & Point::SetMaterialProperty (double value) {
  // conductivity 저장
  this->mp_u = value;
  return *this;
}

/* 종료시각을 저장 */
Point & Point::SetTerminalTime (double value) {
  // 종료시각을 저장
  this->terminalTime = value;
  return *this;
}

/* 현재 시각을 저장  */
Point & Point::SetTime (double value) {
  // 현재 시각을 저장
  this->presentTime = value;
  return *this;
}

/* time step을 저장 */
Point & Point::SetTimeStep (int value) {
  // time step을 저장
  this->timeStep = value;
  return *this;
}

/* dt를 저장 */
Point & Point::SetDt (double value) {
  // dt를 저장
  this->dt = value;
  return *this;
}

/* 현재시각을 update */
bool Point::NextTime () {
  // 현재시각을 update
  this->SetTime (this->Time () + this->TimeStep ());
  // 현재 시각이 종료시각을 넘었으면 거짓값을 return
  if (this->Time () > this->TerminalTime ()) return false;
  // 현재 시각이 종료시각을 넘지 않았으면 참값을 return
  else return true;
}

/* normal vector를 저장 */
Point & Point::SetNormal (char xy, double src) {
  // x-성분의 normal vector를 저장하는 경우
  if (xy == 'X' || xy == 'x') this->Normal ().SetCoordinate ('x', src);
  // y-성분의 normal vector를 저장하는 경우
  if (xy == 'Y' || xy == 'y') this->Normal ().SetCoordinate ('y', src);
  return *this;
}

/* normal vector를 저장 */
Point & Point::SetNormal (Coordinate src) {
  this->Normal () = src;
  return *this;
}

/* 주위점들을 찾는 모듈 */
Point & Point::FindAxialElement (AxialData *adat, Point *pts) {
  // 축선의 index의 초기화
  int axialline = -1;
  // 점을 포함하는 x-축선의 index를 저장
  this->axis[0] = adat->Axial_Index (this->Index (), 'x');
  // 점을 포함하는 y-축선의 index를 저장
  this->axis[1] = adat->Axial_Index (this->Index (), 'y');
  // 오른쪽점이 존재하는 경우, 오른쪽점의 주소를 저장
  if (adat->EWNS_Index (this->Index (), 'E') > -1) this->E = &pts[adat->EWNS_Index (this->Index (), 'E')];
  // 왼쪽점이 존재하는 경우, 왼쪽점의 주소를 저장
  if (adat->EWNS_Index (this->Index (), 'W') > -1) this->W = &pts[adat->EWNS_Index (this->Index (), 'W')];
  // 위쪽점이 존재하는 경우, 위쪽점의 주소를 저장
  if (adat->EWNS_Index (this->Index (), 'N') > -1) this->N = &pts[adat->EWNS_Index (this->Index (), 'N')];
  // 아래쪽점이 존재하는 경우, 아래쪽점의 주소를 저장
  if (adat->EWNS_Index (this->Index (), 'S') > -1) this->S = &pts[adat->EWNS_Index (this->Index (), 'S')];
  // 오른쪽점이 존재하지 않는 경우
  if (this->EWNS ('E', 'E') == NULL) {
    // 가장 가까이에 있는 오른쪽 축선의 index를 찾아서 저장
    axialline = FindEastAxialLine (adat);
    // 찾은 오른쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 위쪽점의 주소를 찾아서 저장
    this->EN = &pts[FindVerticalPoints (adat, axialline, 'N')];
    // 찾은 오른쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 아래쪽점의 주소를 찾아서 저장
    this->ES = &pts[FindVerticalPoints (adat, axialline, 'S')];
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS ('E', 'N') == this->EWNS ('E', 'S')) {
      // 찾은 점을 오른쪽점으로 저장
      this->E = this->EWNS ('E', 'N');
      // 오른쪽의 위쪽점 초기화
      this->EN = NULL;
      // 오른쪽의 아래쪽점 초기화
      this->ES = NULL;
    }
    // 오른쪽의 위쪽점이 존재하지 않고 오른쪽의 아래쪽점이 존재하는 경우, 그 점을 오른쪽점으로 저장
    if (this->EWNS ('E', 'N') == NULL && this->EWNS ('E', 'S') != NULL) this->E = this->EWNS ('E', 'S');
    // 오른쪽의 위쪽점이 존재하고 오른쪽의 아래쪽점이 존재하지 않는 경우, 그 점을 오른쪽점으로 저장
    if (this->EWNS ('E', 'N') != NULL && this->EWNS ('E', 'S') == NULL) this->E = this->EWNS ('E', 'S');
  }
  // 왼쪽점이 존재하지 않는 경우
  if (this->EWNS ('W', 'W') == NULL) {
    // 가장 가까이에 있는 왼쪽 축선의 index를 찾아서 저장
    axialline = FindWestAxialLine (adat);
    // 찾은 왼쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 위쪽점의 주소를 찾아서 저장
    this->WN = &pts[FindVerticalPoints (adat, axialline, 'N')];
    // 찾은 왼쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 아래쪽점의 주소를 찾아서 저장
    this->WS = &pts[FindVerticalPoints (adat, axialline, 'S')];
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS ('W', 'N') == this->EWNS ('W', 'S')) {
      // 찾은 점을 왼쪽점으로 저장
      this->W = this->EWNS ('W', 'N');
      // 왼쪽의 위쪽점 초기화
      this->WN = NULL;
      // 왼쪽의 아래쪽점 초기화
      this->WS = NULL;
    }
    // 왼쪽의 위쪽점이 존재하지 않고 왼쪽의 아래쪽점이 존재하는 경우, 그 점을 왼쪽점으로 저장
    if (this->EWNS ('W', 'N') == NULL && this->EWNS ('W', 'S') != NULL) this->W = this->EWNS ('W', 'S');
    // 왼쪽의 위쪽점이 존재하고 왼쪽의 아래쪽점이 존재하지 않는 경우, 그 점을 왼쪽점으로 저장
    if (this->EWNS ('W', 'N') != NULL && this->EWNS ('W', 'S') == NULL) this->W = this->EWNS ('W', 'N');
  }
  // 위쪽점이 존재하지 않는 경우
  if (this->EWNS ('N', 'N') == NULL) {
    // 가장 가까인에 있는 위쪽 축선의 index를 찾아서 저장
    axialline = FindNorthAxialLine (adat);
    // 찾은 위쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 오른쪽점의 주소를 찾아서 저장
    this->NE = &pts[FindHorizontalPoints (adat, axialline, 'E')];
    // 찾은 위쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 왼쪽점의 주소를 찾아서 저장
    this->NW = &pts[FindHorizontalPoints (adat, axialline, 'W')];
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS ('N', 'E') == this->EWNS ('N', 'W')) {
      // 찾은 점을 위쪽점으로 저장
      this->N = this->EWNS ('N', 'E');
      // 위쪽의 오른쪽점 초기화
      this->NE = NULL;
      // 위쪽의 왼쪽점 초기화
      this->NW = NULL;
    }
    // 위쪽의 오른쪽점이 존재하지 않고 위쪽의 왼쪽점이 존재하는 경우, 그 점을 위쪽점으로 저장
    if (this->EWNS ('N', 'E') == NULL && this->EWNS ('N', 'W') != NULL) this->N = this->EWNS ('N', 'W');
    // 위쪽의 오른쪽점이 존재하고 위쪽의 왼쪽점이 존재하지 않는 경우, 그 점을 위쪽점으로 저장
    if (this->EWNS ('N', 'E') != NULL && this->EWNS ('N', 'W') == NULL) this->N = this->EWNS ('N', 'E');
  }
  // 아래쪽점이 존재하지 않는 경우
  if (this->EWNS ('S', 'S') == NULL) {
    // 가장 가까인에 있는 아래쪽 축선의 index를 찾아서 저장
    axialline = FindSouthAxialLine (adat);
    // 찾은 아래쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 오른쪽점의 주소를 찾아서 저장
    this->SE = &pts[FindHorizontalPoints (adat, axialline, 'E')];
    // 찾은 아래쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 왼쪽점의 주소를 찾아서 저장
    this->SW = &pts[FindHorizontalPoints (adat, axialline, 'W')];
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS ('S', 'E') == this->EWNS ('S', 'W')) {
      // 찾은 점을 아래쪽점으로 저장
      this->S = this->EWNS ('S', 'E');
      // 아래쪽의 오른쪽점 초기화
      this->SE = NULL;
      // 아래쪽의 왼쪽점 초기화
      this->SW = NULL;
    }
    // 아래쪽의 오른쪽점이 존재하지 않고 아래쪽의 왼쪽점이 존재하는 경우, 그 점을 아래쪽점으로 저장
    if (this->EWNS ('S', 'E') == NULL && this->EWNS ('S', 'W') != NULL) this->S = this->EWNS ('S', 'W');
    // 아래쪽의 오른쪽점이 존재하고 아래쪽의 왼쪽점이 존재하지 않는 경우, 그 점을 아래쪽점으로 저장
    if (this->EWNS ('S', 'E') != NULL && this->EWNS ('S', 'W') == NULL) this->S = this->EWNS ('S', 'E');
  }
  return *this;
}

/* 두번째 주위점들을 찾는 모듈 */
Point & Point::Find2ndAxialElement (AxialData *adat, Point *pts) {
  // 축선의 index의 초기화
  int axialline = -1;
  // 오른쪽점이 존재하고 오른쪽점의 오른쪽점이 존재하는 경우, 오른쪽점의 오른쪽점의 주소를 저장
  if (this->EWNS ('E', 'E') != NULL) if (this->EWNS ('E', 'E')->EWNS ('E', 'E') != NULL) this->e = this->EWNS ('E', 'E')->EWNS ('E', 'E');
  // 왼쪽점이 존재하고 왼쪽점의 왼쪽점이 존재하는 경우, 왼쪽점의 왼쪽점의 주소를 저장
  if (this->EWNS ('W', 'W') != NULL) if (this->EWNS ('W', 'W')->EWNS ('W', 'W') != NULL) this->w = this->EWNS ('W', 'W')->EWNS ('W', 'W');
  // 위쪽점이 존재하고 위쪽점의 위쪽점이 존재하는 경우, 위쪽점의 위쪽점의 주소를 저장
  if (this->EWNS ('N', 'N') != NULL) if (this->EWNS ('N', 'N')->EWNS ('N', 'N') != NULL) this->n = this->EWNS ('N', 'N')->EWNS ('N', 'N');
  // 아래쪽점이 존재하고 아래쪽점의 아래쪽점이 존재하는 경우, 아래쪽점의 아래쪽점의 주소를 저장
  if (this->EWNS ('S', 'S') != NULL) if (this->EWNS ('S', 'S')->EWNS ('S', 'S') != NULL) this->s = this->EWNS ('S', 'S')->EWNS ('S', 'S');
  // 오른쪽점의 오른쪽점이 존재하지 않는 경우
  if (this->EWNS2nd ('E', 'E') == NULL) {
    // 두번째로 가까이에 있는 오른쪽 축선의 index를 찾아서 저장
    axialline = FindEast2ndAxialLine (adat);
    // 찾은 오른쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 위쪽점의 주소를 찾아서 저장
    this->en = &pts[FindVerticalPoints (adat, axialline, 'N')];
    // 찾은 오른쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 아래쪽점의 주소를 찾아서 저장
    this->es = &pts[FindVerticalPoints (adat, axialline, 'S')];
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS2nd ('E', 'N') == this->EWNS2nd ('E', 'S')) {
      // 찾은 점을 두번째 오른쪽점으로 저장
      this->e = this->EWNS2nd ('E', 'N');
      // 두번째 오른쪽의 위쪽점 초기화
      this->en = NULL;
      // 두번째 오른쪽의 아래쪽점 초기화
      this->es = NULL;
    }
    // 두번째 오른쪽의 위쪽점이 존재하지 않고 두번째 오른쪽의 아래쪽점이 존재하는 경우, 그 점을 두번재 오른쪽점으로 저장
    if (this->EWNS2nd ('E', 'N') == NULL && this->EWNS2nd ('E', 'S') != NULL) this->e = this->EWNS2nd ('E', 'S');
    // 두번째 오른쪽의 위쪽점이 존재하고 두번째 오른쪽의 아래쪽점이 존재하지 않는 경우, 그 점을 두번째 오른쪽점으로 저장
    if (this->EWNS2nd ('E', 'N') != NULL && this->EWNS2nd ('E', 'S') == NULL) this->e = this->EWNS2nd ('E', 'S');
  }
  // 왼쪽점의 왼쪽점이 존재하지 않는 경우
  if (this->EWNS2nd ('W', 'W') == NULL) {
    // 두번째로 가까이에 있는 왼쪽 축선의 index를 찾아서 저장
    axialline = FindWest2ndAxialLine (adat);
    // 찾은 왼쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 위쪽점의 주소를 찾아서 저장
    this->wn = &pts[FindVerticalPoints (adat, axialline, 'N')];
    // 찾은 왼쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 아래쪽점의 주소를 찾아서 저장
    this->ws = &pts[FindVerticalPoints (adat, axialline, 'S')];
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS2nd ('W', 'N') == this->EWNS2nd ('W', 'S')) {
      // 찾은 점을 두번째 왼쪽점으로 저장
      this->w = this->EWNS2nd ('W', 'N');
      // 두번째 왼쪽의 위쪽점 초기화
      this->wn = NULL;
      // 두번째 왼쪽의 아래쪽점 초기화
      this->ws = NULL;
    }
    // 두번째 왼쪽의 위쪽점이 존재하지 않고 두번째 왼쪽의 아래쪽점이 존재하는 경우, 그 점을 두번째 왼쪽점으로 저장
    if (this->EWNS2nd ('W', 'N') == NULL && this->EWNS2nd ('W', 'S') != NULL) this->w = this->EWNS2nd ('W', 'S');
    // 두번째 왼쪽의 위쪽점이 존재하고 두번째 왼쪽의 아래쪽점이 존재하지 않는 경우, 그 점을 두번째 왼쪽점으로 저장
    if (this->EWNS2nd ('W', 'N') != NULL && this->EWNS2nd ('W', 'S') == NULL) this->w = this->EWNS2nd ('W', 'N');
  }
  // 위쪽점의 위쪽점이 존재하지 않는 경우
  if (this->EWNS2nd ('N', 'N') == NULL) {
    // 두번째로 가까이에 있는 위쪽 축선의 index를 찾아서 저장
    axialline = FindNorth2ndAxialLine (adat);
    // 찾은 위쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 오른쪽점의 주소를 찾아서 저장
    this->ne = &pts[FindHorizontalPoints (adat, axialline, 'E')];
    // 찾은 위쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 왼쪽점의 주소를 찾아서 저장
    this->nw = &pts[FindHorizontalPoints (adat, axialline, 'W')];
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS2nd ('N', 'E') == this->EWNS2nd ('N', 'W')) {
      // 찾은 점을 두번째 위쪽점으로 저장
      this->n = this->EWNS2nd ('N', 'E');
      // 두번째 위쪽의 오른쪽점 초기화
      this->ne = NULL;
      // 두번째 위쪽의 왼쪽점 초기화
      this->nw = NULL;
    }
    // 두번째 위쪽의 오른쪽점이 존재하지 않고 두번째 위쪽의 왼쪽점이 존재하는 경우, 그 점을 두번째 위쪽점으로 저장
    if (this->EWNS2nd ('N', 'E') == NULL && this->EWNS2nd ('N', 'W') != NULL) this->n = this->EWNS2nd ('N', 'W');
    // 두번째 위쪽의 오른쪽점이 존재하고 두번째 위쪽의 왼쪽점이 존재하지 않는 경우, 그 점을 두번째 위쪽점으로 저장
    if (this->EWNS2nd ('N', 'E') != NULL && this->EWNS2nd ('N', 'W') == NULL) this->n = this->EWNS2nd ('N', 'E');
  }
  // 아래쪽점의 아래쪽점이 존재하지 않는 경우
  if (this->EWNS2nd ('S', 'S') == NULL) {
    // 두번째로 가까이에 있는 아래쪽 축선의 index를 찾아서 저장
    axialline = FindSouth2ndAxialLine (adat);
    // 찾은 아래쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 오른쪽점의 주소를 찾아서 저장
    this->se = &pts[FindHorizontalPoints (adat, axialline, 'E')];
    // 찾은 아래쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 왼쪽점의 주소를 찾아서 저장
    this->sw = &pts[FindHorizontalPoints (adat, axialline, 'W')];
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS2nd ('S', 'E') == this->EWNS2nd ('S', 'W')) {
      // 찾은 점을 아래쪽점으로 저장
      this->s= this->EWNS2nd ('S', 'E');
      // 아래쪽의 오른쪽점 초기화
      this->se = NULL;
      // 아래쪽의 왼쪽점 초기화
      this->sw = NULL;
    }
    // 두번째 아래쪽의 오른쪽점이 존재하지 않고 두번째 아래쪽의 왼쪽점이 존재하는 경우, 그 점을 두번째 아래쪽점으로 저장
    if (this->EWNS2nd ('S', 'E') == NULL && this->EWNS2nd ('S', 'W') != NULL) this->s = this->EWNS2nd ('S', 'W');
    // 두번째 아래쪽의 오른쪽점이 존재하고 두번째 아래쪽의 왼쪽점이 존재하지 않는 경우, 그 점을 두번째 아래쪽점으로 저장
    if (this->EWNS2nd ('S', 'E') != NULL && this->EWNS2nd ('S', 'W') == NULL) this->s = this->EWNS2nd ('S', 'E');
  }
  return *this;
}

Point & Point::FindAxialElementIntp (AxialData *adat, Point *pts) {
  // 축선의 index의 초기화
  int axialline = -1;
  // 축선의 index가 -1인지를 확인하는 switch
  int SwitchE = -1, SwitchW = -1, SwitchN = -1, SwitchS = -1;
  // 저장을 하기 위한 Point변수의 주소
  Point *trg = NULL, *trg1 = NULL, *trg2 = NULL;
  // 오른쪽점이 존재하지 않는 경우
  if (this->EWNS ('E', 'E') == NULL) {
    // 가장 가까이에 있는 오른쪽 축선의 index를 찾아서 저장
    axialline = FindEastAxialLine (adat);
    // 찾은 오른쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 위쪽점의 주소를 찾아서 저장
    this->EN = &pts[FindVerticalPoints (adat, axialline, 'N')];
    // 찾은 오른쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 아래쪽점의 주소를 찾아서 저장
    this->ES = &pts[FindVerticalPoints (adat, axialline, 'S')];
    // 축선의 index가 -1인 경우 switch를 1으로 놓느다
    if (axialline == -1) SwitchE = 1;
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS ('E', 'N') == this->EWNS ('E', 'S')) {
      // 찾은 점을 오른쪽점으로 저장
      this->E = this->EWNS ('E', 'N');
      // 오른쪽의 위쪽점 초기화
      this->EN = NULL;
      // 오른쪽의 아래쪽점 초기화
      this->ES = NULL;
    }
    // 오른쪽의 위쪽점이 존재하지 않고 오른쪽의 아래쪽점이 존재하는 경우, 그 점을 오른쪽점으로 저장
    if (this->EWNS ('E', 'N') == NULL && this->EWNS ('E', 'S') != NULL) this->E = this->EWNS ('E', 'S');
    // 오른쪽의 위쪽점이 존재하고 오른쪽의 아래쪽점이 존재하지 않는 경우, 그 점을 오른쪽점으로 저장
    if (this->EWNS ('E', 'N') != NULL && this->EWNS ('E', 'S') == NULL) this->E = this->EWNS ('E', 'S');
  }
  // 왼쪽점이 존재하지 않는 경우
  if (this->EWNS ('W', 'W') == NULL) {
    // 가장 가까이에 있는 왼쪽 축선의 index를 찾아서 저장
    axialline = FindWestAxialLine (adat);
    // 찾은 왼쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 위쪽점의 주소를 찾아서 저장
    this->WN = &pts[FindVerticalPoints (adat, axialline, 'N')];
    // 찾은 왼쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 아래쪽점의 주소를 찾아서 저장
    this->WS = &pts[FindVerticalPoints (adat, axialline, 'S')];
    // 축선의 index가 -1인 경우 switch를 2으로 놓느다
    if (axialline == -1) SwitchW = 2;
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS ('W', 'N') == this->EWNS ('W', 'S')) {
      // 찾은 점을 왼쪽점으로 저장
      this->W = this->EWNS ('W', 'N');
      // 왼쪽의 위쪽점 초기화
      this->WN = NULL;
      // 왼쪽의 아래쪽점 초기화
      this->WS = NULL;
    }
    // 왼쪽의 위쪽점이 존재하지 않고 왼쪽의 아래쪽점이 존재하는 경우, 그 점을 왼쪽점으로 저장
    if (this->EWNS ('W', 'N') == NULL && this->EWNS ('W', 'S') != NULL) this->W = this->EWNS ('W', 'S');
    // 왼쪽의 위쪽점이 존재하고 왼쪽의 아래쪽점이 존재하지 않는 경우, 그 점을 왼쪽점으로 저장
    if (this->EWNS ('W', 'N') != NULL && this->EWNS ('W', 'S') == NULL) this->W = this->EWNS ('W', 'N');
  }
  // 위쪽점이 존재하지 않는 경우
  if (this->EWNS ('N', 'N') == NULL) {
    // 가장 가까인에 있는 위쪽 축선의 index를 찾아서 저장
    axialline = FindNorthAxialLine (adat);
    // 찾은 위쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 오른쪽점의 주소를 찾아서 저장
    this->NE = &pts[FindHorizontalPoints (adat, axialline, 'E')];
    // 찾은 위쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 왼쪽점의 주소를 찾아서 저장
    this->NW = &pts[FindHorizontalPoints (adat, axialline, 'W')];
    // 축선의 index가 -1인 경우 switch를 3으로 놓느다
    if (axialline == -1) SwitchN = 3;
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS ('N', 'E') == this->EWNS ('N', 'W')) {
      // 찾은 점을 위쪽점으로 저장
      this->N = this->EWNS ('N', 'E');
      // 위쪽의 오른쪽점 초기화
      this->NE = NULL;
      // 위쪽의 왼쪽점 초기화
      this->NW = NULL;
    }
    // 위쪽의 오른쪽점이 존재하지 않고 위쪽의 왼쪽점이 존재하는 경우, 그 점을 위쪽점으로 저장
    if (this->EWNS ('N', 'E') == NULL && this->EWNS ('N', 'W') != NULL) this->N = this->EWNS ('N', 'W');
    // 위쪽의 오른쪽점이 존재하고 위쪽의 왼쪽점이 존재하지 않는 경우, 그 점을 위쪽점으로 저장
    if (this->EWNS ('N', 'E') != NULL && this->EWNS ('N', 'W') == NULL) this->N = this->EWNS ('N', 'E');
  }
  // 아래쪽점이 존재하지 않는 경우
  if (this->EWNS ('S', 'S') == NULL) {
    // 가장 가까인에 있는 아래쪽 축선의 index를 찾아서 저장
    axialline = FindSouthAxialLine (adat);
    // 찾은 아래쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 오른쪽점의 주소를 찾아서 저장
    this->SE = &pts[FindHorizontalPoints (adat, axialline, 'E')];
    // 찾은 아래쪽 축선위에 있는 점들중에서 현재점에서 가장 가까운 왼쪽점의 주소를 찾아서 저장
    this->SW = &pts[FindHorizontalPoints (adat, axialline, 'W')];
    // 축선의 index가 -1인 경우 switch를 4으로 놓느다
    if (axialline == -1) SwitchS = 4;
    // 찾은 두 점이 같은 점인 경우
    if (this->EWNS ('S', 'E') == this->EWNS ('S', 'W')) {
      // 찾은 점을 아래쪽점으로 저장
      this->S = this->EWNS ('S', 'E');
      // 아래쪽의 오른쪽점 초기화
      this->SE = NULL;
      // 아래쪽의 왼쪽점 초기화
      this->SW = NULL;
    }
    // 아래쪽의 오른쪽점이 존재하지 않고 아래쪽의 왼쪽점이 존재하는 경우, 그 점을 아래쪽점으로 저장
    if (this->EWNS ('S', 'E') == NULL && this->EWNS ('S', 'W') != NULL) this->S = this->EWNS ('S', 'W');
    // 아래쪽의 오른쪽점이 존재하고 아래쪽의 왼쪽점이 존재하지 않는 경우, 그 점을 아래쪽점으로 저장
    if (this->EWNS ('S', 'E') != NULL && this->EWNS ('S', 'W') == NULL) this->S = this->EWNS ('S', 'E');
  }

  // 오른쪽의 점이 현재점인 경우
  if (SwitchE == 1) {
    // 왼쪽점이 존재하는 경우
    if (this->EWNS ('W', 'W') != NULL) {
      // 왼쪽의 점을 저장
      trg = this->EWNS ('W', 'W');
      // 저장한 점의 오른쪽점이 존재하는 경우
      while (trg->EWNS ('E', 'E') != NULL) {
        // 저장한 점과 오른쪽의 점이 같은 경우 break
        if (trg == trg->EWNS ('E', 'E')) break;
        // 오른쪽점을 저장
        trg = trg->EWNS ('E', 'E');
        // 저장한 점의 x-좌표가 현재점보다 오른쪽인 경우 break
        if (trg->Coord ().Value ('x') > this->Coord ().Value ('x')) break;
      }
      // 저장한 점의 x-좌표가 현재점보다 왼쪽에 있는 경우 에러메시지를 출력하고 종료
      if (trg->Coord ().Value ('x') < this->Coord ().Value ('x') && trg != trg->EWNS ('E', 'E')) PrintError ("FindAxialElementIntp - SwhitchE");
      // 저장한 점을 오른쪽점으로 저장
      this->E = trg;
    }
    // 왼쪽점의 위쪽점과 왼쪽점의 아래쪽점이 존재하는 경우
    else if (this->EWNS ('W', 'N') != NULL && this->EWNS ('W', 'S') != NULL) {
      // 왼쪽점의 위쪽점과 왼쪽점의 아래쪽점을 저장
      trg1 = this->EWNS ('W', 'N');
      trg2 = this->EWNS ('W', 'S');
      // 저장한 점들의 오른쪽점이 존재한는 경우
      while (trg1->EWNS ('E', 'E') != NULL && trg2->EWNS ('E', 'E') != NULL) {
        // 저장한 점들과 오른쪽의 점들이 같은 경우 break
        if (trg1 == trg1->EWNS ('E', 'E') && trg2 == trg2->EWNS ('E', 'E')) break;
        // 오른쪽의 점들을 저장
        trg1 = trg1->EWNS ('E', 'E');
        trg2 = trg2->EWNS ('E', 'E');
        // 저장한 점들의 x-좌표가 현재점보다 오른쪽인 경우 break
        if (trg1->Coord ().Value ('x') > this->Coord ().Value ('x') && trg2->Coord ().Value ('x') > this->Coord ().Value ('x')) break;
      }
      // 저장한 점들의 x-좌표가 현재점보다 왼쪽에 있는 경우 에러메시지를 출력하고 종료
      if (trg1->Coord ().Value ('x') < this->Coord ().Value ('x') && trg2->Coord ().Value ('x') < this->Coord ().Value ('x') && trg1 != trg1->EWNS ('E', 'E') && trg2 != trg2->EWNS ('E', 'E')) PrintError ("FindAxialElementIntp - SwhitchE");
      // 저장한 점들을 오른쪽의 점들로 저장
      this->E = NULL;
      this->EN = trg1;
      this->ES = trg2;;
    }
  }

  // 왼쪽의 점이 현재점인 경우
  if (SwitchW == 2) {
    // 오른쪽점이 존재하는 경우
    if (this->EWNS ('E', 'E') != NULL) {
      // 오른쪽의 점을 저장
      trg = this->EWNS ('E', 'E');
      // 저장한 점의 왼쪽점이 존재하는 경우
      while (trg->EWNS ('W', 'W') != NULL) {
        // 저장한 점과 왼쪽의 점이 같은 경우 break
        if (trg == trg->EWNS ('W', 'W')) break;
        // 왼쪽점을 저장
        trg = trg->EWNS ('W', 'W');
        // 저장한 점의 x-좌표가 현재점보다 왼쪽인 경우 break
        if (trg->Coord ().Value ('x') < this->Coord ().Value ('x')) break;
      }
      // 저장한 점의 x-좌표가 현재점보다 오른쪽에 있는 경우 에러메시지를 출력하고 종료
      if (trg->Coord ().Value ('x') > this->Coord ().Value ('x') && trg != trg->EWNS ('W', 'W')) PrintError ("FindAxialElementIntp - SwhitchW");
      // 저장한 점을 왼쪽점으로 저장
      this->W = trg;
    }
    // 오른쪽점의 위쪽점과 오른쪽점의 아래쪽점이 존재하는 경우
    else if (this->EWNS ('E', 'N') != NULL && this->EWNS ('E', 'S') != NULL) {
      // 오른쪽점의 위쪽점과 오른쪽점의 아래쪽점을 저장
      trg1 = this->EWNS ('E', 'N');
      trg2 = this->EWNS ('E', 'S');
      // 저장한 점들의 왼쪽점이 존재한는 경우
      while (trg1->EWNS ('W', 'W') != NULL && trg2->EWNS ('W', 'W') != NULL) {
        // 저장한 점들과 왼쪽의 점들이 같은 경우 break
        if (trg1 == trg1->EWNS ('W', 'W') && trg2 == trg2->EWNS ('W', 'W')) break;
        // 왼쪽의 점들을 저장
        trg1 = trg1->EWNS ('W', 'W');
        trg2 = trg2->EWNS ('W', 'W');
        // 저장한 점들의 x-좌표가 현재점보다 왼쪽인 경우 break
        if (trg1->Coord ().Value ('x') < this->Coord ().Value ('x') && trg2->Coord ().Value ('x') < this->Coord ().Value ('x')) break;
      }
      // 저장한 점들의 x-좌표가 현재점보다 오른쪽에 있는 경우 에러메시지를 출력하고 종료
      if (trg1->Coord ().Value ('x') > this->Coord ().Value ('x') && trg2->Coord ().Value ('x') > this->Coord ().Value ('x') && trg1 != trg1->EWNS ('W', 'W') && trg2 != trg2->EWNS ('W', 'W')) PrintError ("FindAxialElementIntp - SwhitchW");
      // 저장한 점들을 왼쪽의 점들로 저장
      this->W = NULL;
      this->WN = trg1;
      this->WS = trg2;;
    }
  }

  // 위쪽의 점이 현재점인 경우
  if (SwitchN == 3) {
    // 아래쪽점이 존재하는 경우
    if (this->EWNS ('S', 'S') != NULL) {
      // 아래쪽의 점을 저장
      trg = this->EWNS ('S', 'S');
      // 저장한 점의 위쪽점이 존재하는 경우
      while (trg->EWNS ('N', 'N') != NULL) {
        // 저장한 점과 위쪽의 점이 같은 경우 break
        if (trg == trg->EWNS ('N', 'N')) break;
        // 위쪽점을 저장
        trg = trg->EWNS ('N', 'N');
        // 저장한 점의 y-좌표가 현재점보다 위쪽인 경우 break
        if (trg->Coord ().Value ('y') > this->Coord ().Value ('y')) break;
      }
      // 저장한 점의 y-좌표가 현재점보다 아래쪽에 있는 경우 에러메시지를 출력하고 종료
      if (trg->Coord ().Value ('y') < this->Coord ().Value ('y') && trg != trg->EWNS ('N', 'N')) PrintError ("FindAxialElementIntp - SwhitchN");
      // 저장한 점을 위쪽점으로 저장
      this->N = trg;
    }
    // 아래쪽점의 오른쪽점과 아래쪽점의 왼쪽점이 존재하는 경우
    else if (this->EWNS ('S', 'E') != NULL && this->EWNS ('S', 'W') != NULL) {
      // 아래쪽점의 오른쪽점과 아래쪽점의 왼쪽점을 저장
      trg1 = this->EWNS ('S', 'E');
      trg2 = this->EWNS ('S', 'W');
      // 저장한 점들의 위쪽점이 존재한는 경우
      while (trg1->EWNS ('N', 'N') != NULL && trg2->EWNS ('N', 'N') != NULL) {
        // 저장한 점들과 위쪽의 점들이 같은 경우 break
        if (trg1 == trg1->EWNS ('N', 'N') && trg2 == trg2->EWNS ('N', 'N')) break;
        // 위쪽의 점들을 저장
        trg1 = trg1->EWNS ('N', 'N');
        trg2 = trg2->EWNS ('N', 'N');
        // 저장한 점들의 y-좌표가 현재점보다 위쪽인 경우 break
        if (trg1->Coord ().Value ('y') > this->Coord ().Value ('y') && trg2->Coord ().Value ('y') > this->Coord ().Value ('y')) break;
      }
      // 저장한 점들의 y-좌표가 현재점보다 아래쪽에 있는 경우 에러메시지를 출력하고 종료
      if (trg1->Coord ().Value ('y') < this->Coord ().Value ('y') && trg2->Coord ().Value ('y') < this->Coord ().Value ('y') && trg1 != trg1->EWNS ('N', 'N') && trg2 != trg2->EWNS ('N', 'N')) PrintError ("FindAxialElementIntp - SwhitchN");
      // 저장한 점들을 위쪽의 점들로 저장
      this->N = NULL;
      this->NE = trg1;
      this->NW = trg2;;
    }
  }

  // 아래쪽의 점이 현재점인 경우
  if (SwitchS == 4) {
    // 위쪽점이 존재하는 경우
    if (this->EWNS ('N', 'N') != NULL) {
      // 위쪽의 점을 저장
      trg = this->EWNS ('N', 'N');
      // 저장한 점의 아래쪽점이 존재하는 경우
      while (trg->EWNS ('S', 'S') != NULL) {
        // 저장한 점과 아래쪽의 점이 같은 경우 break
        if (trg == trg->EWNS ('S', 'S')) break;
        // 아래쪽점을 저장
        trg = trg->EWNS ('S', 'S');
        // 저장한 점의 y-좌표가 현재점보다 아래쪽인 경우 break
        if (trg->Coord ().Value ('y') < this->Coord ().Value ('y')) break;
      }
      // 저장한 점의 y-좌표가 현재점보다 위쪽에 있는 경우 에러메시지를 출력하고 종료
      if (trg->Coord ().Value ('y') > this->Coord ().Value ('y') && trg != trg->EWNS ('S', 'S')) PrintError ("FindAxialElementIntp - SwhitchS");
      // 저장한 점을 아래쪽점으로 저장
      this->S = trg;
    }
    // 위쪽점의 오른쪽점과 위쪽점의 왼쪽점이 존재하는 경우
    else if (this->EWNS ('N', 'E') != NULL && this->EWNS ('N', 'W') != NULL) {
      // 위쪽점의 오른쪽점과 위쪽점의 왼쪽점을 저장
      trg1 = this->EWNS ('N', 'E');
      trg2 = this->EWNS ('N', 'W');
      // 저장한 점들의 아래쪽점이 존재한는 경우
      while (trg1->EWNS ('S', 'S') != NULL && trg2->EWNS ('S', 'S') != NULL) {
        // 저장한 점들과 아래쪽의 점들이 같은 경우 break
        if (trg1 == trg1->EWNS ('S', 'S') && trg2 == trg2->EWNS ('S', 'S')) break;
        // 아래쪽의 점들을 저장
        trg1 = trg1->EWNS ('S', 'S');
        trg2 = trg2->EWNS ('S', 'S');
        // 저장한 점들의 y-좌표가 현재점보다 아래쪽인 경우 break
        if (trg1->Coord ().Value ('y') < this->Coord ().Value ('y') && trg2->Coord ().Value ('y') < this->Coord ().Value ('y')) break;
      }
      // 저장한 점들의 y-좌표가 현재점보다 위쪽에 있는 경우 에러메시지를 출력하고 종료
      if (trg1->Coord ().Value ('y') > this->Coord ().Value ('y') && trg2->Coord ().Value ('y') > this->Coord ().Value ('y') && trg1 != trg1->EWNS ('S', 'S') && trg2 != trg2->EWNS ('S', 'S')) PrintError ("FindAxialElementIntp - SwhitchS");
      // 저장한 점들을 아래쪽의 점들로 저장
      this->S = NULL;
      this->SE = trg1;
      this->SW = trg2;;
    }
  }

  if (this->EWNS ('S', 'S') != NULL) {
    if (this->EWNS ('S', 'S')->EWNS ('N', 'N') != NULL) {
      if (this->EWNS ('S', 'S')->EWNS ('N', 'N')->Coord ().Value ('y') < this->Coord ().Value ('y')) {
        this->S = this->EWNS ('S', 'S')->EWNS ('N', 'N');
      }
    }
  }

  if (this->EWNS ('S', 'E') != NULL) {
    if (this->EWNS ('S', 'E')->EWNS ('N', 'N') != NULL) {
      if (this->EWNS ('S', 'E')->EWNS ('N', 'N')->Coord ().Value ('y') < this->Coord ().Value ('y')) {
        this->S = this->EWNS ('S', 'E')->EWNS ('N', 'N');
      }
    }
  }


  if (this->EWNS ('S', 'W') != NULL) {
    if (this->EWNS ('S', 'W')->EWNS ('N', 'N') != NULL) {
      if (this->EWNS ('S', 'W')->EWNS ('N', 'N')->Coord ().Value ('y') < this->Coord ().Value ('y')) {
        this->S = this->EWNS ('S', 'W')->EWNS ('N', 'N');
      }
    }
  }

  // 오른쪽과 왼쪽의 점이 존재하는 경우
  if (this->EWNS ('E', 'E') != NULL && this->EWNS ('W', 'W') != NULL) {
    // 오른쪽과 왼쪽의 점이 같은 경우
    if (this->EWNS ('E', 'E') == this->EWNS ('W', 'W')) {
      // 오른쪽점이 현재점보다 왼쪽에 있는 경우
      if (this->EWNS ('E', 'E')->Coord ().Value ('x') < this->Coord ().Value ('x')) {
        // 오른쪽점을 없앤다
        this->E = NULL;
        // 오른쪽 위의 점을 위의 오른쪽점으로 한다
        this->EN = this->EWNS ('N', 'E');
        // 오른쪽 아래의 점을 아래의 오른쪽점으로 한다
        this->ES = this->EWNS ('S', 'E');
      }

      // 왼쪽점이 현재점보다 오른쪽에 있는 경우
      if (this->EWNS ('W', 'W')->Coord ().Value ('x') > this->Coord ().Value ('x')) {
        // 왼쪽점을 없앤다
        this->W = NULL;
        // 왼쪽 위의 점을 위의 왼쪽점으로 한다
        this->WN = this->EWNS ('N', 'W');
        // 왼쪽 아래의 점을 아래의 왼쪽점으로 한다
        this->WS = this->EWNS ('S', 'W');
      }
    }
  }

  // 오른쪽과 왼쪽의 위쪽점이 존재하는 경우
  if (this->EWNS ('E', 'N') != NULL && this->EWNS ('W', 'N') != NULL) {
    // 오른쪽과 왼쪽의 위쪽점이 같은 경우
    if (this->EWNS ('E', 'N') == this->EWNS ('W', 'N')) {
      // 오른쪽의 위쪽점이 현재점보다 왼쪽에 있는 경우
      if (this->EWNS ('E', 'N')->Coord ().Value ('x') < this->Coord ().Value ('x')) {
        // 오른쪽 위의 점을 위의 오른쪽점으로 한다
        this->EN = this->EWNS ('N', 'E');
      }

      // 왼쪽의 위쪽점이 현재점보다 오른쪽에 있는 경우
      if (this->EWNS ('W', 'N')->Coord ().Value ('x') > this->Coord ().Value ('x')) {
        // 왼쪽 위의 점을 위의 왼쪽점으로 한다
        this->WN = this->EWNS ('N', 'W');
      }
    }
  }

  // 오른쪽과 왼쪽의 아래쪽점이 존재하는 경우
  if (this->EWNS ('E', 'S') != NULL && this->EWNS ('W', 'S') != NULL) {
    // 오른쪽과 왼쪽의 아래쪽점이 같은 경우
    if (this->EWNS ('E', 'S') == this->EWNS ('W', 'S')) {
      // 오른쪽의 아래쪽점이 현재점보다 왼쪽에 있는 경우
      if (this->EWNS ('E', 'S')->Coord ().Value ('x') < this->Coord ().Value ('x')) {
        // 오른쪽 아래의 점을 위의 오른쪽점으로 한다
        this->ES = this->EWNS ('S', 'E');
      }

      // 왼쪽의 아래쪽점이 현재점보다 오른쪽에 있는 경우
      if (this->EWNS ('W', 'N')->Coord ().Value ('x') > this->Coord ().Value ('x')) {
        // 왼쪽 아래의 점을 위의 왼쪽점으로 한다
        this->WS = this->EWNS ('S', 'W');
      }
    }
  }

  // 위쪽과 아래쪽의 점이 존재하는 경우
  if (this->EWNS ('N', 'N') != NULL && this->EWNS ('S', 'S') != NULL) {
    // 위쪽과 아래쪽의 점이 같은 경우
    if (this->EWNS ('N', 'N') == this->EWNS ('S', 'S')) {
      // 위쪽점이 현재점보다 아래쪽에 있는 경우
      if (this->EWNS ('N', 'N')->Coord ().Value ('y') < this->Coord ().Value ('y')) {
        // 위쪽점을 없앤다
        this->N = NULL;
        // 위쪽의 오른쪽 점을 오른쪽의 위쪽점으로 한다
        this->NE = this->EWNS ('E', 'N');
        // 위쪽의 왼쪽 점을 왼쪽의 위쪽점으로 한다
        this->NW = this->EWNS ('W', 'N');
      }

      // 아래쪽점이 현재점보다 위쪽에 있는 경우
      if (this->EWNS ('S', 'S')->Coord ().Value ('y') > this->Coord ().Value ('y')) {
        // 아래쪽점을 없앤다
        this->S = NULL;
        // 아래쪽의 오른쪽의 점을 오른쪽의 아래쪽점으로 한다
        this->SE = this->EWNS ('E', 'S');
        // 아래쪽의 왼쪽의 점을 왼쪽의 아래쪽점으로 한다
        this->SW = this->EWNS ('W', 'S');
      }
    }
  }

  // 위쪽과 아래쪽의 오른쪽 점이 존재하는 경우
  if (this->EWNS ('N', 'E') != NULL && this->EWNS ('S', 'E') != NULL) {
    // 위쪽과 아래쪽의 오른쪽 점이 같은 경우
    if (this->EWNS ('N', 'E') == this->EWNS ('S', 'E')) {
      // 위쪽의 오른쪽점이 현재점보다 아래쪽에 있는 경우
      if (this->EWNS ('N', 'E')->Coord ().Value ('y') < this->Coord ().Value ('y')) {
        // 위쪽의 오른쪽 점을 오른쪽의 위쪽점으로 한다
        this->NE = this->EWNS ('E', 'N');
      }

      // 아래쪽의 오른쪽점이 현재점보다 위쪽에 있는 경우
      if (this->EWNS ('S', 'E')->Coord ().Value ('y') > this->Coord ().Value ('y')) {
        // 아래쪽의 오른쪽의 점을 오른쪽의 아래쪽점으로 한다
        this->SE = this->EWNS ('E', 'S');
      }
    }
  }

  // 위쪽과 아래쪽의 왼쪽 점이 존재하는 경우
  if (this->EWNS ('N', 'W') != NULL && this->EWNS ('S', 'W') != NULL) {
    // 위쪽과 아래쪽의 왼쪽 점이 같은 경우
    if (this->EWNS ('N', 'W') == this->EWNS ('S', 'W')) {
      // 위쪽의 왼쪽점이 현재점보다 아래쪽에 있는 경우
      if (this->EWNS ('N', 'W')->Coord ().Value ('y') < this->Coord ().Value ('y')) {
        // 위쪽의 왼쪽 점을 오른쪽의 위쪽점으로 한다
        this->NW = this->EWNS ('W', 'N');
      }

      // 아래쪽의 왼쪽점이 현재점보다 위쪽에 있는 경우
      if (this->EWNS ('S', 'W')->Coord ().Value ('y') > this->Coord ().Value ('y')) {
        // 아래쪽의 왼쪽의 점을 왼쪽의 아래쪽점으로 한다
        this->SW = this->EWNS ('W', 'S');
      }
    }
  }
  return *this;
}

/* 경계의 주위점들을 찾는 모듈 */
Point & Point::FindBoundaryElement () {
  // 경계점이 아닌 경우에는 아무것도 하지 않는다
  if (this->Condition ()== 'C' || this->Condition () == 'I' || this->Condition () == 'M' || this->Condition () == 'T') return *this;

  char azimuth[4] = {'E', 'W', 'N', 'S'};
  Point **AZI[3];
  Point **azi[3];

  for (const auto &j : azimuth) {
    if (j == 'E') AZI[0] = &this->E, AZI[1] = &this->EN, AZI[2] = &this->ES, azi[0] = &this->e, azi[1] = &this->en, azi[2] = &this->es;
    if (j == 'W') AZI[0] = &this->W, AZI[1] = &this->WN, AZI[2] = &this->WS, azi[0] = &this->w, azi[1] = &this->wn, azi[2] = &this->ws;
    if (j == 'N') AZI[0] = &this->N, AZI[1] = &this->NE, AZI[2] = &this->NW, azi[0] = &this->n, azi[1] = &this->ne, azi[2] = &this->nw;
    if (j == 'S') AZI[0] = &this->S, AZI[1] = &this->SE, AZI[2] = &this->SW, azi[0] = &this->s, azi[1] = &this->se, azi[2] = &this->sw;

    if (this->IsBoundary (j)) {*AZI[0] = this, *azi[0] = this; for (size_t i = 1; i < 3; i++) *AZI[i] = NULL, *azi[i] = NULL;}
  }
  return *this;
}

/* Interface위의 점인 경우, 양쪽의 conductivity의 동일여부 확인 */
Point & Point::IsInterface () {
  // 현재점의 오른쪽의 conductivity
  double East_conductivity = 0.0;
  // 현재점의 왼쪽의 conductivity
  double West_conductivity = 0.0;
  // 현재점의 위쪽의 conductivity
  double North_conductivity = 0.0;
  // 현재점의 아래쪽의 conductivity
  double South_conductivity = 0.0;
  // Interface위의 점이 아니라면 무시
  if (this->Condition () != 'I') return *this;
  // 오른쪽점이 존재한다면 오른쪽점의 conductivity를 East_conductivity에 저장
  if (this->EWNS ('E', 'E') != NULL) East_conductivity = this->EWNS ('E', 'E')->MaterialProperty ();
  // 오른쪽점이 존재하지 않는다면 경우
  else {
    // 오른쪽의 위쪽점이 존재한다면 오른쪽의 위쪽점의 conductivity를 East_conductivity에 저장
    if (this->EWNS ('E', 'N') != NULL) East_conductivity = this->EWNS ('E', 'N')->MaterialProperty ();
    // 오른쪽의 아래쪽점이 존재한다면 오른쪽의 아래쪽점의 conductivity를 East_conductivity에 저장
    if (this->EWNS ('E', 'S') != NULL) East_conductivity = this->EWNS ('E', 'S')->MaterialProperty ();
  }
  // 왼쪽점이 존재한다면 왼쪽점의 conductivity를 West_conductivity에 저장
  if (this->EWNS ('W', 'W') != NULL) West_conductivity = this->EWNS ('W', 'W')->MaterialProperty ();
  // 왼쪽점이 존재하지 않는다면 경우
  else {
    // 왼쪽의 위쪽점이 존재한다면 왼쪽의 위쪽점의 conductivity를 West_conductivity에 저장
    if (this->EWNS ('W', 'N') != NULL) West_conductivity = this->EWNS ('W', 'N')->MaterialProperty ();
    // 왼쪽의 아래쪽점이 존재한다면 왼쪽의 아래쪽점의 conductivity를 West_conductivity에 저장
    if (this->EWNS ('W', 'S') != NULL) West_conductivity = this->EWNS ('W', 'S')->MaterialProperty ();
  }
  // 위쪽점이 존재한다면 위쪽점의 conductivity를 North_conductivity에 저장
  if (this->EWNS ('N', 'N') != NULL) North_conductivity = this->EWNS ('N', 'N')->MaterialProperty ();
  // 위쪽점이 존재하지 않는다면 경우
  else {
    // 위쪽의 오른쪽점이 존재한다면 위쪽의 오른쪽점의 conductivity를 North_conductivity에 저장
    if (this->EWNS ('N', 'E') != NULL) North_conductivity = this->EWNS ('N', 'E')->MaterialProperty ();
    // 위쪽의 왼쪽점이 존재한다면 위쪽의 왼쪽점의 conductivity를 North_conductivity에 저장
    if (this->EWNS ('N', 'W') != NULL) North_conductivity = this->EWNS ('N', 'W')->MaterialProperty ();
  }
  // 아래쪽점이 존재한다면 위쪽점의 conductivity를 South_conductivity에 저장
  if (this->EWNS ('S', 'S') != NULL) South_conductivity = this->EWNS ('S', 'S')->MaterialProperty ();
  // 아래쪽점이 존재하지 않는다면 경우
  else {
    // 아래쪽의 오른쪽점이 존재한다면 아래쪽의 오른쪽점의 conductivity를 South_conductivity에 저장
    if (this->EWNS ('S', 'E') != NULL) South_conductivity = this->EWNS ('S', 'E')->MaterialProperty ();
    // 아래쪽의 왼쪽점이 존재한다면 아래쪽의 왼쪽점의 conductivity를 South_conductivity에 저장
    if (this->EWNS ('S', 'W') != NULL) South_conductivity = this->EWNS ('S', 'W')->MaterialProperty ();
  }
  // 오른쪽, 왼쪽, 위쪽, 아래쪽의 conductivity를 비교해서 어떤점이 결정
  if (IsEqualDouble (East_conductivity, West_conductivity)) if (IsEqualDouble (North_conductivity, South_conductivity)) this->SetCondition ('M');
  return *this;
}

double Point::CalcRightHandSide () {
  // 오른쪽 변의 값을 저장하기 위한 변수
  double Rf = ZeroValue;

  // Rf =  3.0E0 * this->Velocity ('u')->OldValue () * this->OldDiff ('x') + this->Velocity ('v')->OldValue () * this->OldDiff ('y');
  // Rf -= 3.0E0 * this->Velocity ('u')->PreValue () * this->PreDiff ('x') + this->Velocity ('v')->PreValue () * this->PreDiff ('y');


  // Rf  = 4.0E0 * (this->Velocity ('u')->PreValue () * this->PreDiff ('x') + this->Velocity ('v')->PreValue () * this->PreDiff ('y'));
  // Rf -=          this->Velocity ('u')->OldValue () * this->OldDiff ('x') + this->Velocity ('v')->OldValue () * this->OldDiff ('y');
  // Rf += (3.0E0 * this->PreValue () - this->OldValue ()) / this->TimeStep ();
  // if (this->Mark () == 'u') Rf += this->Pressure ()->PreDiff ('x');
  // if (this->Mark () == 'v') Rf += this->Pressure ()->PreDiff ('y');

  return Rf;
}

/* 가장 가까운 오른쪽의 y-축선의 index를 찾는 모듈 */
int Point::FindEastAxialLine (AxialData *adat) {
  // y-축선의 index를 y-축선의 개수 + 1로 초기화
  int axialIndex = adat->XXYYaxial_Num ('y') + 1;
  // y-축선의 x-좌표를 양의 무한대로 초기화
  double axialCoord = numeric_limits<double>::max ();
  // y-축선의 index를 0부터 y-축선의 개수만큼 반복
  for (size_t i = 0; i < adat->XXYYaxial_Num ('y'); i++) {
    // 반복하는 y-축선의 index가 현재점을 포함하는 y-축선의 index와 같으면 무시
    if (this->Mark ().compare ("p") == 0) if (i == adat->Axial_Index (this->Index (), 'y')) continue;
    // y-축선의 양끝점의 y-좌표 사이에 현재점의 y-좌표가 있는 경우
    if ((adat->XYaxial ('y', i, 1) - this->Coord ().Value ('y')) * (adat->XYaxial ('y', i, 2) - this->Coord ().Value ('y')) < NearZero) {
      // y-축선의 x-좌표가 현재점의 x-좌표보다 큰 경우
      if (adat->XYaxial ('y', i, 0) > this->Coord ().Value ('x') + NearZero) {
        // y-축선의 x-좌표가 이전에 찾은 y-축선의 x-좌표보다 더 작은 경우 (현재점과 더 가까운 경우))
        if (adat->XYaxial ('y', i, 0) < axialCoord) {
          // 찾은 y-축선의 index를 저장
          axialIndex = i;
          // 찾은 y-축선의 x-좌표를 저장
          axialCoord = adat->XYaxial ('y', axialIndex, 0);
        }
      }
    }
  }
  // y-축선을 찾지 못한 경우
  if (axialIndex > adat->XXYYaxial_Num ('y')) return -1;
  // y-축선을 return
  return axialIndex;
}

/*가장 가까운 왼쪽의 y-축선의 index를 찾는 모듈 */
int Point::FindWestAxialLine (AxialData *adat) {
  // y-축선의 index를 y-축선의 개수 + 1로 초기화
  int axialIndex = adat->XXYYaxial_Num ('y') + 1;
  // y-축선의 x-좌표를 음의 무한대로 초기화
  double axialCoord = -numeric_limits<double>::max ();
  // y-축선의 index를 0부터 y-축선의 개수만큼 반복
  for (size_t i = 0; i < adat->XXYYaxial_Num ('y'); i++) {
    // 반복하는 y-축선의 index가 현재점을 포함하는 y-축선의 index와 같으면 무시
    if (this->Mark ().compare ("p") == 0) if (i == adat->Axial_Index (this->Index (), 'y')) continue;
    // y-축선의 양끝점의 y-좌표 사이에 현재점의 y-좌표가 있는 경우
    if ((adat->XYaxial ('y', i, 1) - this->Coord ().Value ('y')) * (adat->XYaxial ('y', i, 2) - this->Coord ().Value ('y')) < NearZero) {
      // y-축선의 x-좌표가 현재점의 x-좌표보다 작은 경우
      if (adat->XYaxial ('y', i, 0) < this->Coord ().Value ('x') - NearZero) {
        // y-축선의 x-좌표가 이전에 찾은 y-축선의 x-좌표보다 더 큰 경우 (현재점과 더 가까운 경우))
        if (adat->XYaxial ('y', i, 0) > axialCoord) {
          // 찾은 y-축선의 index를 저장
          axialIndex = i;
          // 찾은 y-축선의 x-좌표를 저장
          axialCoord = adat->XYaxial ('y', axialIndex, 0);
        }
      }
    }
  }
  // y-축선을 찾지 못한 경우
  if (axialIndex > adat->XXYYaxial_Num ('y')) return -1;
  // y-축선을 return
  return axialIndex;
}

/*가장 가까운 위쪽의 x-축선의 index를 찾는 모듈 */
int Point::FindNorthAxialLine (AxialData *adat) {
  // x-축선의 index를 x-축선의 개수 + 1로 초기화
  int axialIndex = adat->XXYYaxial_Num ('x') + 1;
  // x-축선의 y-좌표를 양의 무한대로 초기화
  double axialCoord = numeric_limits<double>::max ();
  // x-축선의 index를 0부터 x-축선의 개수만큼 반복
  for (size_t i = 0; i < adat->XXYYaxial_Num ('x'); i++) {
    // 반복하는 x-축선의 index가 현재점을 포함하는 x-축선의 index와 같으면 무시
    if (this->Mark ().compare ("p") == 0) if (i == adat->Axial_Index (this->Index (), 'x')) continue;
    // x-축선의 양끝점의 x-좌표 사이에 현재점의 x-좌표가 있는 경우
    if ((adat->XYaxial ('x', i, 1) - this->Coord ().Value ('x')) * (adat->XYaxial ('x', i, 2) - this->Coord ().Value ('x')) < NearZero) {
      // x-축선의 y-좌표가 현재점의 y-좌표보다 큰 경우
      if (adat->XYaxial ('x', i, 0) > this->Coord ().Value ('y') + NearZero) {
        // x-축선의 y-좌표가 이전에 찾은 x-축선의 y-좌표보다 더 작은 경우 (현재점과 더 가까운 경우))
        if (adat->XYaxial ('x', i, 0) < axialCoord) {
          // 찾은 x-축선의 index를 저장
          axialIndex = i;
          // 찾은 x-축선의 y-좌표를 저장
          axialCoord = adat->XYaxial ('x', axialIndex, 0);
        }
      }
    }
  }
  // x-축선을 찾지 못한 경우
  if (axialIndex > adat->XXYYaxial_Num ('x')) return -1;
  // y-축선을 return
  return axialIndex;
}

/*가장 가까운 아래쪽의 x-축선의 index를 찾는 모듈 */
int Point::FindSouthAxialLine (AxialData *adat) {
  // x-축선의 index를 x-축선의 개수 + 1로 초기화
  int axialIndex = adat->XXYYaxial_Num ('x') + 1;
  // x-축선의 y-좌표를 음의 무한대로 초기화
  double axialCoord = -numeric_limits<double>::max ();
  // x-축선의 index를 0부터 x-축선의 개수만큼 반복
  for (size_t i = 0; i < adat->XXYYaxial_Num ('x'); i++) {
    // 반복하는 x-축선의 index가 현재점을 포함하는 x-축선의 index와 같으면 무시
    if (this->Mark ().compare ("p") == 0) if (i == adat->Axial_Index (this->Index (), 'x')) continue;
    // x-축선의 양끝점의 x-좌표 사이에 현재점의 x-좌표가 있는 경우
    if ((adat->XYaxial ('x', i, 1) - this->Coord ().Value ('x')) * (adat->XYaxial ('x', i, 2) - this->Coord ().Value ('x')) < NearZero) {
      // x-축선의 y-좌표가 현재점의 y-좌표보다 작은 경우
      if (adat->XYaxial ('x', i, 0) < this->Coord ().Value ('y') - NearZero) {
        // x-축선의 y-좌표가 이전에 찾은 x-축선의 y-좌표보다 더 큰 경우 (현재점과 더 가까운 경우))
        if (adat->XYaxial ('x', i, 0) > axialCoord) {
          // 찾은 x-축선의 index를 저장
          axialIndex = i;
          // 찾은 x-축선의 y-좌표를 저장
          axialCoord = adat->XYaxial ('x', axialIndex, 0);
        }
      }
    }
  }
  // x-축선을 찾지 못한 경우
  if (axialIndex > adat->XXYYaxial_Num ('x')) return -1;
  // y-축선을 return
  return axialIndex;
}

/* 두번재로 가까운 오른쪽의 y-축선의 index를 찾는 모듈 */
int Point::FindEast2ndAxialLine (AxialData *adat) {
  // 가장 가까운 y-축선의 index를 y-축선의 개수 + 1로 초기화
  int axialIndex1 = adat->XXYYaxial_Num ('y') + 1;
  // 두번째로 가까운 y-축선의 index를 y-축선의 개수 + 1로 초기화
  int axialIndex2 = adat->XXYYaxial_Num ('y') + 1;
  // 가장 가까운 y-축선의 x-좌표를 양의 무한대로 초기화
  double axialCoord1 = numeric_limits<double>::max ();
  // 두번째로 가까운 y-축선의 x-좌표를 양의 무한대로 초기화
  double axialCoord2 = numeric_limits<double>::max ();
  // y-축선의 index를 0부터 y-축선의 개수만큼 반복
  for (size_t i = 0; i < adat->XXYYaxial_Num ('y'); i++) {
    // 반복하는 y-축선의 index가 현재점을 포함하는 y-축선의 index와 같으면 무시
    if (this->Mark ().compare ("p") == 0) if (this->Mark ().compare ("p") == 0) if (i == adat->Axial_Index (this->Index (), 'y')) continue;
    // y-축선의 양끝점의 y-좌표 사이에 현재점의 y-좌표가 있는 경우
    if ((adat->XYaxial ('y', i, 1) - this->Coord ().Value ('y')) * (adat->XYaxial ('y', i, 2) - this->Coord ().Value ('y')) < NearZero) {
      // y-축선의 x-좌표가 현재점의 x-좌표보다 큰 경우
      if (adat->XYaxial ('y', i, 0) > this->Coord ().Value ('x') + NearZero) {
        // y-축선의 x-좌표가 이전에 찾은 가장 가까운 y-축선의 x-좌표보다 더 작은 경우 (현재점과 더 가까운 경우))
        if (adat->XYaxial ('y', i, 0) < axialCoord1) {
          // 가장 가까운 y-축선의 index를 두번째로 가까운 y-축선의 index에 저장
          axialIndex2 = axialIndex1;
          // 가장 가까운 y-축선의 x-좌표를 두번째로 가까운 y-축선의 x-좌표에 저장
          axialCoord2 = axialCoord1;
          // 찾은 y-축선의 index를 가장 가까운 y-축선의 index에 저장
          axialIndex1 = i;
          // 찾은 y-축선의 x-좌표를 가장 가까운 y-축선의 x-좌표에 저장
          axialCoord1 = adat->XYaxial ('y', axialIndex1, 0);
        }
        // y-축선의 x-좌표가 이전에 찾은 가장 가까운 y-축선과 두번째로 가까운 y-축선의 x-좌표 사이에 있는 경우
        else if (adat->XYaxial ('y', i, 0) < axialCoord2) {
          // 찾은 y-축선의 index를 두번째로 가까운 y-축선의 index에 저장
          axialIndex2 = i;
          // 찾은 y-축선의 x-좌표를 두번째로 가까운 y-축선의 x-좌표에 저장
          axialCoord2 = adat->XYaxial ('y', axialIndex2, 0);
        }
      }
    }
  }
  // 가장 가까운 y-축선을 찾지 못한 경우
  if (axialIndex1 > adat->XXYYaxial_Num ('y')) return -1;
  // 두번째로 가까운 y-축선을 찾지 못한 경우
  if (axialIndex2 > adat->XXYYaxial_Num ('y')) return -1;
  // 두번째로 가까운 y-축선을 return
  return axialIndex2;
}

/* 두번재로 가까운 왼쪽의 y-축선의 index를 찾는 모듈 */
int Point::FindWest2ndAxialLine (AxialData *adat) {
  // 가장 가까운 y-축선의 index를 y-축선의 개수 + 1로 초기화
  int axialIndex1 = adat->XXYYaxial_Num ('y') + 1;
  // 두번째로 가까운 y-축선의 index를 y-축선의 개수 + 1로 초기화
  int axialIndex2 = adat->XXYYaxial_Num ('y') + 1;
  // 가장 가까운 y-축선의 x-좌표를 음의 무한대로 초기화
  double axialCoord1 = -numeric_limits<double>::max ();
  // 두번째로 가까운 y-축선의 x-좌표를 음의 무한대로 초기화
  double axialCoord2 = -numeric_limits<double>::max ();
  // y-축선의 index를 0부터 y-축선의 개수만큼 반복
  for (size_t i = 0; i < adat->XXYYaxial_Num ('y'); i++) {
    // 반복하는 y-축선의 index가 현재점을 포함하는 y-축선의 index와 같으면 무시
    if (this->Mark ().compare ("p") == 0) if (i == adat->Axial_Index (this->Index (), 'y')) continue;
    // y-축선의 양끝점의 y-좌표 사이에 현재점의 y-좌표가 있는 경우
    if ((adat->XYaxial ('y', i, 1) - this->Coord ().Value ('y')) * (adat->XYaxial ('y', i, 2) - this->Coord ().Value ('y')) < NearZero) {
      // y-축선의 x-좌표가 현재점의 x-좌표보다 작은 경우
      if (adat->XYaxial ('y', i, 0) < this->Coord ().Value ('x') - NearZero) {
        // y-축선의 x-좌표가 이전에 찾은 가장 가까운 y-축선의 x-좌표보다 더 큰 경우 (현재점과 더 가까운 경우))
        if (adat->XYaxial ('y', i, 0) > axialCoord1) {
          // 가장 가까운 y-축선의 index를 두번째로 가까운 y-축선의 index에 저장
          axialIndex2 = axialIndex1;
          // 가장 가까운 y-축선의 x-좌표를 두번째로 가까운 y-축선의 x-좌표에 저장
          axialCoord2 = axialCoord1;
          // 찾은 y-축선의 index를 가장 가까운 y-축선의 index에 저장
          axialIndex1 = i;
          // 찾은 y-축선의 x-좌표를 가장 가까운 y-축선의 x-좌표에 저장
          axialCoord1 = adat->XYaxial ('y', axialIndex1, 0);
        }
        // y-축선의 x-좌표가 이전에 찾은 가장 가까운 y-축선과 두번째로 가까운 y-축선의 x-좌표 사이에 있는 경우
        else if (adat->XYaxial ('y', i, 0) > axialCoord2) {
          // 찾은 y-축선의 index를 두번째로 가까운 y-축선의 index에 저장
          axialIndex2 = i;
          // 찾은 y-축선의 x-좌표를 두번째로 가까운 y-축선의 x-좌표에 저장
          axialCoord2 = adat->XYaxial ('y', axialIndex2, 0);
        }
      }
    }
  }
  // 가장 가까운 y-축선을 찾지 못한 경우
  if (axialIndex1 > adat->XXYYaxial_Num ('y')) return -1;
  // 두번째로 가까운 y-축선을 찾지 못한 경우
  if (axialIndex2 > adat->XXYYaxial_Num ('y')) return -1;
  // 두번째로 가까운 y-축선을 return
  return axialIndex2;
}

/* 두번재로 가까운 위쪽의 x-축선의 index를 찾는 모듈 */
int Point::FindNorth2ndAxialLine (AxialData *adat) {
  // 가장 가까운 x-축선의 index를 x-축선의 개수 + 1로 초기화
  int axialIndex1 = adat->XXYYaxial_Num ('x') + 1;
  // 두번째로 가까운 x-축선의 index를 x-축선의 개수 + 1로 초기화
  int axialIndex2 = adat->XXYYaxial_Num ('x') + 1;
  // 가장 가까운 x-축선의 y-좌표를 양의 무한대로 초기화
  double axialCoord1 = numeric_limits<double>::max ();
  // 두번째로 가까운 x-축선의 y-좌표를 양의 무한대로 초기화
  double axialCoord2 = numeric_limits<double>::max ();
  // x-축선의 index를 0부터 x-축선의 개수만큼 반복
  for (size_t i = 0; i < adat->XXYYaxial_Num ('x'); i++) {
    // 반복하는 x-축선의 index가 현재점을 포함하는 x-축선의 index와 같으면 무시
    if (this->Mark ().compare ("p") == 0) if (i == adat->Axial_Index (this->Index (), 'x')) continue;
    // x-축선의 양끝점의 x-좌표 사이에 현재점의 x-좌표가 있는 경우
    if ((adat->XYaxial ('x', i, 1) - this->Coord ().Value ('x')) * (adat->XYaxial ('x', i, 2) - this->Coord ().Value ('x')) < NearZero) {
      // x-축선의 y-좌표가 현재점의 y-좌표보다 큰 경우
      if (adat->XYaxial ('x', i, 0) > this->Coord ().Value ('y') + NearZero) {
        // x-축선의 y-좌표가 이전에 찾은 가장 가까운 x-축선의 y-좌표보다 더 작은 경우 (현재점과 더 가까운 경우))
        if (adat->XYaxial ('x', i, 0) < axialCoord1) {
          // 가장 가까운 x-축선의 index를 두번째로 가까운 x-축선의 index에 저장
          axialIndex2 = axialIndex1;
          // 가장 가까운 x-축선의 y-좌표를 두번째로 가까운 x-축선의 y-좌표에 저장
          axialCoord2 = axialCoord1;
          // 찾은 x-축선의 index를 가장 가까운 x-축선의 index에 저장
          axialIndex1 = i;
          // 찾은 x-축선의 y-좌표를 가장 가까운 x-축선의 y-좌표에 저장
          axialCoord1 = adat->XYaxial ('x', axialIndex1, 0);
        }
        // x-축선의 y-좌표가 이전에 찾은 가장 가까운 x-축선과 두번째로 가까운 x-축선의 y-좌표 사이에 있는 경우
        else if (adat->XYaxial ('x', i, 0) < axialCoord2) {
          // 찾은 x-축선의 index를 두번째로 가까운 x-축선의 index에 저장
          axialIndex2 = i;
          // 찾은 x-축선의 y-좌표를 두번째로 가까운 x-축선의 y-좌표에 저장
          axialCoord2 = adat->XYaxial ('x', axialIndex2, 0);
        }
      }
    }
  }
  // 가장 가까운 x-축선을 찾지 못한 경우
  if (axialIndex1 > adat->XXYYaxial_Num ('x')) return -1;
  // 두번째로 가까운 x-축선을 찾지 못한 경우
  if (axialIndex2 > adat->XXYYaxial_Num ('x')) return -1;
  // 두번째로 가까운 x-축선을 return
  return axialIndex2;
}

/* 두번재로 가까운 아래쪽의 x-축선의 index를 찾는 모듈 */
int Point::FindSouth2ndAxialLine (AxialData *adat) {
  // 가장 가까운 x-축선의 index를 x-축선의 개수 + 1로 초기화
  int axialIndex1 = adat->XXYYaxial_Num ('x') + 1;
  // 두번째로 가까운 x-축선의 index를 x-축선의 개수 + 1로 초기화
  int axialIndex2 = adat->XXYYaxial_Num ('x') + 1;
  // 가장 가까운 x-축선의 y-좌표를 음의 무한대로 초기화
  double axialCoord1 = -numeric_limits<double>::max ();
  // 두번째로 가까운 x-축선의 y-좌표를 음의 무한대로 초기화
  double axialCoord2 = -numeric_limits<double>::max ();
  // x-축선의 index를 0부터 x-축선의 개수만큼 반복
  for (size_t i = 0; i < adat->XXYYaxial_Num ('x'); i++) {
    // 반복하는 x-축선의 index가 현재점을 포함하는 x-축선의 index와 같으면 무시
    if (this->Mark ().compare ("p") == 0) if (i == adat->Axial_Index (this->Index (), 'x')) continue;
    // x-축선의 양끝점의 x-좌표 사이에 현재점의 x-좌표가 있는 경우
    if ((adat->XYaxial ('x', i, 1) - this->Coord ().Value ('x')) * (adat->XYaxial ('x', i, 2) - this->Coord ().Value ('x')) < NearZero) {
      // x-축선의 y-좌표가 현재점의 y-좌표보다 작은 경우
      if (adat->XYaxial ('x', i, 0) < this->Coord ().Value ('y') - NearZero) {
        // x-축선의 y-좌표가 이전에 찾은 가장 가까운 x-축선의 y-좌표보다 더 큰 경우 (현재점과 더 가까운 경우))
        if (adat->XYaxial ('x', i, 0) > axialCoord1) {
          // 가장 가까운 x-축선의 index를 두번째로 가까운 x-축선의 index에 저장
          axialIndex2 = axialIndex1;
          // 가장 가까운 x-축선의 y-좌표를 두번째로 가까운 x-축선의 y-좌표에 저장
          axialCoord2 = axialCoord1;
          // 찾은 x-축선의 index를 가장 가까운 x-축선의 index에 저장
          axialIndex1 = i;
          // 찾은 x-축선의 y-좌표를 가장 가까운 x-축선의 y-좌표에 저장
          axialCoord1 = adat->XYaxial ('x', axialIndex1, 0);
        }
        // x-축선의 y-좌표가 이전에 찾은 가장 가까운 x-축선과 두번째로 가까운 x-축선의 y-좌표 사이에 있는 경우
        else if (adat->XYaxial ('x', i, 0) > axialCoord2) {
          // 찾은 x-축선의 index를 두번째로 가까운 x-축선의 index에 저장
          axialIndex2 = i;
          // 찾은 x-축선의 y-좌표를 두번째로 가까운 x-축선의 y-좌표에 저장
          axialCoord2 = adat->XYaxial ('x', axialIndex2, 0);
        }
      }
    }
  }
  // 가장 가까운 x-축선을 찾지 못한 경우
  if (axialIndex1 > adat->XXYYaxial_Num ('x')) return -1;
  // 두번째로 가까운 x-축선을 찾지 못한 경우
  if (axialIndex2 > adat->XXYYaxial_Num ('x')) return -1;
  // 두번째로 가까운 x-축선을 return
  return axialIndex2;
}

/* y-축선위에 있는 점들중에서 현재점을 사이에 두는 가장 가까운 점을 찾는 모듈 */
int Point::FindVerticalPoints (AxialData *adat, int axialIndex, char EWNS) {
  // 참조한 y-축선위 시작점
  int startPoint = adat->XXYYaxial_Index ('y', axialIndex);
  // y-축선의 다음 y-축선의 시작점으로
  int endPoint   = adat->XXYYaxial_Index ('y', axialIndex + 1);
  // 찾을 점의 index를 -1로 초기화
  int index = -1;
  // 검색하는 점의 index
  int i;
  // 현재점과 가장 가까운 아래쪽점의 y-좌표를 음의 무한대로 초기화
  double ym = -numeric_limits<double>::max ();
  // 현재점의 y-좌표
  double yb = this->Coord ().Value ('y');
  // 현재점과 가장 가까운 위쪽점의 y-좌표를 양의 무한대로 초기화
  double yp = numeric_limits<double>::max ();
  // 검색하는 점의 y-좌표를 저장하기위한 임시변수
  double tmp;
  // 참조하는 y-축선의 index가 -1인 경우, 현재점의 index를 return
  if (axialIndex == -1) return this->Index ();
  // 현재점과 가장 가까운 위쪽점을 찾는 경우
  if (EWNS == 'N') {
    // 참조한 y-축선의 시작점에서 끝점까지 반복
    for (size_t itr = startPoint; itr < endPoint; itr++) {
      // 검색하는 점의 index를 저장
      i = adat->XYaxial_Index ('y', itr);
      // 검색하는 점의 y-좌표를 저장
      tmp = adat->Pts (i, 'y');
      // 검색하는 점의 y-좌표가 현재점의 y-좌표보다 더 큰 경우
      if ((tmp + NearZero) > yb) {
        // 검색하는 점의 y-좌표가 찾은 가장 가까운 위쪽점의 y-좌표보다 더 작은 경우
        if (tmp < yp) {
          // 검색한 점의 y-좌표를 저장
          yp = tmp;
          // 검색한 점의 index를 저장
          index = i;
        }
      }
    }
    // 찾은 점의 index를 return
    return index;
  }
  // 현재점과 가장 가까운 아래쪽점을 찾는 경우
  if (EWNS == 'S') {
    // 참조한 y-축선의 시작점에서 끝점까지 반복
    for (size_t itr = startPoint; itr < endPoint; itr++) {
      // 검색하는 점의 index를 저장
      i = adat->XYaxial_Index ('y', itr);
      // 검색하는 점의 y-좌표를 저장
      tmp = adat->Pts (i, 'y');
      // 검색하는 점의 y-좌표가 현재점의 y-좌표보다 더 작은 경우
      if ((tmp - NearZero) < yb) {
        // 검색하는 점의 y-좌표가 찾은 가장 가까운 아래쪽점의 y-좌표보다 더 큰 경우
        if (tmp > ym) {
          // 검색한 점의 y-좌표를 저장
          ym = tmp;
          // 검색한 점의 index를 저장
          index = i;
        }
      }
    }
    // 찾은 점의 index를 return
    return index;
  }
  // 찾못된 참조의 경우 에러메시지를 출력하고 종료
  PrintError ("Point::FindVerticalPoints");
  exit (1);
}

/* x-축선위에 있는 점들중에서 현재점을 사이에 두는 가장 가까운 점을 찾는 모듈 */
int Point::FindHorizontalPoints (AxialData *adat, int axialIndex, char EWNS) {
  // 참조한 x-축선위 시작점
  int startPoint = adat->XXYYaxial_Index ('x', axialIndex);
  // x-축선의 다음 x-축선의 시작점으로
  int endPoint   = adat->XXYYaxial_Index ('x', axialIndex + 1);
  // 찾을 점의 index를 -1로 초기화
  int index = -1;
  // 검색하는 점의 index
  int i;
  // 현재점과 가장 가까운 왼쪽점의 x-좌표를 음의 무한대로 초기화
  double xm = -numeric_limits<double>::max ();
  // 현재점의 x-좌표
  double xb = this->Coord ().Value ('x');
  // 현재점과 가장 가까운 오른쪽점의 x-좌표를 양의 무한대로 초기화
  double xp = numeric_limits<double>::max ();
  // 검색하는 점의 x-좌표를 저장하기위한 임시변수
  double tmp;
  // 참조하는 x-축선의 index가 -1인 경우, 현재점의 index를 return
  if (axialIndex == -1) return this->Index ();
  // 현재점과 가장 가까운 오른쪽점을 찾는 경우
  if (EWNS == 'E') {
    // 참조한 x-축선의 시작점에서 끝점까지 반복
    for (int itr = startPoint; itr < endPoint; itr++) {
      // 검색하는 점의 index를 저장
      i = adat->XYaxial_Index ('x', itr);
      // 검색하는 점의 x-좌표를 저장
      tmp = adat->Pts (i, 'x');
      // 검색하는 점의 x-좌표가 현재점의 x-좌표보다 더 큰 경우
      if ((tmp + NearZero) > xb) {
        // 검색하는 점의 x-좌표가 찾은 가장 가까운 오른쪽점의 x-좌표보다 더 작은 경우
        if (tmp < xp) {
          // 검색한 점의 x-좌표를 저장
          xp = tmp;
          // 검색한 점의 index를 저장
          index = i;
        }
      }
    }
    // 찾은 점의 index를 return
    return index;
  }
  // 현재점과 가장 가까운 왼쪽점을 찾는 경우
  if (EWNS == 'W') {
    // 참조한 x-축선의 시작점에서 끝점까지 반복
    for (int itr = startPoint; itr < endPoint; itr++) {
      // 검색하는 점의 index를 저장
      i = adat->XYaxial_Index ('x', itr);
      // 검색하는 점의 x-좌표를 저장
      tmp = adat->Pts (i, 'x');
      // 검색하는 점의 x-좌표가 현재점의 x-좌표보다 더 작은 경우
      if ((tmp - NearZero) < xb) {
        // 검색하는 점의 x-좌표가 찾은 가장 가까운 왼쪽점의 x-좌표보다 더 큰 경우
        if (tmp > xm) {
          // 검색한 점의 x-좌표를 저장
          xm = tmp;
          // 검색한 점의 index를 저장
          index = i;
        }
      }
    }
    // 찾은 점의 index를 return
    return index;
  }
  // 찾못된 참조의 경우 에러메시지를 출력하고 종료
  PrintError ("Point::FindHorizontalPoints");
  exit (1);
}

/* 현재점의 오른쪽, 왼쪽, 위쪽, 아래쪽의 점의 주소를 받아서 저장하는 모듈 */
Point & Point::SetEWNS (char EWNS, char ewns, Point *src) {
  // 현재점의 오른쪽에 있는 점을 저장하는 경우
  if (EWNS == 'E' || EWNS == 'e') {
    // 현재점의 오른쪽점을 저장하는 경우
    if (ewns == 'E' || ewns == 'e') {this->E = src; return *this;}
    // 잘못된 위치 참조, 에러메시지를 출력하고 종료
    if (ewns == 'W' || ewns == 'w') {PrintError ("Point::SetEWNS, EW");}
    // 현재점의 오른쪽의 위쪽점을 저장하는 경우
    if (ewns == 'N' || ewns == 'n') {this->EN = src; return *this;}
    // 현재점의 오른쪽의 아래쪽점을 저장하는 경우
    if (ewns == 'S' || ewns == 's') {this->ES = src; return *this;}
  }
  // 현재점의 왼쪽에 있는 점을 저장하는 경우
  if (EWNS == 'W' || EWNS == 'w') {
    // 잘못된 위치 참조, 에러메시지를 출력하고 종료
    if (ewns == 'E' || ewns == 'e') {PrintError ("Point::SetEWNS, WE");}
    // 현재점의 왼쪽점을 저장하는 경우
    if (ewns == 'W' || ewns == 'w') {this->W = src; return *this;}
    // 현재점의 왼쪽의 위쪽점을 저장하는 경우
    if (ewns == 'N' || ewns == 'n') {this->WN = src; return *this;}
    // 현재점의 왼쪽의 아래쪽점을 저장하는 경우
    if (ewns == 'S' || ewns == 's') {this->WS = src; return *this;}
  }
  // 현재점의 위쪽에 있는 점을 저장하는 경우
  if (EWNS == 'N' || EWNS == 'n') {
    // 현재점의 위쪽의 오른쪽점을 저장하는 경우
    if (ewns == 'E' || ewns == 'e') {this->NE = src; return *this;}
    // 현재점의 위쪽의 왼쪽점을 저장하는 경우
    if (ewns == 'W' || ewns == 'w') {this->NW = src; return *this;}
    // 현재점의 위쪽점을 저장하는 경우
    if (ewns == 'N' || ewns == 'n') {this->N = src; return *this;}
    // 잘못된 위치 참조, 에러메시지를 출력하고 종료
    if (ewns == 'S' || ewns == 's') {PrintError ("Point::SetEWNS, NS");}
  }
  // 현재점의 아래쪽에 있는 점을 저장하는 경우
  if (EWNS == 'S' || EWNS == 's') {
    // 현재점의 아래쪽의 오른쪽점을 저장하는 경우
    if (ewns == 'E' || ewns == 'e') {this->SE = src; return *this;}
    // 현재점의 아래쪽의 왼쪽점을 저장하는 경우
    if (ewns == 'W' || ewns == 'w') {this->SW = src; return *this;}
    // 잘못된 위치 참조, 에러메시지를 출력하고 종료
    if (ewns == 'N' || ewns == 'n') {PrintError ("Point::SetEWNS, SN");}
    // 현재점의 아래쪽점을 저장하는 경우
    if (ewns == 'S' || ewns == 's') {this->S = src; return *this;}
  }
  // 잘못된 위치 참조, 에러메시지를 출력하고 종료
  PrintError ("Point::SetEWNS");
  exit (1);
}

/* 현재점의 두번째 오른쪽, 왼쪽, 위쪽, 아래쪽의 점의 주소를 받아서 저장하는 모듈 */
Point & Point::SetEWNS2nd (char EWNS, char ewns, Point *src) {
  // 현재점의 두번째 오른쪽에 있는 점을 저장하는 경우
  if (EWNS == 'E' || EWNS == 'e') {
    // 현재점의 두번째 오른쪽점을 저장하는 경우
    if (ewns == 'E' || ewns == 'e') {this->e = src; return *this;}
    // 잘못된 위치 참조, 에러메시지를 출력하고 종료
    if (ewns == 'W' || ewns == 'w') {PrintError ("Point::SetEWNS2nd, EW");}
    // 현재점의 두번째 오른쪽의 위쪽점을 저장하는 경우
    if (ewns == 'N' || ewns == 'n') {this->en = src; return *this;}
    // 현재점의 두번째 오른쪽의 아래쪽점을 저장하는 경우
    if (ewns == 'S' || ewns == 's') {this->es = src; return *this;}
  }
  // 현재점의 두번째 왼쪽에 있는 점을 저장하는 경우
  if (EWNS == 'W' || EWNS == 'w') {
    // 잘못된 위치 참조, 에러메시지를 출력하고 종료
    if (ewns == 'E' || ewns == 'e') {PrintError ("Point::SetEWNS2nd, WE");}
    // 현재점의 두번째 왼쪽점을 저장하는 경우
    if (ewns == 'W' || ewns == 'w') {this->w = src; return *this;}
    // 현재점의 두번째 왼쪽의 위쪽점을 저장하는 경우
    if (ewns == 'N' || ewns == 'n') {this->wn = src; return *this;}
    // 현재점의 두번째 왼쪽의 아래쪽점을 저장하는 경우
    if (ewns == 'S' || ewns == 's') {this->ws = src; return *this;}
  }
  // 현재점의 두번째 위쪽에 있는 점을 저장하는 경우
  if (EWNS == 'N' || EWNS == 'n') {
    // 현재점의 두번째 위쪽의 오른쪽점을 저장하는 경우
    if (ewns == 'E' || ewns == 'e') {this->ne = src; return *this;}
    // 현재점의 두번째 위쪽의 왼쪽점을 저장하는 경우
    if (ewns == 'W' || ewns == 'w') {this->nw = src; return *this;}
    // 현재점의 두번째 위쪽점을 저장하는 경우
    if (ewns == 'N' || ewns == 'n') {this->n = src; return *this;}
    // 잘못된 위치 참조, 에러메시지를 출력하고 종료
    if (ewns == 'S' || ewns == 's') {PrintError ("Point::SetEWNS2nd, NS");}
  }
  // 현재점의 두번째 아래쪽에 있는 점을 저장하는 경우
  if (EWNS == 'S' || EWNS == 's') {
    // 현재점의 두번째 아래쪽의 오른쪽점을 저장하는 경우
    if (ewns == 'E' || ewns == 'e') {this->se = src; return *this;}
    // 현재점의 두번째 아래쪽의 왼쪽점을 저장하는 경우
    if (ewns == 'W' || ewns == 'w') {this->sw = src; return *this;}
    // 잘못된 위치 참조, 에러메시지를 출력하고 종료
    if (ewns == 'N' || ewns == 'n') {PrintError ("Point::SetEWNS2nd, SN");}
    // 현재점의 두번째 아래쪽점을 저장하는 경우
    if (ewns == 'S' || ewns == 's') {this->s = src; return *this;}
  }
  // 잘못된 위치 참조, 에러메시지를 출력하고 종료
  PrintError ("Point::SetEWNS2nd");
  exit (1);
}

/* 현재점의 오른쪽, 왼쪽, 위쪽, 아래쪽의 점의 index를 내보내는 모듈 */
Point & Point::ExportEWNS (FILE* fp) {
  // 현재점의 index
  fprintf (fp,"%d\t", this->Index ());
  // 현재점의 오른쪽점이 존재하는 경우, 오른쪽점의 index를 쓰기
  if (this->EWNS ('E', 'E') != NULL) fprintf (fp, "%9d\t", this->EWNS ('E', 'E')->Index ());
  // 현재점의 오른쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 왼쪽점이 존재하는 경우, 왼쪽점의 index를 쓰기
  if (this->EWNS ('W', 'W') != NULL) fprintf (fp, "%9d\t", this->EWNS ('W', 'W')->Index ());
  // 현재점의 왼쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 위쪽점이 존재하는 경우, 위쪽점의 index를 쓰기
  if (this->EWNS ('N', 'N') != NULL) fprintf (fp, "%9d\t", this->EWNS ('N', 'N')->Index ());
  // 현재점의 위쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 아래쪽점이 존재하는 경우, 아래쪽점의 index를 쓰기
  if (this->EWNS ('S', 'S') != NULL) fprintf (fp, "%9d\t", this->EWNS ('S', 'S')->Index ());
  // 현재점의 아래쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 오른쪽의 위쪽점이 존재하는 경우, 오른쪽의 위쪽점의 index를 쓰기
  if (this->EWNS ('E', 'N') != NULL) fprintf (fp,"%9d\t", this->EWNS ('E', 'N')->Index ());
  // 현재점의 오른쪽의 위쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 오른쪽의 아래쪽점이 존재하는 경우, 오른쪽의 아래쪽점의 index를 쓰기
  if (this->EWNS ('E', 'S') != NULL) fprintf (fp,"%9d\t", this->EWNS ('E', 'S')->Index ());
  // 현재점의 오른쪽의 아래쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 왼쪽의 위쪽점이 존재하는 경우, 왼쪽의 위쪽점의 index를 쓰기
  if (this->EWNS ('W', 'N') != NULL) fprintf (fp,"%9d\t", this->EWNS ('W', 'N')->Index ());
  // 현재점의 왼쪽의 위쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 왼쪽의 아래쪽점이 존재하는 경우, 왼쪽의 아래쪽점의 index를 쓰기
  if (this->EWNS ('W', 'S') != NULL) fprintf (fp,"%9d\t", this->EWNS ('W', 'S')->Index ());
  // 현재점의 왼쪽의 아래쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 위쪽의 오른쪽점이 존재하는 경우, 위쪽의 오른쪽점의 index를 쓰기
  if (this->EWNS ('N', 'E') != NULL) fprintf (fp,"%9d\t", this->EWNS ('N', 'E')->Index ());
  // 현재점의 위쪽의 오른쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 위쪽의 왼쪽점이 존재하는 경우, 위쪽의 왼쪽점의 index를 쓰기
  if (this->EWNS ('N', 'W') != NULL) fprintf (fp,"%9d\t", this->EWNS ('N', 'W')->Index ());
  // 현재점의 위쪽의 왼쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 아래쪽의 오른쪽점이 존재하는 경우, 아래쪽의 오른쪽점의 index를 쓰기
  if (this->EWNS ('S', 'E') != NULL) fprintf (fp,"%9d\t", this->EWNS ('S', 'E')->Index ());
  // 현재점의 아래쪽의 오른쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 아래쪽의 왼쪽점이 존재하는 경우, 아래쪽의 왼쪽점의 index를 쓰기
  if (this->EWNS ('S', 'W') != NULL) fprintf (fp,"%9d\n", this->EWNS ('S', 'W')->Index ());
  // 현재점의 아래쪽의 왼쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\n", -1);
  return *this;
}

/* 현재점의 두번째 오른쪽, 왼쪽, 위쪽, 아래쪽의 점의 index를 내보내는 모듈 */
Point & Point::ExportEWNS2nd (FILE* fp) {
  // 현재점의 index
  fprintf (fp,"%d\t", this->Index ());
  // 현재점의 오른쪽점의 오른쪽점이 존재하는 경우, 오른쪽점의 index를 쓰기
  if (this->EWNS2nd ('E', 'E') != NULL) fprintf (fp, "%9d\t", this->EWNS2nd ('E', 'E')->Index ());
  // 현재점의 오른쪽점의 오른쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 왼쪽점의 왼쪽점이 존재하는 경우, 왼쪽점의 index를 쓰기
  if (this->EWNS2nd ('W', 'W') != NULL) fprintf (fp, "%9d\t", this->EWNS2nd ('W', 'W')->Index ());
  // 현재점의 왼쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 위쪽점의 왼쪽점의 위쪽점이 존재하는 경우, 위쪽점의 index를 쓰기
  if (this->EWNS2nd ('N', 'N') != NULL) fprintf (fp, "%9d\t", this->EWNS2nd ('N', 'N')->Index ());
  // 현재점의 위쪽점의 위쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 아래쪽점의 아래쪽점이 존재하는 경우, 아래쪽점의 index를 쓰기
  if (this->EWNS2nd ('S', 'S') != NULL) fprintf (fp, "%9d\t", this->EWNS2nd ('S', 'S')->Index ());
  // 현재점의 아래쪽점의 아래쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 오른쪽점의 오른쪽의 위쪽점이 존재하는 경우, 오른쪽의 위쪽점의 index를 쓰기
  if (this->EWNS2nd ('E', 'N') != NULL) fprintf (fp,"%9d\t", this->EWNS2nd ('E', 'N')->Index ());
  // 현재점의 오른쪽점의 오른쪽의 위쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 오른쪽점의 오른쪽의 아래쪽점이 존재하는 경우, 오른쪽의 아래쪽점의 index를 쓰기
  if (this->EWNS2nd ('E', 'S') != NULL) fprintf (fp,"%9d\t", this->EWNS2nd ('E', 'S')->Index ());
  // 현재점의 오른쪽점의 오른쪽의 아래쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 왼쪽점의 왼쪽의 위쪽점이 존재하는 경우, 왼쪽의 위쪽점의 index를 쓰기
  if (this->EWNS2nd ('W', 'N') != NULL) fprintf (fp,"%9d\t", this->EWNS2nd ('W', 'N')->Index ());
  // 현재점의 왼쪽점의 왼쪽의 위쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 왼쪽점의 왼쪽의 아래쪽점이 존재하는 경우, 왼쪽의 아래쪽점의 index를 쓰기
  if (this->EWNS2nd ('W', 'S') != NULL) fprintf (fp,"%9d\t", this->EWNS2nd ('W', 'S')->Index ());
  // 현재점의 왼쪽점의 왼쪽의 아래쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 위쪽점의 위쪽의 오른쪽점이 존재하는 경우, 위쪽의 오른쪽점의 index를 쓰기
  if (this->EWNS2nd ('N', 'E') != NULL) fprintf (fp,"%9d\t", this->EWNS2nd ('N', 'E')->Index ());
  // 현재점의 위쪽점의 위쪽의 오른쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 위쪽점의 위쪽의 왼쪽점이 존재하는 경우, 위쪽의 왼쪽점의 index를 쓰기
  if (this->EWNS2nd ('N', 'W') != NULL) fprintf (fp,"%9d\t", this->EWNS2nd ('N', 'W')->Index ());
  // 현재점의 위쪽점의 위쪽의 왼쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 아래쪽점의 아래쪽의 오른쪽점이 존재하는 경우, 아래쪽의 오른쪽점의 index를 쓰기
  if (this->EWNS2nd ('S', 'E') != NULL) fprintf (fp,"%9d\t", this->EWNS2nd ('S', 'E')->Index ());
  // 현재점의 아래쪽점의 아래쪽의 오른쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\t", -1);
  // 현재점의 아래쪽점의 아래쪽의 왼쪽점이 존재하는 경우, 아래쪽의 왼쪽점의 index를 쓰기
  if (this->EWNS2nd ('S', 'W') != NULL) fprintf (fp,"%9d\n", this->EWNS2nd ('S', 'W')->Index ());
  // 현재점의 아래쪽점의 아래쪽의 왼쪽점이 존재하지 않는 경우, -1을 쓰기
  else                               fprintf (fp, "%9d\n", -1);
  return *this;
}

bool Point::IsBoundary () {
  char azimuth[4] = {'E', 'W', 'N' ,'S'};

  for (const auto &i : azimuth) if (this->IsBoundary (i)) return true;
  return false;
}

bool Point::IsBoundary (const char azimuth) {

  if (this->Condition () == 'C' || this->Condition () == 'I' || this->Condition () == 'M') return false;

  char axis = ' ';
  char errorMassage[1024];
  double sign = 1.0E0;

  if      (azimuth == 'E' || azimuth == 'W') axis = 'x';
  else if (azimuth == 'N' || azimuth == 'S') axis = 'y';
  else sprintf (errorMassage, "Error in Point::IsBoundary, azimuth = %c, please check azimuth\n", azimuth), PrintError (errorMassage), exit (1);

  if (IsEqualDouble (this->Normal ().Value (axis), ZeroValue)) return false;
  if (azimuth == 'W' || azimuth == 'S') sign = -1.0E0;

  if      (this->Pressure ()->Normal ().Value (axis) * sign > 0) return true;
  else if (this->Pressure ()->Normal ().Value (axis) * sign < 0) return false;

  sprintf (errorMassage, "Error in Point::IsBoundary, Something is wrong\n"), PrintError (errorMassage);
  exit (1);
}

double Point::PrintDebuggingData (const string verbose, xData *xdat, yData *ydat, bool is_sol) {

  int verboseCoord = int (verbose[0]) - '0';
  int verboseBound = int (verbose[1]) - '0';
  int verboseValue = int (verbose[2]) - '0';
  int verboseazimu = int (verbose[3]) - '0';
  int verboseRepre = int (verbose[4]) - '0';

  char azimuth[4] = {'E', 'W', 'N', 'S'};
  char sub_azimuth[] = {'E', 'E', 'E'};

  char massage[128];

  if (verboseCoord) printf ("\nPrintDebuggingData (verbose = %s): Point[%d] (%s)\n", verbose.c_str (), this->Index (), this->Mark ().c_str ());

  if (verboseCoord) HeadVelocity (this)->Coord ().ShowContents ("Coodinate = ");

  if (verboseBound) {
    if (this->Condition () != 'N') printf ("Boundary condition = %c, boundary value = %23.16e\n", HeadVelocity (this)->Condition (), HeadVelocity (this)->Boundaryvalue ());
    else                           printf ("Boundary condition = %c, boundary value = %23.16e\n", HeadVelocity (this)->Condition (), HeadVelocity (this)->Boundaryvalue ()), HeadVelocity (this)->Normal ().ShowContents ("Normal = ");
  }
  if (verboseValue) printf ("Value = %23.16e\n", this->Value ());

  if (verboseazimu) {
    printf ("\nisboundary = ");
    for (const auto &i : azimuth) if (HeadVelocity (this)->IsBoundary (i)) printf ("%c", i);
    printf ("\n");
  }

  if (verboseazimu) for (const auto &i : azimuth) {
    if      (i == 'E') sub_azimuth[0] = 'E', sub_azimuth[1] = 'N', sub_azimuth[2] = 'S';
    else if (i == 'W') sub_azimuth[0] = 'W', sub_azimuth[1] = 'N', sub_azimuth[2] = 'S';
    else if (i == 'N') sub_azimuth[0] = 'E', sub_azimuth[1] = 'W', sub_azimuth[2] = 'N';
    else if (i == 'S') sub_azimuth[0] = 'E', sub_azimuth[1] = 'W', sub_azimuth[2] = 'S';

    for (const auto &j : sub_azimuth) if (HeadVelocity (this)->EWNS (i, j)) sprintf (massage, "this->EWNS (%c, %c)", i, j), HeadVelocity (this)->EWNS (i, j)->PrintDebuggingData (verbose, massage);
  }

  double returnValue = ZeroValue;

  if (verboseRepre) {
    CalcRepresenCoef (this, xdat, ydat, is_sol);
    returnValue = HeadVelocity (this)->CheckRepresentationFormula (this, xdat, ydat, is_sol);
  }

  return returnValue;
}

Point & Point::PrintDebuggingData (const string verbose, const char *massage) {

  int verboseCoord = int (verbose[0]);
  int verboseBound = int (verbose[1]);
  int verboseValue = int (verbose[2]);

  printf ("\nPrintDebuggingData (verbose = %s): Point[%d] (%s)\n", verbose.c_str (), this->Index (), massage);

  if (verboseCoord) this->Coord ().ShowContents ("Coodinate = ");

  if (verboseBound) {
    if (this->Condition () != 'N') printf ("Boundary condition = %c, boundary value = %23.16e\n", this->Condition (), this->Boundaryvalue ());
    else                           printf ("Boundary condition = %c, boundary value = %23.16e\n", this->Condition (), this->Boundaryvalue ()), this->Normal ().ShowContents ("Normal = ");
  }
  if (verboseValue) printf ("Value = %23.16e\n", this->Value ());

  return *this;
}

double Point::CheckRepresentationFormula (Point *pt, xData *xdat, yData *ydat, bool is_sol) {

  double result[2] = {ZeroValue, ZeroValue};
  double ent_x[14] = {
    xdat->Cu, xdat->Eu, xdat->Wu, xdat->ENu, xdat->ESu, xdat->WNu, xdat->WSu,
    xdat->Cphi, xdat->Ephi, xdat->Wphi, xdat->ENphi, xdat->ESphi, xdat->WNphi, xdat->WSphi
  };

  double ent_y[14] = {
    ydat->Cu, ydat->Nu, ydat->Su, ydat->NEu, ydat->NWu, ydat->SEu, ydat->SWu,
    ydat->Cphi, ydat->Nphi, ydat->Sphi, ydat->NEphi, ydat->NWphi, ydat->SEphi, ydat->SWphi
  };

  Point *pt_x[14], *pt_y[14];

  pt_x[0] = this, pt_x[7] = this->Phi ();
  if (this->EWNS ('E', 'E')) pt_x[ 1] = this->EWNS ('E', 'E'), pt_x[ 8] = this->EWNS ('E', 'E')->Phi (); else pt_x[ 1] = NULL, pt_x[ 8] = NULL;
  if (this->EWNS ('W', 'W')) pt_x[ 2] = this->EWNS ('W', 'W'), pt_x[ 9] = this->EWNS ('W', 'W')->Phi (); else pt_x[ 2] = NULL, pt_x[ 9] = NULL;
  if (this->EWNS ('E', 'N')) pt_x[ 3] = this->EWNS ('E', 'N'), pt_x[10] = this->EWNS ('E', 'N')->Phi (); else pt_x[ 3] = NULL, pt_x[10] = NULL;
  if (this->EWNS ('E', 'S')) pt_x[ 4] = this->EWNS ('E', 'S'), pt_x[11] = this->EWNS ('E', 'S')->Phi (); else pt_x[ 4] = NULL, pt_x[11] = NULL;
  if (this->EWNS ('W', 'N')) pt_x[ 5] = this->EWNS ('W', 'N'), pt_x[12] = this->EWNS ('W', 'N')->Phi (); else pt_x[ 5] = NULL, pt_x[12] = NULL;
  if (this->EWNS ('W', 'S')) pt_x[ 6] = this->EWNS ('W', 'S'), pt_x[13] = this->EWNS ('W', 'S')->Phi (); else pt_x[ 6] = NULL, pt_x[13] = NULL;

  pt_y[0] = this, pt_y[7] = this->Phi ();
  if (this->EWNS ('N', 'N')) pt_y[ 1] = this->EWNS ('N', 'N'), pt_y[ 8] = this->EWNS ('N', 'N')->Phi (); else pt_y[ 1] = NULL, pt_y[ 8] = NULL;
  if (this->EWNS ('S', 'S')) pt_y[ 2] = this->EWNS ('S', 'S'), pt_y[ 9] = this->EWNS ('S', 'S')->Phi (); else pt_y[ 2] = NULL, pt_y[ 9] = NULL;
  if (this->EWNS ('N', 'E')) pt_y[ 3] = this->EWNS ('N', 'E'), pt_y[10] = this->EWNS ('N', 'E')->Phi (); else pt_y[ 3] = NULL, pt_y[10] = NULL;
  if (this->EWNS ('N', 'W')) pt_y[ 4] = this->EWNS ('N', 'W'), pt_y[11] = this->EWNS ('N', 'W')->Phi (); else pt_y[ 4] = NULL, pt_y[11] = NULL;
  if (this->EWNS ('S', 'E')) pt_y[ 5] = this->EWNS ('S', 'E'), pt_y[12] = this->EWNS ('S', 'E')->Phi (); else pt_y[ 5] = NULL, pt_y[12] = NULL;
  if (this->EWNS ('S', 'W')) pt_y[ 6] = this->EWNS ('S', 'W'), pt_y[13] = this->EWNS ('S', 'W')->Phi (); else pt_y[ 6] = NULL, pt_y[13] = NULL;

  if (*this == 'D' && !is_sol) {
    unordered_map <char, Point**> p, q, r, P, Q, R;
    unordered_map <char, char> o;
    p['E'] = &pt_x[0], p['W'] = &pt_x[0], p['N'] = &pt_y[0], p['S'] = &pt_y[0];
    q['E'] = &pt_x[1], q['W'] = &pt_x[2], q['N'] = &pt_y[1], q['S'] = &pt_y[2];
    r['E'] = &pt_x[2], r['W'] = &pt_x[1], r['N'] = &pt_y[2], r['S'] = &pt_y[1];
    P['E'] = &pt_x[7], P['W'] = &pt_x[7], P['N'] = &pt_y[7], P['S'] = &pt_y[7];
    Q['E'] = &pt_x[8], Q['W'] = &pt_x[9], Q['N'] = &pt_y[8], Q['S'] = &pt_y[9];
    R['E'] = &pt_x[9], R['W'] = &pt_x[8], R['N'] = &pt_y[9], R['S'] = &pt_y[8];
    o['E'] = 'W', o['W'] = 'E', o['N'] = 'S', o['S'] = 'N';

    unordered_map<char, double*> m;
    unordered_map<char, char> coord;
    unordered_map<char, string> axis;

    coord['E'] = 'x', coord['W'] = 'x', coord['N'] = 'y', coord['S'] = 'y';
    axis['x'] = "Y", axis['y'] = "X";

    for (const auto &i : {'W', 'E', 'S', 'N'})
    if (pt->Axis (coord[i]) > -1)
    if (pt->IsBoundary (o[i]))
    if (pt->EWNS (i, i))
    if (pt->EWNS (i, i)->EWNS (i, i))
    if (pt->EWNS (i, i)->EWNS (i, i)->Condition () == 'C') {
      *p[i] = pt->EWNS (i, i)->EWNS (i, i);
      *q[i] = (*p[i])->EWNS (i, i);
      *r[i] = this;
      *P[i] = (*p[i])->Phi ();
      *Q[i] = (*q[i])->Phi ();
      *R[i] = (*r[i])->Phi ();

      printf ("xxb, yyb = %23.16e, %23.16e\n", (*p[i])->Coord() [0], (*p[i])->Coord() [1]);
      printf ("xxp, yyb = %23.16e, %23.16e\n", (*p[i])->Coord() [0], (*p[i])->Coord() [1]);
      printf ("xxm, yyb = %23.16e, %23.16e\n", (*p[i])->Coord() [0], (*p[i])->Coord() [1]);
      printf ("xxb, yyp = %23.16e, %23.16e\n", (*r[i])->Coord() [0], (*r[i])->Coord() [1]);
      printf ("xxb, yym = %23.16e, %23.16e\n", (*q[i])->Coord() [0], (*q[i])->Coord() [1]);
      break;
    }
  }


  if (*this == 'N' && !is_sol) {
    unordered_map<char, unordered_map<char, Point**>> PT;
    unordered_map<char, unordered_map<char, Point**>> PHI;
    PT['E']['E'] = &pt_x[1], PT['E']['N'] = &pt_x[3], PT['E']['S'] = &pt_x[4];
    PT['W']['W'] = &pt_x[2], PT['W']['N'] = &pt_x[5], PT['W']['S'] = &pt_x[6];

    PT['N']['N'] = &pt_y[1], PT['N']['E'] = &pt_y[3], PT['N']['W'] = &pt_y[4];
    PT['S']['S'] = &pt_y[2], PT['S']['E'] = &pt_y[5], PT['S']['W'] = &pt_y[6];

    PHI['E']['E'] = &pt_x[8], PHI['E']['N'] = &pt_x[10], PHI['E']['S'] = &pt_x[11];
    PHI['W']['W'] = &pt_x[9], PHI['W']['N'] = &pt_x[12], PHI['W']['S'] = &pt_x[13];

    PHI['N']['N'] = &pt_y[8], PHI['N']['E'] = &pt_y[10], PHI['N']['W'] = &pt_y[11];
    PHI['S']['S'] = &pt_y[9], PHI['S']['E'] = &pt_y[12], PHI['S']['W'] = &pt_y[13];

    unordered_map<char, char> opposite_azimuth;
    opposite_azimuth['E'] = 'W', opposite_azimuth['W'] = 'E', opposite_azimuth['N'] = 'S', opposite_azimuth['S'] = 'N';

    bool sh = true;

    double XM = pt->MinMaxCoordinate ('x', 'm'), XP = pt->MinMaxCoordinate ('x', 'p');
    double YM = pt->MinMaxCoordinate ('y', 'm'), YP = pt->MinMaxCoordinate ('y', 'p');
    unordered_map<char, double*> m;
    m['E'] = &XP, m['W'] = &XM, m['N'] = &YP, m['S'] = &YM;

    for (const auto &i : {'N', 'S'})
    if (pt->IsBoundary (opposite_azimuth[i]))
    for (const auto &j : {i, 'E', 'W'}) {
      if (pt->EWNS2nd (i, j)) printf ("(i, j) = (%c, %c)\n", i, j), *PT[i][j] = pt->EWNS2nd (i, j), *PHI[i][j] = pt->EWNS2nd (i, j)->Phi (), *m[i] = pt->EWNS2nd (i, j)->Coord ().Value ('y'), sh = false;
      else *PT[i][j] = NULL, *PHI[i][j] = NULL;
    }

    if (IsEqualDouble (XP - XM, YP - YM) || sh)
    for (const auto &i : {'E', 'W'})
    if (pt->IsBoundary (opposite_azimuth[i]))
    for (const auto &j : {i, 'N', 'S'}) {
      if (pt->EWNS2nd (i, j)) *PT[i][j] = pt->EWNS2nd (i, j), *PHI[i][j] = pt->EWNS2nd (i, j)->Phi ();
      else *PT[i][j] = NULL, *PHI[i][j] = NULL;
    }
  }

  printf ("CheckRepresentationFormula (%s)\n", pt->Mark ().c_str ());

  if (pt->Mark ().length () == 1) { // solution cases
    if (this->Condition () == 'N') ydat->F += this->Boundaryvalue ();

    printf ("Boundaryvalue = %23.16e\n", this->Boundaryvalue ());

    for (size_t i = 0; i < 14; i++) if (pt_x[i]) {if (i < 7) result[0] += ent_x[i] * u_ftn_Dirichlet (pt_x[i]); else result[0] += ent_x[i] * u_ftn_Dirichlet (pt_x[i]);}
    for (size_t i = 0; i < 14; i++) if (pt_y[i]) {if (i < 7) result[1] += ent_y[i] * u_ftn_Dirichlet (pt_y[i]); else result[1] += ent_y[i] * u_ftn_Dirichlet (pt_y[i]);}

    for (size_t i = 0; i < 14; i++) if (pt_x[i]) { if (i < 7) printf ("i = %02zu, value = %23.16e, ent_x = %23.16e, u_ftn = %23.16e\n", i, ent_x[i] *  u_ftn_Dirichlet (pt_x[i]), ent_x[i], u_ftn_Dirichlet (pt_x[i]));
                                                   else       printf ("i = %02zu, value = %23.16e, ent_x = %23.16e, u_ftn = %23.16e\n", i, ent_x[i] *  u_ftn_Dirichlet (pt_x[i]), ent_x[i], u_ftn_Dirichlet (pt_x[i])); }
    for (size_t i = 0; i < 14; i++) if (pt_y[i]) { if (i < 7) printf ("i = %02zu, value = %23.16e, ent_y = %23.16e, u_ftn = %23.16e\n", i, ent_y[i] *  u_ftn_Dirichlet (pt_y[i]), ent_y[i], u_ftn_Dirichlet (pt_y[i]));
                                                   else       printf ("i = %02zu, value = %23.16e, ent_y = %23.16e, u_ftn = %23.16e\n", i, ent_y[i] *  u_ftn_Dirichlet (pt_y[i]), ent_y[i], u_ftn_Dirichlet (pt_y[i])); }

    printf ("result[0] = %23.16e\n", result[0]);
    printf ("result[1] = %23.16e\n", result[1]);

    printf ("xdat->F = %23.16e\n", xdat->F);
    printf ("ydat->F = %23.16e\n", ydat->F);
    if (this->Condition () == 'N')  printf ("Representation formula for Nuemann boundary condition = %23.16e\n", result[0] + result[1] - xdat->F - ydat->F);
    else printf ("Representation formula for x = %23.16e\n", result[0] - xdat->F), printf ("Representation formula for y = %23.16e\n", result[1] - ydat->F);

    // printf ("Approximation formula for x = %23.16e\n", sol_x_approx (pt)), printf ("Approximation formula for y = %23.16e\n", sol_y_approx (pt));
    if (this->Condition () == 'N') printf ("Approximation formula for Nuemann boundary condition = %23.16e\n", sol_diff_approx (pt));

    if (this->Condition () == 'N')
    printf ("abs (ux - xValue) = %23.16e\n", fabs (u_ftn_Dirichlet (this->Diff ('x')) * this->Normal ()[0]) - (result[0] - xdat->F)),
    printf ("abs (uy - yValue) = %23.16e\n", fabs (u_ftn_Dirichlet (this->Diff ('y')) * this->Normal ()[1]) - (result[1] - ydat->F + this->Boundaryvalue ()));

    // printf ("this->F = %23.16e\n", this->F ());

    if (*this == 'N')      return fabs (result[0] + result[1] - ydat->F);
    else if (*this == 'D') return max (fabs (result[0]), fabs (result[1] - ydat->F));
    else                   return result[0] + result[1] - ydat->F;
  }

  else if (pt->Mark ().length () == 2) {

    for (size_t i = 1; i < 14; i++) if (pt_x[i]) {if (i < 7) result[0] += ent_x[i] * u_ftn_Dirichlet (pt_x[i]); else result[0] += ent_x[i] *  u_ftn_Dirichlet (pt_x[i])                                  ;}
    for (size_t i = 1; i < 14; i++) if (pt_y[i]) {if (i < 7) result[1] += ent_y[i] * u_ftn_Dirichlet (pt_y[i]); else result[1] += ent_y[i] * (u_ftn_Dirichlet (pt_y[i]) - f_ftn (HeadVelocity (pt_y[i])));}

    for (size_t i = 1; i < 14; i++) if (pt_x[i]) { if (i < 7) printf ("i = %02zu, value = %23.16e, ent_x = %23.16e, u_ftn = %23.16e\n", i, ent_x[i] *  u_ftn_Dirichlet (pt_x[i]),                                   ent_x[i],  u_ftn_Dirichlet (pt_x[i])                                  );
                                                   else       printf ("i = %02zu, value = %23.16e, ent_x = %23.16e, u_ftn = %23.16e\n", i, ent_x[i] *  u_ftn_Dirichlet (pt_x[i]),                                   ent_x[i],  u_ftn_Dirichlet (pt_x[i])                                  ); }
    for (size_t i = 1; i < 14; i++) if (pt_y[i]) { if (i < 7) printf ("i = %02zu, value = %23.16e, ent_y = %23.16e, u_ftn = %23.16e\n", i, ent_y[i] *  u_ftn_Dirichlet (pt_y[i]),                                   ent_y[i],  u_ftn_Dirichlet (pt_y[i])                                  );
                                                   else       printf ("i = %02zu, value = %23.16e, ent_y = %23.16e, u_ftn = %23.16e\n", i, ent_y[i] * (u_ftn_Dirichlet (pt_y[i]) - f_ftn (HeadVelocity (pt_y[i]))), ent_y[i], (u_ftn_Dirichlet (pt_y[i]) - f_ftn (HeadVelocity (pt_y[i])))); }


    printf ("result[0] = %23.16e\n", result[0]);
    printf ("result[1] = %23.16e\n", result[1]);

    printf ("xdat->F = %23.16e\n", xdat->F);
    printf ("ydat->F = %23.16e\n", ydat->F);
    printf ("Representation formula for x = %23.16e\n", result[0] - u_ftn_Dirichlet (pt)), printf ("Representation formula for y = %23.16e\n", result[1] - ydat->F - u_ftn_Dirichlet (pt));

  }

  else if (pt->Mark ().length () == 3) {
    char xy = pt->Mark ()[2];

    result[0] -= 1 / this->MaterialProperty () * (5.0E-1 * this->F () + u_ftn_Dirichlet (this->Phi ()) - 5.0E-1 * u_ftn_Dirichlet (this) / this->Dt () + 5.0E-1 * this->Pre ()->Value () / this->Dt ());
    result[1] -= 1 / this->MaterialProperty () * (5.0E-1 * this->F () - u_ftn_Dirichlet (this->Phi ()) - 5.0E-1 * u_ftn_Dirichlet (this) / this->Dt () + 5.0E-1 * this->Pre ()->Value () / this->Dt ());

    printf ("result[0] = %23.16e\n", result[0]);
    printf ("result[1] = %23.16e\n", result[1]);

    if (xy == 'x') for (size_t i = 1; i < 14; i++) if (pt_x[i]) {if (i < 7) result[0] += ent_x[i] * u_ftn_Dirichlet (pt_x[i]->Diff (xy)); else result[0] += ent_x[i] * u_ftn_Dirichlet (pt_x[i]);}
    if (xy == 'y') for (size_t i = 1; i < 14; i++) if (pt_y[i]) {if (i < 7) result[1] += ent_y[i] * u_ftn_Dirichlet (pt_y[i]->Diff (xy)); else result[1] += ent_y[i] * u_ftn_Dirichlet (pt_y[i]);}

    for (size_t i = 0; i < 14; i++) if (pt_x[i]) { if (i < 7) printf ("i = %02zu, value = %23.16e, ent_x = %23.16e, u_ftn = %23.16e\n", i, ent_x[i] * u_ftn_Dirichlet (pt_x[i]->Diff (xy)), ent_x[i], u_ftn_Dirichlet (pt_x[i]->Diff (xy)));
                                                   else       printf ("i = %02zu, value = %23.16e, ent_x = %23.16e, u_ftn = %23.16e\n", i, ent_x[i] * u_ftn_Dirichlet (pt_x[i]),            ent_x[i], u_ftn_Dirichlet (pt_x[i]));}
    for (size_t i = 0; i < 14; i++) if (pt_y[i]) { if (i < 7) printf ("i = %02zu, value = %23.16e, ent_y = %23.16e, u_ftn = %23.16e\n", i, ent_y[i] * u_ftn_Dirichlet (pt_y[i]->Diff (xy)), ent_y[i], u_ftn_Dirichlet (pt_y[i]->Diff (xy)));
                                                   else       printf ("i = %02zu, value = %23.16e, ent_y = %23.16e, u_ftn = %23.16e\n", i, ent_y[i] * u_ftn_Dirichlet (pt_y[i]),            ent_y[i], u_ftn_Dirichlet (pt_y[i]));}

    printf ("result[0] = %23.16e\n", result[0]);
    printf ("result[1] = %23.16e\n", result[1]);

    printf ("xdat->F = %23.16e\n", xdat->F);
    printf ("ydat->F = %23.16e\n", ydat->F);
    printf ("Representation formula for x = %23.16e\n", result[0] - 5.0E-1 * xdat->F - u_ftn_Dirichlet (pt)), printf ("Representation formula for y = %23.16e\n", result[1] - 5.0E-1 * ydat->F - u_ftn_Dirichlet (pt));

    printf ("Approximation formula for x = %23.16e\n", sol_x_diff_diff_approx (HeadVelocity (pt)));

    if (xy == 'x')      return fabs (result[0] - u_ftn_Dirichlet (pt) - 5.0E-1 * xdat->F);
    else if (xy == 'y') return fabs (result[1] - u_ftn_Dirichlet (pt) - 5.0E-1 * ydat->F);
  }

  exit (0);
}

bool Point::operator== (const int idx) {
  if (this->Index () == idx) return true;
  else return false;
}

bool Point::operator== (const char condition) {
  if (this->Condition () == condition) return true;
  else return false;
}

bool Point::operator== (const Coordinate &src) {
  if (this->Coord () == src) return true;
  else return false;
}


#endif
