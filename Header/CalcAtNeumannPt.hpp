#ifndef CALCNEUMANNPT
#define CALCNEUMANNPT

#include "MatrixProcess.hpp"
/* Neumann 경계조건이 주어진 점에서의 해의 값을 계산하는 모듈 */
double CalcAtNeumannpt (char xy, Point *pt, xData *xdat, yData *ydat) {
  // x-축선에서의 Neumann 경계조건이 주어진 점에서의 해의 값을 계산
  if (xy == 'X' || xy == 'x') {
    // 오른쪽점이 존재하는 경우
    if (pt->EWNS ('E', 'E') != NULL) {
      // 오른쪽점의 경계조건이 내부점 또는 주위의 conductivity가 같은 interface위의 점인 경우
      if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') {
        // Neumann 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산
        CalcNeumannCoef (pt, xdat, ydat);
        // 계산된 계수들과 계산한 해의 값을 이용해서 Neumann 경계조건이 주어진 점에서의 해의 값을 계산
        return - (xdat->Eu * pt->EWNS ('E', 'E')->Value () + xdat->Wu * pt->Boundaryvalue () + xdat->Cphi * pt->Phi ()->Value () + xdat->Ephi * pt->EWNS ('E', 'E')->Phi ()->Value () - xdat->F) / xdat->Cu;
      }
    }
    // 왼쪽점이 존재하는 경우
    if (pt->EWNS ('W', 'W') != NULL) {
      // 왼쪽점의 경계조건이 내부점 또는 주위의 conductivity가 같은 interface위의 점인 경우
      if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') {
        // Neumann 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산
        CalcNeumannCoef (pt, xdat, ydat);
        // 계산된 계수들과 계산한 해의 값을 이용해서 Neumann 경계조건이 주어진 점에서의 해의 값을 계산
        return - (xdat->Eu * pt->Boundaryvalue () + xdat->Wu * pt->EWNS ('W', 'W')->Value () + xdat->Cphi * pt->Phi ()->Value () + xdat->Wphi * pt->EWNS ('W', 'W')->Phi ()->Value () - xdat->F) / xdat->Cu;
      }
    }
  }
  // y-축선에서의 Neumann 경계조건이 주어진 점에서의 해의 값을 계산
  if (xy == 'Y' || xy == 'y') {
    // 위쪽점이 존재하는 경우
    if (pt->EWNS ('N', 'N') != NULL) {
      // 위쪽점의 경계조건이 내부점 또는 주위의 conductivity가 같은 interface위의 점인 경우
      if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') {
        // Neumann 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산
        CalcNeumannCoef (pt, xdat, ydat);
        // 계산된 계수들과 계산한 해의 값을 이용해서 Neumann 경계조건이 주어진 점에서의 해의 값을 계산
        return - (ydat->Nu * pt->EWNS ('N', 'N')->Value () + ydat->Su * pt->Boundaryvalue () + ydat->Cphi * pt->Phi ()->Value () + ydat->Nphi * pt->EWNS ('N', 'N')->Phi ()->Value () - ydat->F) / ydat->Cu;
      }
    }
    // 아래쪽점이 존재하는 경우
    if (pt->EWNS ('S', 'S') != NULL) {
      // 아래쪽점의 경계조건이 내부점 또는 주위의 conductivity가 같은 interface위의 점인 경우
      if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') {
        // Neumann 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산
        CalcNeumannCoef (pt, xdat, ydat);
        // 계산된 계수들과 계산한 해의 값을 이용해서 Neumann 경계조건이 주어진 점에서의 해의 값을 계산
        return - (ydat->Nu * pt->Boundaryvalue () + ydat->Su * pt->EWNS ('S', 'S')->Value () + ydat->Cphi * pt->Phi ()->Value () + ydat->Sphi * pt->EWNS ('S', 'S')->Phi ()->Value () - ydat->F) / ydat->Cu;
      }
    }
  }
  // 잘못된 참조를 한 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcAtNeumannpt");
  return 0.0;
}

/* Neumann 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산하는 모듈 */
void CalcNeumannCoef (Point *pt, xData *xdat, yData *ydat) {
  // nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  int nD = 2;
  // 사용하는 그린함수의 형태는 항상 Dirichlet - Dirichlet
  int bdx = 0, bdy = 0;
  // 국소 x-축선의 시작점과 끝점의 x-좌표와 현재점의 x-좌표
  double xm = pt->MinMaxCoordinate ('x', 'm'), xb = pt->Coord ().Value ('x'), xp = pt->MinMaxCoordinate ('x', 'p');
  // 국소 y-축선의 시작점과 끝점의 y-좌표와 현재점의 y-좌표
  double ym = pt->MinMaxCoordinate ('y', 'm'), yb = pt->Coord ().Value ('y'), yp = pt->MinMaxCoordinate ('y', 'p');
  // 국소 x-축선의 시작점의 conductivity
  double mpx1 = pt->MaterialProperty ();
  // 국소 x-축선의 끝점의 conductivity
  double mpx2 = pt->MaterialProperty ();
  // 국소 y-축선의 시작점의 conductivity
  double mpy1 = pt->MaterialProperty ();
  // 국소 y-축선의 끝점의 conductivity
  double mpy2 = pt->MaterialProperty ();
  // 참조하는 점이 주의의 conductivity가 다른 Interface위의 점인경우
  if (pt->Condition () == 'I') {
    // 오른쪽점이 존재한는 경우, 오른쪽점의 conductivity를 국소 x-축선의 끝점의 conductivity에 저장
    if (pt->EWNS ('E', 'E') != NULL) mpx2 = pt->EWNS ('E', 'E')->MaterialProperty ();
    // 오른쪽의 위쪽점이 존재한는 경우, 오른쪽의 위쪽점의 conductivity를 국소 x-축선의 끝점의 conductivity에 저장
    else if (pt->EWNS ('E', 'N') != NULL) mpx2 = pt->EWNS ('E', 'N')->MaterialProperty ();
    // 오른쪽의 아래쪽점이 존재한는 경우, 오른쪽의 아래쪽점의 conductivity를 국소 x-축선의 끝점의 conductivity에 저장
    else if (pt->EWNS ('E', 'S') != NULL) mpx2 = pt->EWNS ('E', 'S')->MaterialProperty ();
    // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
    else PrintError ("Interface has no East Points!!");
    // 왼쪽점이 존재한는 경우, 왼쪽점의 conductivity를 국소 x-축선의 시작점의 conductivity에 저장
    if (pt->EWNS ('W', 'W') != NULL) mpx1 = pt->EWNS ('W', 'W')->MaterialProperty ();
    // 왼쪽의 위쪽점이 존재한는 경우, 왼쪽의 위쪽점의 conductivity를 국소 x-축선의 시작점의 conductivity에 저장
    else if (pt->EWNS ('W', 'N') != NULL) mpx1 = pt->EWNS ('W', 'N')->MaterialProperty ();
    // 왼쪽의 아래쪽점이 존재한는 경우, 왼쪽의 아래쪽점의 conductivity를 국소 x-축선의 시작점의 conductivity에 저장
    else if (pt->EWNS ('W', 'S') != NULL) mpx1 = pt->EWNS ('W', 'S')->MaterialProperty ();
    // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
    else PrintError ("Interface has no West Points!!");
    // 위쪽점이 존재한는 경우, 위쪽점의 conductivity를 국소 y-축선의 끝점의 conductivity에 저장
    if (pt->EWNS ('N', 'N') != NULL) mpy2 = pt->EWNS ('N', 'N')->MaterialProperty ();
    // 위쪽의 오른쪽점이 존재한는 경우, 위쪽의 오른쪽점의 conductivity를 국소 y-축선의 끝점의 conductivity에 저장
    else if (pt->EWNS ('N', 'E') != NULL) mpy2 = pt->EWNS ('N', 'E')->MaterialProperty ();
    // 위쪽의 왼쪽점이 존재한는 경우, 위쪽의 왼쪽점의 conductivity를 국소 y-축선의 끝점의 conductivity에 저장
    else if (pt->EWNS ('N', 'W') != NULL) mpy2 = pt->EWNS ('N', 'W')->MaterialProperty ();
    // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
    else PrintError ("Interface has no North Points!!");
    // 아래쪽점이 존재한는 경우, 아래쪽점의 conductivity를 국소 y-축선의 시작점의 conductivity에 저장
    if (pt->EWNS ('S', 'S') != NULL) mpy1 = pt->EWNS ('S', 'S')->MaterialProperty ();
    // 아래쪽의 오른쪽점이 존재한는 경우, 아래쪽의 오른쪽점의 conductivity를 국소 y-축선의 시작점의 conductivity에 저장
    else if (pt->EWNS ('S', 'E') != NULL) mpy1 = pt->EWNS ('S', 'E')->MaterialProperty ();
    // 아래쪽의 왼쪽점이 존재한는 경우, 아래쪽의 왼쪽점의 conductivity를 국소 y-축선의 시작점의 conductivity에 저장
    else if (pt->EWNS ('S', 'W') != NULL) mpy1 = pt->EWNS ('S', 'W')->MaterialProperty ();
    // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
    else PrintError ("Interface has no South Points!!");
  }
  // 점이 왼쪽경계점인 경우
  if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') bdx = 1;
  // 점이 오른쪽경계점인 경우
  if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') bdx = 2;
  // 점이 아래쪽경계점인 경우
  if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') bdy = 1;
  // 점이 위쪽경계점인 경우
  if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') bdy = 2;
  //  Neumann 경계조건이 아닌 경우, 에러메시지를 출력한다.
  if (bdx == 0 && bdy == 0 && pt->Condition () != 'N') PrintError ("it does not have Neuamnn boundary condition");
  // x-축선에서의 해의 표현식의 계수들을 0으로 초기화
  xdat->F  = 0.0;
  xdat->Cu = 0.0, xdat->Cphi = 0.0;
  xdat->Eu = 0.0, xdat->ENu = 0.0, xdat->ESu = 0.0;
  xdat->Wu = 0.0, xdat->WNu = 0.0, xdat->WSu = 0.0;
  xdat->Ephi = 0.0, xdat->ENphi = 0.0, xdat->ESphi = 0.0;
  xdat->Wphi = 0.0, xdat->WNphi = 0.0, xdat->WSphi = 0.0;
  xdat->Cf = 0.0; xdat->Ef = 0.0; xdat->Wf = 0.0;
  // y-축선에서의 해의 표현식의 계수들을 0으로 초기화
  ydat->F = 0.0;
  ydat->Cu = 0.0, ydat->Cphi = 0.0;
  ydat->Nu = 0.0, ydat->NEu = 0.0, ydat->NWu = 0.0;
  ydat->Su = 0.0, ydat->SEu = 0.0, ydat->SWu = 0.0;
  ydat->Nphi = 0.0, ydat->NEphi = 0.0, ydat->NWphi = 0.0;
  ydat->Sphi = 0.0, ydat->SEphi = 0.0, ydat->SWphi = 0.0;
  ydat->Cf = 0.0; ydat->Nf = 0.0; ydat->Sf = 0.0;
  // 국소 x-축선에서 u의 현재점에서의 계수
  xdat->Cu = - 1.0;
  // 국소 x-축선에서 u의 국소 x-축선에서 끝점에서의 계수
  xdat->Eu = greens_coefficient_t (xp, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 u의 국소 x-축선에서 시작점에서의 계수
  xdat->Wu = greens_coefficient_t (xm, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 phi의 현재점에서의 계수 (점이 왼쪽 경계점인 경우)
  if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') xdat->Cphi = greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 phi의 현재점에서의 계수 (점이 오른쪽 경계점인 경우)
  if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') xdat->Cphi = greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 phi의 국소 x-축선에서 끝점에서의 계수 (점이 왼쪽 경계점인 경우)
  if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') xdat->Ephi = greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 phi의 국소 x-축선에서 시작점에서의 계수 (점이 오른쪽 경계점인 경우)
  if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') xdat->Wphi = greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 f의 현재점에서의 계수 (점이 왼쪽 경계점인 경우)
  if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') xdat->Cf = greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 f의 현재점에서의 계수 (점이 오른쪽 경계점인 경우)
  if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') xdat->Cf = greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 f의 국소 x-축선에서 끝점에서의 계수 (점이 왼쪽 경계점인 경우)
  if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') xdat->Ef = greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 f의 국소 x-축선에서 시작점에서의 계수 (점이 오른쪽 경계점인 경우)
  if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') xdat->Wf = greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 y-축선에서 u의 현재점에서의 계수
  ydat->Cu = - 1.0;
  // 국소 y-축선에서 u의 국소 y-축선에서 끝점에서의 계수
  ydat->Nu = greens_coefficient_t (yp, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 u의 국소 y-축선에서 시작점에서의 계수
  ydat->Su = greens_coefficient_t (ym, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 3차원 축대칭이고, y좌표가 0인 경우, u의 계수를 0으로 한다.
  if (nD == 3 && IsEqualDouble (yb, 0.0)) ydat->Su = 0.0;
  // 국소 y-축선에서 phi의 현재점에서의 계수 (점이 아래쪽 경계점인 경우)
  if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') ydat->Cphi = - greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 phi의 현재점에서의 계수 (3차원 축대칭이고, y좌표가 0인 경우)
  if (nD == 3 && IsEqualDouble (yb, 0.0)) ydat->Cphi = - (yp * (3.0 / 4.0)) / mpy2;
  // 국소 y-축선에서 phi의 현재점에서의 계수 (점이 위쪽 경계점인 경우)
  if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') ydat->Cphi = - greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 phi의 국소 y-축선에서 끝점에서의 계수 (점이 아래쪽 경계점인 경우)
  if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') ydat->Nphi = - greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 phi의 국소 y-축선에서 끝점에서의 계수 (3차원 축대칭이고, y좌표가 0인 경우)
  if (nD == 3 && IsEqualDouble (yb, 0.0)) ydat->Nphi = - (yp * (1.0 / 4.0)) / mpy2;
  // 국소 y-축선에서 phi의 국소 y-축선에서 시작점에서의 계수 (점이 위쪽 경계점인 경우)
  if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') ydat->Sphi = - greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 f의 현재점에서의 계수 (점이 아래쪽 경계점인 경우)
  if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') ydat->Cf = greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 f의 현재점에서의 계수 (3차원 축대칭이고, y좌표가 0인 경우)
  if (nD == 3 && IsEqualDouble (yb, 0.0)) ydat->Cf = (yp * (3.0 / 4.0)) / mpy2;
  // 국소 y-축선에서 f의 현재점에서의 계수 (점이 위쪽 경계점인 경우)
  if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') ydat->Cf = greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 f의 국소 y-축선에서 끝점에서의 계수 (점이 아래쪽 경계점인 경우)
  if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') ydat->Nf = greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 f의 국소 y-축선에서 끝점에서의 계수 (3차원 축대칭이고, y좌표가 0인 경우)
  if (nD == 3 && IsEqualDouble (yb, 0.0)) ydat->Nf = (yp * (1.0 / 4.0)) / mpy2;
  // 국소 y-축선에서 f의 국소 y-축선에서 시작점에서의 계수 (점이 위쪽 경계점인 경우)
  if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') ydat->Sf = greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 우변의 항 f를 계산
  RightHandSide (pt, xdat, ydat);
}

/* Infinity 경계조건이 주어진 점에서의 해의 값을 계산하는 모듈 */
double CalcAtInfinitypt (char xy, Point *pt, xData *xdat, yData *ydat) {
  // x-축선에서의 Infinity 경계조건이 주어진 점에서의 해의 값을 계산
  if (xy == 'X' || xy == 'x') {
    // 오른쪽점이 존재하는 경우
    if (pt->EWNS ('E', 'E') != NULL) {
      // 오른쪽점의 경계조건이 내부점 또는 주위의 conductivity가 같은 interface위의 점인 경우
      if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') {
        // Infinity 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산
        CalcInfinityCoef (pt, xdat, ydat);
        // 계산된 계수들과 계산한 해의 값을 이용해서 Infinity 경계조건이 주어진 점에서의 해의 값을 계산
        return - (xdat->Eu * pt->EWNS ('E', 'E')->Value () + xdat->Cphi * pt->Phi ()->Value () + xdat->Ephi * pt->EWNS ('E', 'E')->Phi ()->Value () - xdat->F) / xdat->Cu;
      }
    }
    // 왼쪽점이 존재하는 경우
    if (pt->EWNS ('W', 'W') != NULL) {
      // 왼쪽점의 경계조건이 내부점 또는 주위의 conductivity가 같은 interface위의 점인 경우
      if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') {
        // Infinity 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산
        CalcInfinityCoef (pt, xdat, ydat);
        // 계산된 계수들과 계산한 해의 값을 이용해서 Infinity 경계조건이 주어진 점에서의 해의 값을 계산
        return - (xdat->Wu * pt->EWNS ('W', 'W')->Value () + xdat->Cphi * pt->Phi ()->Value () + xdat->Wphi * pt->EWNS ('W', 'W')->Phi ()->Value () - xdat->F) / xdat->Cu;
      }
    }
  }
  // y-축선에서의 Infinity 경계조건이 주어진 점에서의 해의 값을 계산
  if (xy == 'Y' || xy == 'y') {
    // 위쪽점이 존재하는 경우
    if (pt->EWNS ('N', 'N') != NULL) {
      // 위쪽점의 경계조건이 내부점 또는 주위의 conductivity가 같은 interface위의 점인 경우
      if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') {
        // Infinity 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산
        CalcInfinityCoef (pt, xdat, ydat);
        // 계산된 계수들과 계산한 해의 값을 이용해서 Infinity 경계조건이 주어진 점에서의 해의 값을 계산
        return - (ydat->Nu * pt->EWNS ('N', 'N')->Value () + ydat->Cphi * pt->Phi ()->Value () + ydat->Nphi * pt->EWNS ('N', 'N')->Phi ()->Value () - ydat->F) / ydat->Cu;
      }
    }
    // 아래쪽점이 존재하는 경우
    if (pt->EWNS ('S', 'S') != NULL) {
      // 아래쪽점의 경계조건이 내부점 또는 주위의 conductivity가 같은 interface위의 점인 경우
      if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') {
        // Infinity 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산
        CalcInfinityCoef (pt, xdat, ydat);
        // 계산된 계수들과 계산한 해의 값을 이용해서 Infinity 경계조건이 주어진 점에서의 해의 값을 계산
        return - (ydat->Su * pt->EWNS ('S', 'S')->Value () + ydat->Cphi * pt->Phi ()->Value () + ydat->Sphi * pt->EWNS ('S', 'S')->Phi ()->Value () - ydat->F) / ydat->Cu;
      }
    }
  }
  // 잘못된 참조를 한 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcAtNeumannpt");
  return 0.0;
}

/* Infinity 경계조건이 주어진 점에서의 해의 값을 계산하기 위해서 해의 표현식의 계수를 계산하는 모듈 */
void CalcInfinityCoef (Point *pt, xData *xdat, yData *ydat) {
  // nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  int nD = 2;
  // 사용하는 그린함수의 형태는 항상 Dirichlet - Dirichlet
  int bdx = 0, bdy = 0;
  // 국소 x-축선의 시작점과 끝점의 x-좌표와 현재점의 x-좌표
  double xm = pt->MinMaxCoordinate ('x', 'm'), xb = pt->Coord ().Value ('x'), xp = pt->MinMaxCoordinate ('x', 'p');
  // 국소 y-축선의 시작점과 끝점의 y-좌표와 현재점의 y-좌표
  double ym = pt->MinMaxCoordinate ('y', 'm'), yb = pt->Coord ().Value ('y'), yp = pt->MinMaxCoordinate ('y', 'p');
  // 국소 x-축선의 시작점의 conductivity
  double mpx1 = pt->MaterialProperty ();
  // 국소 x-축선의 끝점의 conductivity
  double mpx2 = pt->MaterialProperty ();
  // 국소 y-축선의 시작점의 conductivity
  double mpy1 = pt->MaterialProperty ();
  // 국소 y-축선의 끝점의 conductivity
  double mpy2 = pt->MaterialProperty ();
  // 참조하는 점이 주의의 conductivity가 다른 Interface위의 점인경우
  if (pt->Condition () == 'I') {
    // 오른쪽점이 존재한는 경우, 오른쪽점의 conductivity를 국소 x-축선의 끝점의 conductivity에 저장
    if (pt->EWNS ('E', 'E') != NULL) mpx2 = pt->EWNS ('E', 'E')->MaterialProperty ();
    // 오른쪽의 위쪽점이 존재한는 경우, 오른쪽의 위쪽점의 conductivity를 국소 x-축선의 끝점의 conductivity에 저장
    else if (pt->EWNS ('E', 'N') != NULL) mpx2 = pt->EWNS ('E', 'N')->MaterialProperty ();
    // 오른쪽의 아래쪽점이 존재한는 경우, 오른쪽의 아래쪽점의 conductivity를 국소 x-축선의 끝점의 conductivity에 저장
    else if (pt->EWNS ('E', 'S') != NULL) mpx2 = pt->EWNS ('E', 'S')->MaterialProperty ();
    // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
    else PrintError ("Interface has no East Points!!");
    // 왼쪽점이 존재한는 경우, 왼쪽점의 conductivity를 국소 x-축선의 시작점의 conductivity에 저장
    if (pt->EWNS ('W', 'W') != NULL) mpx1 = pt->EWNS ('W', 'W')->MaterialProperty ();
    // 왼쪽의 위쪽점이 존재한는 경우, 왼쪽의 위쪽점의 conductivity를 국소 x-축선의 시작점의 conductivity에 저장
    else if (pt->EWNS ('W', 'N') != NULL) mpx1 = pt->EWNS ('W', 'N')->MaterialProperty ();
    // 왼쪽의 아래쪽점이 존재한는 경우, 왼쪽의 아래쪽점의 conductivity를 국소 x-축선의 시작점의 conductivity에 저장
    else if (pt->EWNS ('W', 'S') != NULL) mpx1 = pt->EWNS ('W', 'S')->MaterialProperty ();
    // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
    else PrintError ("Interface has no West Points!!");
    // 위쪽점이 존재한는 경우, 위쪽점의 conductivity를 국소 y-축선의 끝점의 conductivity에 저장
    if (pt->EWNS ('N', 'N') != NULL) mpy2 = pt->EWNS ('N', 'N')->MaterialProperty ();
    // 위쪽의 오른쪽점이 존재한는 경우, 위쪽의 오른쪽점의 conductivity를 국소 y-축선의 끝점의 conductivity에 저장
    else if (pt->EWNS ('N', 'E') != NULL) mpy2 = pt->EWNS ('N', 'E')->MaterialProperty ();
    // 위쪽의 왼쪽점이 존재한는 경우, 위쪽의 왼쪽점의 conductivity를 국소 y-축선의 끝점의 conductivity에 저장
    else if (pt->EWNS ('N', 'W') != NULL) mpy2 = pt->EWNS ('N', 'W')->MaterialProperty ();
    // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
    else PrintError ("Interface has no North Points!!");
    // 아래쪽점이 존재한는 경우, 아래쪽점의 conductivity를 국소 y-축선의 시작점의 conductivity에 저장
    if (pt->EWNS ('S', 'S') != NULL) mpy1 = pt->EWNS ('S', 'S')->MaterialProperty ();
    // 아래쪽의 오른쪽점이 존재한는 경우, 아래쪽의 오른쪽점의 conductivity를 국소 y-축선의 시작점의 conductivity에 저장
    else if (pt->EWNS ('S', 'E') != NULL) mpy1 = pt->EWNS ('S', 'E')->MaterialProperty ();
    // 아래쪽의 왼쪽점이 존재한는 경우, 아래쪽의 왼쪽점의 conductivity를 국소 y-축선의 시작점의 conductivity에 저장
    else if (pt->EWNS ('S', 'W') != NULL) mpy1 = pt->EWNS ('S', 'W')->MaterialProperty ();
    // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
    else PrintError ("Interface has no South Points!!");
  }
  // 점이 왼쪽경계점인 경우
  if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') bdx = 3;
  // 점이 오른쪽경계점인 경우
  if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') bdx = 4;
  // 점이 아래쪽경계점인 경우
  if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') bdy = 3;
  // 점이 위쪽경계점인 경우
  if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') bdy = 4;
  //  Neumann 경계조건이 아닌 경우, 에러메시지를 출력한다.
  if (bdx == 0 && bdy == 0 && pt->Condition () != 'F') PrintError ("it does not have Infinity boundary condition");
  // x-축선에서의 해의 표현식의 계수들을 0으로 초기화
  xdat->F  = 0.0;
  xdat->Cu = 0.0, xdat->Cphi = 0.0;
  xdat->Eu = 0.0, xdat->ENu = 0.0, xdat->ESu = 0.0;
  xdat->Wu = 0.0, xdat->WNu = 0.0, xdat->WSu = 0.0;
  xdat->Ephi = 0.0, xdat->ENphi = 0.0, xdat->ESphi = 0.0;
  xdat->Wphi = 0.0, xdat->WNphi = 0.0, xdat->WSphi = 0.0;
  xdat->Cf = 0.0; xdat->Ef = 0.0; xdat->Wf = 0.0;
  // y-축선에서의 해의 표현식의 계수들을 0으로 초기화
  ydat->F = 0.0;
  ydat->Cu = 0.0, ydat->Cphi = 0.0;
  ydat->Nu = 0.0, ydat->NEu = 0.0, ydat->NWu = 0.0;
  ydat->Su = 0.0, ydat->SEu = 0.0, ydat->SWu = 0.0;
  ydat->Nphi = 0.0, ydat->NEphi = 0.0, ydat->NWphi = 0.0;
  ydat->Sphi = 0.0, ydat->SEphi = 0.0, ydat->SWphi = 0.0;
  ydat->Cf = 0.0; ydat->Nf = 0.0; ydat->Sf = 0.0;
  // 국소 x-축선에서 u의 현재점에서의 계수
  xdat->Cu = - 1.0;
  // 국소 x-축선에서 u의 국소 x-축선에서 끝점에서의 계수
  xdat->Eu = greens_coefficient_t (xp, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 u의 국소 x-축선에서 시작점에서의 계수
  xdat->Wu = greens_coefficient_t (xm, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 phi의 현재점에서의 계수 (점이 왼쪽 경계점인 경우)
  if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') xdat->Cphi = greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 phi의 현재점에서의 계수 (점이 오른쪽 경계점인 경우)
  if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') xdat->Cphi = greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 phi의 국소 x-축선에서 끝점에서의 계수 (점이 왼쪽 경계점인 경우)
  if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') xdat->Ephi = greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 phi의 국소 x-축선에서 시작점에서의 계수 (점이 오른쪽 경계점인 경우)
  if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') xdat->Wphi = greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 f의 현재점에서의 계수 (점이 왼쪽 경계점인 경우)
  if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') xdat->Cf = greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 f의 현재점에서의 계수 (점이 오른쪽 경계점인 경우)
  if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') xdat->Cf = greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 f의 국소 x-축선에서 끝점에서의 계수 (점이 왼쪽 경계점인 경우)
  if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->Condition () == 'C' || pt->EWNS ('E', 'E')->Condition () == 'M') xdat->Ef = greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // 국소 x-축선에서 f의 국소 x-축선에서 시작점에서의 계수 (점이 오른쪽 경계점인 경우)
  if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->Condition () == 'C' || pt->EWNS ('W', 'W')->Condition () == 'M') xdat->Wf = greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  // x-축선의 양끝점중 한 개의 점의 경계조건이 Infinity 혹은 Singularity일 때, 정확한 적분을 계산
  ExactIntegration ('x', xdat, ydat, nD, bdx, bdy, xm, xb, xp, ym, yb, yp);
  // 국소 y-축선에서 u의 현재점에서의 계수
  ydat->Cu = - 1.0;
  // 국소 y-축선에서 u의 국소 y-축선에서 끝점에서의 계수
  ydat->Nu = greens_coefficient_t (yp, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 u의 국소 y-축선에서 시작점에서의 계수
  ydat->Su = greens_coefficient_t (ym, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 phi의 현재점에서의 계수 (점이 아래쪽 경계점인 경우)
  if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') ydat->Cphi = - greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 phi의 현재점에서의 계수 (점이 위쪽 경계점인 경우)
  if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') ydat->Cphi = - greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 phi의 국소 y-축선에서 끝점에서의 계수 (점이 아래쪽 경계점인 경우)
  if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') ydat->Nphi = - greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 phi의 국소 y-축선에서 시작점에서의 계수 (점이 위쪽 경계점인 경우)
  if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') ydat->Sphi = - greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 f의 현재점에서의 계수 (점이 아래쪽 경계점인 경우)
  if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') ydat->Cf = greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 f의 현재점에서의 계수 (점이 위쪽 경계점인 경우)
  if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') ydat->Cf = greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 f의 국소 y-축선에서 끝점에서의 계수 (점이 아래쪽 경계점인 경우)
  if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->Condition () == 'C' || pt->EWNS ('N', 'N')->Condition () == 'M') ydat->Nf = greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // 국소 y-축선에서 f의 국소 y-축선에서 시작점에서의 계수 (점이 위쪽 경계점인 경우)
  if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->Condition () == 'C' || pt->EWNS ('S', 'S')->Condition () == 'M') ydat->Sf = greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  // y-축선의 양끝점중 한 개의 점의 경계조건이 Infinity 혹은 Singularity일 때, 정확한 적분을 계산
  ExactIntegration ('y', xdat, ydat, nD, bdx, bdy, xm, xb, xp, ym, yb, yp);
  // 우변의 항 f를 계산
  RightHandSide (pt, xdat, ydat);
}

#endif
