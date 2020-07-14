#ifndef CALCREPRESENSOL_H
#define CALCREPRESENSOL_H

#include "Point.hpp"

/* 해의 표현식의 계수를 계산 */
void CalcRepresenCoef (Point *pt, xData *xdat, yData *ydat, bool is_sol) {
  /*
  nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  2차원의 무한대에서의 점근적함수는 u (x, y) = log (x^2 + y^2)
  3차원축대칭의 무한대에서의 점근적함수는 u (r, z) = 1 / (r^2 + z^2)
  */
  int nD = 2;
  /*
  사용하는 그린함수의 형태
  bdx: x-축선
  bdy: y-축선
  0: Dirichlet - Dirichlet
  1: Neumann - Dirichlet
  2: Dirichlet - Neumann
  3: Infinity - Dirichlet
  4: Dirichlet - Infinity
  */
  int bdx = 0, bdy = 0;
  // 국소 x-축선의 시작점과 끝점의 x-좌표와 현재점의 x-좌표
  double xm = pt->Pressure ()->MinMaxCoordinate ('x', 'm'), xb = pt->Pressure ()->Coord ().Value ('x'), xp = pt->Pressure ()->MinMaxCoordinate ('x', 'p');
  // 국소 y-축선의 시작점과 끝점의 y-좌표와 현재점의 y-좌표
  double ym = pt->Pressure ()->MinMaxCoordinate ('y', 'm'), yb = pt->Pressure ()->Coord ().Value ('y'), yp = pt->Pressure ()->MinMaxCoordinate ('y', 'p');
  // 국소 x-축선의 시작점의 conductivity
  double mpx1 = pt->Pressure ()->MaterialProperty ();
  // 국소 x-축선의 끝점의 conductivity
  double mpx2 = pt->Pressure ()->MaterialProperty ();
  // 국소 y-축선의 시작점의 conductivity
  double mpy1 = pt->Pressure ()->MaterialProperty ();
  // 국소 y-축선의 끝점의 conductivity
  double mpy2 = pt->Pressure ()->MaterialProperty ();
  // 참조하는 점이 주의의 conductivity가 다른 Interface위의 점인경우
  if (pt->Condition () == 'I') SettingMaterialProperties (pt->Pressure (), &mpx1, &mpx2, &mpy1, &mpy2);
  //
  // SettingGreenShape (pt->Pressure (), &bdx, &bdy);
  //
  InitializationCoef (xdat, ydat, "XY");
  //
  if (pt->Mark ().length () == 1) {
    if      (pt->Condition () == 'N') SettingNeumannCoefficient   (pt, xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1, mpy2, is_sol);
    else if (pt->Condition () == 'D') SettingDirichletCoefficient (pt, xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1, mpy2, nD, is_sol);
    else                              SettingCoefficient          (pt, xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1, mpy2, nD, is_sol);
  }
  else if (pt->Mark ().length () == 2) SettingDiffCoefficient     (pt, xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1, mpy2);
  else if (pt->Mark ().length () == 3) SettingDiffDiffCoefficient (pt, xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1, mpy2);

  // Interface와 x-축선이 만나서 생기는 점인 경우, y-축선에서의 표현식을 만들지 않음 (phi를 계산하지 않기 때문에 식이 1개만 필요)
  if (pt->Pressure ()->Condition () == 'I' && pt->Pressure ()->Axis ('y') < 0) InitializationCoef (xdat, ydat, "Y");
  // Interface와 y-축선이 만나서 생기는 점인 경우, x-축선에서의 표현식을 만들지 않음 (phi를 계산하지 않기 때문에 식이 1개만 필요)
  if (pt->Pressure ()->Condition () == 'I' && pt->Pressure ()->Axis ('x') < 0) InitializationCoef (xdat, ydat, "X");

  if (pt->Mark ().length () == 3) ApproximateDiffDiffPt  (HeadVelocity (pt), xdat, ydat, xm, xb, xp, ym, yb, yp, mpx1, mpx2, mpy1, mpy2, is_sol);
  else                            ApproximateInterfacePt (HeadVelocity (pt), xdat, ydat, xm, xb, xp, ym, yb, yp, mpx1, mpx2, mpy1, mpy2, is_sol);

  if (pt->Mark ().length () == 1) {
      if (pt->Condition() != 'D') RightHandSide(pt, xdat, ydat);

    if (::SolverType) {
      if      (pt->Condition () == 'N') SettingNeumannTimeIntegrationRightHand (pt, xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1, mpy2, is_sol);
      else if (pt->Condition () != 'D') SettingTimeIntegrationRightHand (pt, xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1, mpy2, is_sol);
    }
  }
  else if (pt->Mark ().length () == 2) {
    RightHandSide (HeadVelocity (pt), xdat, ydat);
    if (::SolverType) SettingDiffTimeIntegrationRightHand (HeadVelocity (pt), xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1, mpy2);
  }
  else if (pt->Mark ().length () == 3) {
    RightHandSide (HeadVelocity (pt), xdat, ydat);
    if (::SolverType) SettingDiffDiffTimeIntegrationRightHand (HeadVelocity (pt), xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1, mpy2);
  }
}

/* 우변의 항 f를 계산 */
void RightHandSide (Point *pt, xData *xdat, yData *ydat) {
  // 현재점에서의 우변의 f를 저장
  xdat->F += xdat->Cf * pt->F ();
  ydat->F += ydat->Cf * pt->F ();
  // 오른쪽점이 존재하는 경우, 오른쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('E', 'E') != NULL) xdat->F += xdat->Ef * pt->EWNS ('E', 'E')->F ();
  // 왼쪽점이 존재하는 경우, 왼쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('W', 'W') != NULL) xdat->F += xdat->Wf * pt->EWNS ('W', 'W')->F ();
  // 위쪽점이 존재하는 경우, 위쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('N', 'N') != NULL) ydat->F += ydat->Nf * pt->EWNS ('N', 'N')->F ();
  // 아래쪽점이 존재하는 경우, 아래쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('S', 'S') != NULL) ydat->F += ydat->Sf * pt->EWNS ('S', 'S')->F ();
  // 오른쪽의 위쪽점이 존재하는 경우, 오른쪽의 위쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('E', 'N') != NULL) xdat->F += xdat->ENf * pt->EWNS ('E', 'N')->F ();
  // 오른쪽의 아래쪽점이 존재하는 경우, 오른쪽의 아래쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('E', 'S') != NULL) xdat->F += xdat->ESf * pt->EWNS ('E', 'S')->F ();
  // 왼쪽의 위쪽점이 존재하는 경우, 왼쪽의 위쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('W', 'N') != NULL) xdat->F += xdat->WNf * pt->EWNS ('W', 'N')->F ();
  // 왼쪽의 아래쪽점이 존재하는 경우, 왼쪽의 아래쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('W', 'S') != NULL) xdat->F += xdat->WSf * pt->EWNS ('W', 'S')->F ();
  // 위쪽의 오른쪽점이 존재하는 경우, 위쪽의 오른쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('N', 'E') != NULL) ydat->F += ydat->NEf * pt->EWNS ('N', 'E')->F ();
  // 위쪽의 왼쪽점이 존재하는 경우, 위쪽의 왼쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('N', 'W') != NULL) ydat->F += ydat->NWf * pt->EWNS ('N', 'W')->F ();
  // 아래쪽의 오른쪽점이 존재하는 경우, 아래쪽의 오른쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('S', 'E') != NULL) ydat->F += ydat->SEf * pt->EWNS ('S', 'E')->F ();
  // 아래쪽의 왼쪽점이 존재하는 경우, 아래쪽의 왼쪽점에서의 우변의 f를 저장
  if (pt->EWNS ('S', 'W') != NULL) ydat->F += ydat->SWf * pt->EWNS ('S', 'W')->F ();
  // f가 우변이기 때문에 -1을 곱함
  xdat->F = - 5.0E-1 * xdat->F;
  ydat->F = - 5.0E-1 * ydat->F;
}

void RightHandSideDirichlet (Point *pt, Point *calc_pt, xData *xdat, yData *ydat) {

  unordered_map<char, double> val;
  unordered_map<char, char> opposite_azimuth;
  val['E'] = xdat->Ef, val['W'] = xdat->Wf;
  val['N'] = ydat->Nf, val['S'] = ydat->Sf;
  opposite_azimuth['E'] = 'W', opposite_azimuth['W'] = 'E';
  opposite_azimuth['N'] = 'S', opposite_azimuth['S'] = 'N';

  xdat->F += xdat->Cf * calc_pt->F ();
  ydat->F += ydat->Cf * calc_pt->F ();
  for (const auto &i : {'E', 'W'})
  if (calc_pt->EWNS (i, i))
  if (calc_pt->EWNS (i, i)->EWNS (i, i))
  if (calc_pt->EWNS (i, i)->EWNS (i, i) == pt) {
    xdat->F += val[i] * pt->F ();
    if (calc_pt->EWNS (opposite_azimuth[i], opposite_azimuth[i])) xdat->F += val[opposite_azimuth[i]] * calc_pt->EWNS (opposite_azimuth[i], opposite_azimuth[i])->F ();
    else PrintError ("RightHandSideDirichlet");
  }

  for (const auto &i : {'N', 'S'})
  if (calc_pt->EWNS (i, i))
  if (calc_pt->EWNS (i, i)->EWNS (i, i))
  if (calc_pt->EWNS (i, i)->EWNS (i, i) == pt) {
    ydat->F += val[i] * pt->F ();
    if (calc_pt->EWNS (opposite_azimuth[i], opposite_azimuth[i])) ydat->F += val[opposite_azimuth[i]] * calc_pt->EWNS (opposite_azimuth[i], opposite_azimuth[i])->F ();
    else PrintError ("RightHandSideDirichlet");
  }
  // f가 우변이기 때문에 -1을 곱함
  xdat->F = - 5.0E-1 * xdat->F;
  ydat->F = - 5.0E-1 * ydat->F;

}

void RightHandSideNeumann (Point *pt, xData *xdat, yData *ydat) {
  char errorMassage[256];
  unordered_map<char, double> eVal;
  unordered_map<char, double> wVal;
  unordered_map<char, double> nVal;
  unordered_map<char, double> sVal;
  unordered_map<char, char> opposite_azimuth;
  unordered_map<char, unordered_map<char, double>> val;
  opposite_azimuth['E'] = 'W', opposite_azimuth['W'] = 'E';
  opposite_azimuth['N'] = 'S', opposite_azimuth['S'] = 'N';
  val['E'] = eVal, val['W'] = wVal;
  val['N'] = nVal, val['S'] = sVal;
  val['E']['E'] = xdat->Ef, val['E']['N'] = xdat->ENf, val['E']['S'] = xdat->ESf;
  val['W']['W'] = xdat->Wf, val['W']['N'] = xdat->WNf, val['W']['S'] = xdat->WSf;
  val['N']['N'] = ydat->Nf, val['N']['E'] = ydat->NEf, val['N']['W'] = ydat->NWf;
  val['S']['S'] = ydat->Sf, val['S']['E'] = ydat->SEf, val['S']['W'] = ydat->SWf;

  xdat->F += xdat->Cf * pt->F ();
  ydat->F += ydat->Cf * pt->F ();

  bool idx = false;

  for (auto &i : {'N', 'S'})
  if (pt->IsBoundary (opposite_azimuth[i]))
  for (auto &j : {i, 'E', 'W'}) {
    if (pt->EWNS2nd (i, j)) ydat->F += val[i][j] * pt->EWNS2nd (i, j)->F ();
    else if (i == j) idx = true;
    else if (idx) sprintf (errorMassage, "RightHandSideNeumann, pt->EWNS2nd (%c, %c) does not exist.", i, j), PrintError (errorMassage);
  }

  // f가 우변이기 때문에 -1을 곱함
  xdat->F = - xdat->F;
  ydat->F = - ydat->F;
}

/* 경계근처에서의 해의 표현식을 처리 */
void TransposeBoundaryData (Point *pt, xData *xdat, yData *ydat, bool is_sol) {
  if      (pt->Condition () == 'N') TransposeNeumannBoundaryData (pt, xdat, ydat);
  else if (pt->Condition () != 'D') TransposeOtherBoundaryData (pt, xdat, ydat);
}

/* 경계에서의 phi값을 부여 */
void AssignPhivalue (AxialData *adat, Point *pt) {
  for (size_t i = 0; i < adat->Pts_Num (); i++) AssignPhi (&pt[i]);
}

int AssignPhi (Point *pt) {
  char azimu[4] = {'E', 'W', 'N', 'S'};
  char bound[3] = {'D', 'N', 'F'};
  char inner[3] = {'C', 'I', 'M'};

  for (const auto &i : bound) if (pt->Condition () == i)
  for (const auto &j : azimu) if (pt->EWNS (j, j))
  for (const auto &k : inner) if (pt->EWNS (j, j)->Condition () == k) {
    pt->Phi ()->SetValue (pt->EWNS (j, j)->Phi ()->Value ());
    return 0;
  }

  return 1;
}

/* 경계와 Interface가 만나는지 확인하는 모듈 */
bool SetInterfaceBoundary (Point *pt, char EWNS) {
  // 현재점이 Interface위의 점인 경우 거짓값을 return
  if (pt->Condition () != 'I' && pt->Condition () != 'M') return false;
  // 오른쪽을 참조
  if (EWNS == 'E' || EWNS == 'e') {
    // 오른쪽점이 경계점인 경우
    if (pt->EWNS ('E', 'E') == pt) {
      // 위쪽점이 존재하고 위쪽점의 오른쪽점이 존재하고 위쪽점의 오른쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->EWNS ('E', 'E') != NULL) if (pt->EWNS ('N', 'N')->EWNS ('E', 'E')->Condition () == 'N') return true;
      // 아래쪽점이 존재하고 아래쪽점의 오른쪽점이 존재하고 아래쪽점의 오른쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->EWNS ('E', 'E') != NULL) if (pt->EWNS ('S', 'S')->EWNS ('E', 'E')->Condition () == 'N') return true;
      // 위쪽의 오른쪽점이 존재하고 위쪽의 오른쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('N', 'E') != NULL) if (pt->EWNS ('N', 'E')->Condition () == 'N') return true;
      // 아래쪽의 오른쪽점이 존재하고 아래쪽의 오른쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('S', 'E') != NULL) if (pt->EWNS ('S', 'E')->Condition () == 'N') return true;
    }
  }
  // 왼쪽을 참조
  if (EWNS == 'W' || EWNS == 'w') {
    // 왼쪽점이 경계점인 경우
    if (pt->EWNS ('W', 'W') == pt) {
      // 위쪽점이 존재하고 위쪽점의 왼쪽점이 존재하고 위쪽점의 왼쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('N', 'N') != NULL) if (pt->EWNS ('N', 'N')->EWNS ('W', 'W') != NULL) if (pt->EWNS ('N', 'N')->EWNS ('W', 'W')->Condition () == 'N') return true;
      // 아래쪽점이 존재하고 아래쪽점의 왼쪽점이 존재하고 아래쪽점의 왼쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('S', 'S') != NULL) if (pt->EWNS ('S', 'S')->EWNS ('W', 'W') != NULL) if (pt->EWNS ('S', 'S')->EWNS ('W', 'W')->Condition () == 'N') return true;
      // 위쪽의 왼쪽점이 존재하고 위쪽의 왼쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('N', 'W') != NULL) if (pt->EWNS ('N', 'W')->Condition () == 'N') return true;
      // 아래쪽의 왼쪽점이 존재하고 아래쪽의 왼쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('S', 'W') != NULL) if (pt->EWNS ('S', 'W')->Condition () == 'N') return true;
    }
  }
  // 위쪽을 참조
  if (EWNS == 'N' || EWNS == 'n') {
    // 위쪽점이 경계점인 경우
    if (pt->EWNS ('N', 'N') == pt) {
      // 오른쪽점이 존재하고 오른쪽점의 위쪽점이 존재하고 오른쪽점의 위쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->EWNS ('N', 'N') != NULL) if (pt->EWNS ('E', 'E')->EWNS ('N', 'N')->Condition () == 'N') return true;
      // 왼쪽점이 존재하고 왼쪽점의 위쪽점이 존재하고 왼쪽점의 위쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->EWNS ('N', 'N') != NULL) if (pt->EWNS ('W', 'W')->EWNS ('N', 'N')->Condition () == 'N') return true;
      // 오른쪽의 위쪽점이 존재하고 오른쪽의 위쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('E', 'N') != NULL) if (pt->EWNS ('E', 'N')->Condition () == 'N') return true;
      // 왼쪽의 위쪽점이 존재하고 왼쪽의 위쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('W', 'N') != NULL) if (pt->EWNS ('W', 'N')->Condition () == 'N') return true;
    }
  }
  // 아래쪽을 참조
  if (EWNS == 'S' || EWNS == 's') {
    // 아래쪽점이 경계점인 경우
    if (pt->EWNS ('S', 'S') == pt) {
      // 오른쪽점이 존재하고 오른쪽점의 아래쪽점이 존재하고 오른쪽점의 아래쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('E', 'E') != NULL) if (pt->EWNS ('E', 'E')->EWNS ('S', 'S') != NULL) if (pt->EWNS ('E', 'E')->EWNS ('S', 'S')->Condition () == 'N') return true;
      // 왼쪽점이 존재하고 왼쪽점의 아래쪽점이 존재하고 왼쪽점의 아래쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('W', 'W') != NULL) if (pt->EWNS ('W', 'W')->EWNS ('S', 'S') != NULL) if (pt->EWNS ('W', 'W')->EWNS ('S', 'S')->Condition () == 'N') return true;
      // 오른쪽의 아래쪽점이 존재하고 오른쪽의 아래쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('E', 'S') != NULL) if (pt->EWNS ('E', 'S')->Condition () == 'N') return true;
      // 왼쪽의 아래쪽점이 존재하고 왼쪽의 아래쪽점의 경계조건이 Neumann 조건이 주어진 경우, 참값을 return
      if (pt->EWNS ('W', 'S') != NULL) if (pt->EWNS ('W', 'S')->Condition () == 'N') return true;
    }
  }
  return 0;
}

/* 축선의 양끝점중 한 개의 점의 경계조건이 Infinity 혹은 Singularity일 때, 정확한 적분을 계산하는 모듈 */
void ExactIntegration(char xy, xData *xdat, yData *ydat, int nD, int bdx, int bdy, double xm, double xb, double xp, double ym, double yb, double yp) {
  /*
  AsymptoticType: asymptotic behavior의 타입
  1: u (x, y) ~ log (sqrt (x^2 + y^2)) + c * 1

  축대칭의 경우:
  1: u (z, r) ~ c / sqrt (z^2 + r^2)
  */
  int AsymptoticType = 2;
  // x-축선을 참조하는 경우
  if (xy == 'X' || xy == 'x') {
    // x-축선 위에서 Infinity - Dirichlet 그린함수를 사용하는 경우
    if (bdx == 3) {
      // 2차원문제인 경우
      if (nD == 2) {
        // 점근적함수가 log (sqrt (x^2 + y^2)) + c * 1인 경우
        if (AsymptoticType == 1) {
          // yb가 0인 경우
          if (IsEqualDouble (yb, ZeroValue)) {
            // x-축선위에서 점근적함수의 적분 계산
            xdat->F += (xb*(2.0/3.0))/xp-2.0/3.0;

            // x-축선위에서 축선의 오른쪽점에서의 u 계수인 점근적함수의 적분 계산
            xdat->Eu += (xb*2.0)/xp-2.0;

            // x-축선위에서  x < xi인 경우의 그린함수와 phi1의 적분의 계산
            xdat->F += 1.0/(xb*xb*xb)*(xp*xp)*(xb-xp)*(1.0/3.0);

            // x-축선위에서  x < xi인 경우의 그린함수와 phi2의 적분의 계산 (f의 계수)
            xdat->F += 0.0;

            // x-축선위에서  x < xi인 경우의 그린함수와 phi2의 적분의 계산 (Eu의 계수)
            xdat->Eu += 0.0;
          }
          // yb가 0이 아닌 경우
          else {
            // x-축선위에서 점근적함수의 적분 계산
            xdat->F += xp*1.0/(yb*yb*yb)*(xb-xp)*(yb*2.0-xp*3.141592653589793+xp*atan(xp/yb)*2.0);

            // x-축선위에서 축선의 오른쪽점에서의 u 계수인 점근적함수의 적분 계산
            xdat->Eu += (xb*2.0)/xp-2.0;

            // x-축선위에서  x < xi인 경우의 그린함수와 phi1의 적분의 계산
            xdat->F += -((xp*xp)*1.0/(yb*yb*yb)*(xb-xp)*((xb*xb*xb)*atan(xb/yb)*2.0-(xb*xb*xb)*3.141592653589793+(xb*xb)*yb*2.0+yb*yb*yb-xb*(yb*yb)*3.141592653589793+xb*(yb*yb)*atan(xb/yb)*2.0))/(xb*(xb*xb+yb*yb));

            // x-축선위에서  x < xi인 경우의 그린함수와 phi2의 적분의 계산 (f의 계수)
            xdat->F += 0.0;

            // x-축선위에서  x < xi인 경우의 그린함수와 phi2의 적분의 계산 (Eu의 계수)
            xdat->Eu += 0.0;
          }
        }

        // 점근적함수가 c * 1.0E0 / sqrt (x ^ 2 + y ^ 2)인 경우
        if (AsymptoticType == 2) {
          // yb가 0인 경우
          if (IsEqualDouble (yb, ZeroValue)) {
            // x-축선위에서 점근적함수의 적분 계산
            xdat->F += 0.0;

            // x-축선위에서 축선의 오른쪽점에서의 u 계수인 점근적함수의 적분 계산
            xdat->Eu += (xb*(3.0/2.0)-xp*(3.0/2.0))/xp;

            // x-축선위에서  x < xi인 경우의 그린함수와 phi1의 적분의 계산
            xdat->F += 0.0;

            // x-축선위에서  x < xi인 경우의 그린함수와 phi2의 적분의 계산 (f의 계수)
            xdat->F += 0.0;

            // x-축선위에서  x < xi인 경우의 그린함수와 phi2의 적분의 계산 (Eu의 계수)
            xdat->Eu += 1.0/(xb*xb*xb*xb)*(xp*xp*xp)*(xb-xp)*(-1.0/2.0);
          }
          // yb가 0이 아닌 경우
          else {
            // x-축선위에서 점근적함수의 적분 계산
            xdat->F += 0.0;

            // x-축선위에서 축선의 오른쪽점에서의 u 계수인 점근적함수의 적분 계산
            xdat->Eu += (1.0/(yb*yb*yb*yb)*(xb-xp)*((xp*xp)*(yb*yb)+(xp*xp*xp)*sqrt(xp*xp+yb*yb)*2.0+(xp*xp*xp*xp)*2.0-yb*yb*yb*yb)*-2.0)/xp;

            // x-축선위에서  x < xi인 경우의 그린함수와 phi1의 적분의 계산
            xdat->F += 0.0;

            // x-축선위에서  x < xi인 경우의 그린함수와 phi2의 적분의 계산 (f의 계수)
            xdat->F += 0.0;

            // x-축선위에서  x < xi인 경우의 그린함수와 phi2의 적분의 계산 (Eu의 계수)
            xdat->Eu += ((xp*xp)*1.0/(yb*yb*yb*yb)*1.0/pow(xb*xb+yb*yb,3.0/2.0)*sqrt(xp*xp+yb*yb)*(xb-xp)*((xb*xb)*(yb*yb)*6.0+xb*pow(xb*xb+yb*yb,3.0/2.0)*4.0+(xb*xb*xb*xb)*4.0+yb*yb*yb*yb))/xb;
          }
        }
      }
      // 축대칭문제인 경우
      if (nD == 3) {
        // 점금적함수가 u (x, y) ~ c / sqrt (z^2 + r^2) 인 경우
        if (AsymptoticType == 1) {
          // x-축선위에서 점근적함수의 적분 계산
          xdat->F += 0.0;
          // x-축선위에서 축선의 오른쪽점에서의 u 계수인 점근적함수의 적분 계산
          xdat->Eu += (1.0/(yb*yb*yb*yb)*(xb-xp)*((xp*xp)*(yb*yb)+(xp*xp*xp)*sqrt(xp*xp+yb*yb)*2.0+(xp*xp*xp*xp)*2.0-yb*yb*yb*yb)*-2.0)/xp;
          // x-축선위에서 축선의 오른쪽점에서의 u 계수인 점근적함수의 적분 계산 (phi의 적분)
          xdat->Eu += ((xp*xp)*1.0/(yb*yb*yb*yb)*1.0/pow(xb*xb+yb*yb,3.0/2.0)*sqrt(xp*xp+yb*yb)*(xb-xp)*((xb*xb)*(yb*yb)*6.0+xb*pow(xb*xb+yb*yb,3.0/2.0)*4.0+(xb*xb*xb*xb)*4.0+yb*yb*yb*yb))/xb;
          // x-축선위에서  x < xi인 경우의 그린함수와 phi의 적분의 계산
          // xdat->F += ((xp*xp)*1.0/(yb*yb*yb*yb)*1.0/pow(xb*xb+yb*yb,3.0/2.0)*(xb-xp)*((xb*xb)*(yb*yb)*6.0+xb*pow(xb*xb+yb*yb,3.0/2.0)*4.0+(xb*xb*xb*xb)*4.0+yb*yb*yb*yb))/xb;
          // x-축선위에서  x > xi인 경우의 그린함수와 phi의 적분의 계산
          xdat->F += 0.0;
        }
      }
    }
    // x-축선 위에서 Dirichlet - Infinity 그린함수를 사용하는 경우
    if (bdx == 4) {
      // 2차원문제인 경우
      if (nD == 2) {
        // 점근적함수가 log (sqrt (x^2 + y^2)) + c * 1인 경우
        if (AsymptoticType == 1) {
          // yb가 0인 경우
          if (IsEqualDouble (yb, ZeroValue)) {
            // x-축선위에서 점근적함수의 적분 계산
            xdat->F += (xb*(2.0/3.0))/xm-2.0/3.0;

            // x-축선위에서 축선의 왼쪽점에서의 u 계수인 점근적함수의 적분 계산
            xdat->Wu += (xb*2.0)/xm-2.0;

            // x-축선위에서  x > xi인 경우의 그린함수와 phi1의 적분의 계산
            xdat->F += 1.0/(xb*xb*xb)*(xm*xm)*(xb-xm)*(1.0/3.0);

            // x-축선위에서  x > xi인 경우의 그린함수와 phi2의 적분의 계산 (f의 계수)
            xdat->F += 0.0;

            // x-축선위에서  x > xi인 경우의 그린함수와 phi2의 적분의 계산
            xdat->Wu += 0.0;
          }
          // yb가 0이 아닌 경우
          else {
            // x-축선위에서 점근적함수의 적분 계산
            xdat->F += xm*1.0/(yb*yb*yb)*(xb-xm)*(yb*2.0-xm*3.141592653589793+xm*atan(xm/yb)*2.0);

            // x-축선위에서 축선의 왼쪽점에서의 u 계수인 점근적함수의 적분 계산
            xdat->Wu += (xb*2.0)/xm-2.0;

            // x-축선위에서  x > xi인 경우의 그린함수와 phi1의 적분의 계산
            xdat->F += -((xm*xm)*1.0/(yb*yb*yb)*(xb-xm)*((xb*xb*xb)*atan(xb/yb)*2.0-(xb*xb*xb)*3.141592653589793+(xb*xb)*yb*2.0+yb*yb*yb-xb*(yb*yb)*3.141592653589793+xb*(yb*yb)*atan(xb/yb)*2.0))/(xb*(xb*xb+yb*yb));

            // x-축선위에서  x > xi인 경우의 그린함수와 phi2의 적분의 계산 (f의 계수)
            xdat->F += 0.0;

            // x-축선위에서  x > xi인 경우의 그린함수와 phi2의 적분의 계산
            xdat->Wu += 0.0;
          }
        }

        // 점근적함수가 c * 1.0E0 / sqrt (x ^ 2 + y ^ 2)인 경우
        if (AsymptoticType == 2) {
          // yb가 0인 경우
          if (IsEqualDouble (yb, ZeroValue)) {
            // x-축선위에서 점근적함수의 적분 계산
            xdat->F += 0.0;

            // x-축선위에서 축선의 왼쪽점에서의 u 계수인 점근적함수의 적분 계산
            xdat->Wu += (xb*(3.0/2.0))/xm-3.0/2.0;

            // x-축선위에서  x > xi인 경우의 그린함수와 phi1의 적분의 계산
            xdat->F += 0.0;

            // x-축선위에서  x > xi인 경우의 그린함수와 phi2의 적분의 계산 (f의 계수)
            xdat->F += 0.0;

            // x-축선위에서  x > xi인 경우의 그린함수와 phi2의 적분의 계산
            xdat->Wu += 1.0/(xb*xb*xb*xb)*(xm*xm*xm)*(xb-xm)*(-1.0/2.0);
          }
          // yb가 0이 아닌 경우
          else {
            // x-축선위에서 점근적함수의 적분 계산
            xdat->F += 0.0;

            // x-축선위에서 축선의 왼쪽점에서의 u 계수인 점근적함수의 적분 계산
            xdat->Wu += (1.0/(yb*yb*yb*yb)*(xb-xm)*((xm*xm)*(yb*yb)-(xm*xm*xm)*sqrt(xm*xm+yb*yb)*2.0+(xm*xm*xm*xm)*2.0-yb*yb*yb*yb)*-2.0)/xm;

            // x-축선위에서  x > xi인 경우의 그린함수와 phi1의 적분의 계산
            xdat->F += 0.0;

            // x-축선위에서  x > xi인 경우의 그린함수와 phi2의 적분의 계산 (f의 계수)
            xdat->F += 0.0;

            // x-축선위에서  x > xi인 경우의 그린함수와 phi2의 적분의 계산
            xdat->Wu += ((xm*xm)*1.0/(yb*yb*yb*yb)*1.0/pow(xb*xb+yb*yb,3.0/2.0)*sqrt(xm*xm+yb*yb)*(xb-xm)*((xb*xb)*(yb*yb)*6.0-xb*pow(xb*xb+yb*yb,3.0/2.0)*4.0+(xb*xb*xb*xb)*4.0+yb*yb*yb*yb))/xb;
          }
        }
      }
      // 축대칭문제인 경우
      if (nD == 3) {
        // 점금적함수가 u (x, y) ~ c / sqrt (z^2 + r^2) 인 경우
        if (AsymptoticType == 1) {
          // x-축선위에서 점근적함수의 적분 계산
          xdat->F += 0.0;
          // x-축선위에서 축선의 왼쪽점에서의 u 계수인 점근적함수의 적분 계산
          xdat->Wu += (1.0/(yb*yb*yb*yb)*(xb-xm)*((xm*xm)*(yb*yb)-(xm*xm*xm)*sqrt(xm*xm+yb*yb)*2.0+(xm*xm*xm*xm)*2.0-yb*yb*yb*yb)*-2.0)/xm;
          // x-축선위에서 축선의 왼쪽점에서의 u 계수인 점근적함수의 적분 계산 (phi의 적분)
          xdat->Wu += ((xm*xm)*1.0/(yb*yb*yb*yb)*1.0/pow(xb*xb+yb*yb,3.0/2.0)*sqrt(xm*xm+yb*yb)*(xb-xm)*((xb*xb)*(yb*yb)*6.0-xb*pow(xb*xb+yb*yb,3.0/2.0)*4.0+(xb*xb*xb*xb)*4.0+yb*yb*yb*yb))/xb;
          // x-축선위에서  x < xi인 경우의 그린함수와 phi의 적분의 계산
          xdat->F += 0.0;
          // x-축선위에서  x > xi인 경우의 그린함수와 phi의 적분의 계산
          // xdat->F += ((xm*xm)*1.0/(yb*yb*yb*yb)*1.0/pow(xb*xb+yb*yb,3.0/2.0)*(xb-xm)*((xb*xb)*(yb*yb)*6.0-xb*pow(xb*xb+yb*yb,3.0/2.0)*4.0+(xb*xb*xb*xb)*4.0+yb*yb*yb*yb))/xb;
        }
      }
    }
  }
  // y-축선을 참조하는 경우
  if (xy == 'Y' || xy == 'y') {
    // y-축선 위에서 Infinity - Dirichlet 그린함수를 사용하는 경우
    if (bdy == 3) {
      // 2차원문제인 경우
      if (nD == 2) {
        // 점근적함수가 log (sqrt (x^2 + y^2)) + c * 1인 경우
        if (AsymptoticType == 1) {
          // xb가 0인 경우
          if (IsEqualDouble (xb, ZeroValue)) {
            // y-축선위에서 점근적 함수의 적분 계산
            ydat->F += (yb*(2.0/3.0))/yp-2.0/3.0;

            // y-축선위에서 축선의 위쪽점에서의 u 계수인 점근적함수의 적분 계산
            ydat->Nu += (yb*2.0)/yp-2.0;

            // y-축선위에서 y < eta인 경우의 그린함수와 f-phi1의 적분의 계산
            ydat->F += 1.0/(yb*yb*yb)*(yp*yp)*(yb-yp)*(1.0/3.0);

            // y-축선위에서 y < eta인 경우의 그린함수와 -phi2의 적분의 계산 (f의 계수)
            ydat->F += 0.0;

            // y-축선위에서 y < eta인 경우의 그린함수와 -phi2의 적분의 계산 (Nu의 계수)
            ydat->Nu += 0.0;
          }
          // xb가 0이 아닌 경우
          else {
            // y-축선위에서 점근적 함수의 적분 계산
            ydat->F += 1.0/(xb*xb*xb)*yp*(xb*2.0-yp*atan(xb/yp)*2.0)*(yb-yp);

            // y-축선위에서 축선의 위쪽점에서의 u 계수인 점근적함수의 적분 계산
            ydat->Nu += (yb*2.0)/yp-2.0;

            // y-축선위에서 y < eta인 경우의 그린함수와 f-phi1의 적분의 계산
            ydat->F += (1.0/(xb*xb*xb)*(yp*yp)*(yb-yp)*((yb*yb*yb)*atan(xb/yb)*2.0-xb*(yb*yb)*2.0-xb*xb*xb+(xb*xb)*yb*atan(xb/yb)*2.0))/(yb*(xb*xb+yb*yb));

            // y-축선위에서 y < eta인 경우의 그린함수와 -phi2의 적분의 계산 (f의 계수)
            ydat->F += 0.0;

            // y-축선위에서 y < eta인 경우의 그린함수와 -phi2의 적분의 계산 (Nu의 계수)
            ydat->Nu += 0.0;
          }
        }

        // 점근적함수가 c * 1.0E0 / sqrt (x ^ 2 + y ^ 2)인 경우
        if (AsymptoticType == 2) {
          // xb가 0인 경우
          if (IsEqualDouble (xb, ZeroValue)) {
            // y-축선위에서 점근적 함수의 적분 계산
            ydat->F += 0.0;

            // y-축선위에서 축선의 위쪽점에서의 u 계수인 점근적함수의 적분 계산
            ydat->Nu += (yb*(3.0/2.0)-yp*(3.0/2.0))/yp;

            // y-축선위에서 y < eta인 경우의 그린함수와 f-phi1의 적분의 계산
            ydat->F += 1.0/(yb*yb*yb*yb)*(yp*yp)*(yb-yp)*(1.0/4.0);

            // y-축선위에서 y < eta인 경우의 그린함수와 -phi2의 적분의 계산 (f의 계수)
            ydat->F += 0.0;

            // y-축선위에서 y < eta인 경우의 그린함수와 -phi2의 적분의 계산 (Nu의 계수)
            ydat->Nu += 1.0/(yb*yb*yb*yb)*(yp*yp*yp)*(yb-yp)*(-1.0/4.0);
          }
          // xb가 0이 아닌 경우
          else {
            // y-축선위에서 점근적 함수의 적분 계산
            ydat->F += 0.0;

            // y-축선위에서 축선의 위쪽점에서의 u 계수인 점근적함수의 적분 계산
            ydat->Nu += (1.0/(xb*xb*xb*xb)*(yb-yp)*((xb*xb)*(yp*yp)+(yp*yp*yp)*sqrt(xb*xb+yp*yp)*2.0-xb*xb*xb*xb+(yp*yp*yp*yp)*2.0)*-2.0)/yp;

            // y-축선위에서 y < eta인 경우의 그린함수와 f-phi1의 적분의 계산
            ydat->F += (1.0/(xb*xb*xb*xb)*(yp*yp)*1.0/sqrt(xb*xb+yb*yb)*(yb-yp)*(yb*sqrt(xb*xb+yb*yb)*2.0+xb*xb+(yb*yb)*2.0))/yb;

            // y-축선위에서 y < eta인 경우의 그린함수와 -phi2의 적분의 계산 (f의 계수)
            ydat->F += 0.0;

            // y-축선위에서 y < eta인 경우의 그린함수와 -phi2의 적분의 계산 (Nu의 계수)
            ydat->Nu += (1.0/(xb*xb*xb*xb)*(yp*yp)*1.0/pow(xb*xb+yb*yb,3.0/2.0)*sqrt(xb*xb+yp*yp)*(yb-yp)*((xb*xb)*(yb*yb)*9.0+yb*pow(xb*xb+yb*yb,3.0/2.0)*6.0+(xb*xb*xb*xb)*2.0+(yb*yb*yb*yb)*6.0))/yb;
          }
        }
      }
    }
    // y-축선 위에서 Dirichlet - Infinity 그린함수를 사용하는 경우
    if (bdy == 4) {
      // 2차원문제인 경우
      if (nD == 2) {
        // 점근적함수가 log (sqrt (x^2 + y^2)) + c * 1인 경우
        if (AsymptoticType == 1) {
          // xb가 0인 경우
          if (IsEqualDouble (xb, ZeroValue)) {
            // y-축선위에서 점근적함수의 적분 계산
            ydat->F += (yb*(2.0/3.0))/ym-2.0/3.0;

            // y-축선위에서 축선의 아래쪽점에서의 u 계수인 점근적함수의 적분 계산
            ydat->Su += (yb*2.0)/ym-2.0;

            // y-축선위에서  y > eta인 경우의 그린함수와 f-phi1의 적분의 계산
            ydat->F += 1.0/(yb*yb*yb)*(ym*ym)*(yb-ym)*(1.0/3.0);

            // y-축선위에서  y > eta인 경우의 그린함수와 -phi2의 적분의 계산 (f의 계수)
            ydat->F += 0.0;

            // y-축선위에서  y > eta인 경우의 그린함수와 -phi2의 적분의 계산 (Su의 계수)
            ydat->Su += 0.0;
          }
          // xb가 0이 아닌 경우
          else {
            // y-축선위에서 점근적함수의 적분 계산
            ydat->F += 1.0/(xb*xb*xb)*ym*(xb*2.0-ym*atan(xb/ym)*2.0)*(yb-ym);

            // y-축선위에서 축선의 아래쪽점에서의 u 계수인 점근적함수의 적분 계산
            ydat->Su += (yb*2.0)/ym-2.0;

            // y-축선위에서  y > eta인 경우의 그린함수와 f-phi1의 적분의 계산
            ydat->F += (1.0/(xb*xb*xb)*(ym*ym)*(yb-ym)*((yb*yb*yb)*atan(xb/yb)*2.0-xb*(yb*yb)*2.0-xb*xb*xb+(xb*xb)*yb*atan(xb/yb)*2.0))/(yb*(xb*xb+yb*yb));

            // y-축선위에서  y > eta인 경우의 그린함수와 -phi2의 적분의 계산 (f의 계수)
            ydat->F += 0.0;

            // y-축선위에서  y > eta인 경우의 그린함수와 -phi2의 적분의 계산 (Su의 계수)
            ydat->Su += 0.0;
          }
        }

        // 점근적함수가 c * 1.0E0 / sqrt (x ^ 2 + y ^ 2)인 경우
        if (AsymptoticType == 2) {
          // xb가 0인 경우
          if (IsEqualDouble (xb, ZeroValue)) {
            // y-축선위에서 점근적함수의 적분 계산
            ydat->F += 0.0;

            // y-축선위에서 축선의 아래쪽점에서의 u 계수인 점근적함수의 적분 계산
            ydat->Su += (yb*(3.0/2.0))/ym-3.0/2.0;

            // y-축선위에서  y > eta인 경우의 그린함수와 f-phi1의 적분의 계산
            ydat->F += 1.0/(yb*yb*yb*yb)*(ym*ym)*(yb-ym)*(-1.0/4.0);

            // y-축선위에서  y > eta인 경우의 그린함수와 -phi2의 적분의 계산 (f의 계수)
            ydat->F += 0.0;

            // y-축선위에서  y > eta인 경우의 그린함수와 -phi2의 적분의 계산 (Su의 계수)
            ydat->Su += 1.0/(yb*yb*yb*yb)*(ym*ym*ym)*(yb-ym)*(-1.0/4.0);
          }
          // xb가 0이 아닌 경우
          else {
            // y-축선위에서 점근적함수의 적분 계산
            ydat->F += 0.0;

            // y-축선위에서 축선의 아래쪽점에서의 u 계수인 점근적함수의 적분 계산
            ydat->Su += (1.0/(xb*xb*xb*xb)*(yb-ym)*((xb*xb)*(ym*ym)-(ym*ym*ym)*sqrt(xb*xb+ym*ym)*2.0-xb*xb*xb*xb+(ym*ym*ym*ym)*2.0)*-2.0)/ym;

            // y-축선위에서  y > eta인 경우의 그린함수와 f-phi1의 적분의 계산
            ydat->F += -(1.0/(xb*xb*xb*xb)*(ym*ym)*1.0/sqrt(xb*xb+yb*yb)*(yb-ym)*(yb*sqrt(xb*xb+yb*yb)*-2.0+xb*xb+(yb*yb)*2.0))/yb;

            // y-축선위에서  y > eta인 경우의 그린함수와 -phi2의 적분의 계산 (f의 계수)
            ydat->F += 0.0;

            // y-축선위에서  y > eta인 경우의 그린함수와 -phi2의 적분의 계산 (Su의 계수)
            ydat->Su += (1.0/(xb*xb*xb*xb)*(ym*ym)*1.0/pow(xb*xb+yb*yb,3.0/2.0)*sqrt(xb*xb+ym*ym)*(yb-ym)*((xb*xb)*(yb*yb)*9.0-yb*pow(xb*xb+yb*yb,3.0/2.0)*6.0+(xb*xb*xb*xb)*2.0+(yb*yb*yb*yb)*6.0))/yb;
          }
        }
      }
      // 축대칭문제인 경우
      if (nD == 3) {
        if (AsymptoticType == 1) {
          // y-축선의 x-좌표가 0인 경우
          if (IsEqualDouble (xb, 0.0)) {
            // y-축선위에서 점근적함수의 적분 계산
            ydat->F += 0.0;
            // y-축선위에서 축선의 아래쪽점에서의 u 계수인 점근적함수의 적분 계산
            ydat->Su += log(fabs(yb))*(9.0/4.0)-log(fabs(ym))*(9.0/4.0);
            // y-축선위에서 축선의 아래쪽점에서의 u 계수인 점근적함수의 적분 계산 (phi의 적분)
            ydat->Su += 1.0/(yb*yb*yb*yb)*(ym*ym*ym*ym)*(log(fabs(yb))-log(fabs(ym)))*(-1.0/4.0);
            // y-축선위에서  y < eta인 경우의 그린함수와 -phi의 적분의 계산
            // ydat->F += 0.0;
            // y-축선위에서  y > eta인 경우의 그린함수와 -phi의 적분의 계산
            // ydat->F += 1.0/(yb*yb*yb*yb)*(ym*ym*ym)*(log(fabs(yb))-log(fabs(ym)))*(-1.0/4.0);
          }
          // y-축선의 x-좌표가 0이 아닌 경우
          else {
            // y-축선위에서 점근적함수의 적분 계산
            ydat->F += 0.0;
            // y-축선위에서 축선의 아래쪽점에서의 u 계수인 점근적함수의 적분 계산
            ydat->Su += 1.0/(xb*xb*xb*xb)*(log(fabs(yb))-log(fabs(ym)))*((xb*xb)*(ym*ym)-(ym*ym*ym)*sqrt(xb*xb+ym*ym)*2.0-xb*xb*xb*xb+(ym*ym*ym*ym)*2.0)*-3.0;
            // y-축선위에서 축선의 아래쪽점에서의 u 계수인 점근적함수의 적분 계산 (phi의 적분)
            ydat->Su += (1.0/(xb*xb*xb*xb)*(ym*ym*ym)*1.0/pow(xb*xb+yb*yb,3.0/2.0)*sqrt(xb*xb+ym*ym)*(log(fabs(yb))-log(fabs(ym)))*((xb*xb)*(yb*yb)*9.0-yb*pow(xb*xb+yb*yb,3.0/2.0)*6.0+(xb*xb*xb*xb)*2.0+(yb*yb*yb*yb)*6.0))/yb;
            // y-축선위에서  y < eta인 경우의 그린함수와 -phi의 적분의 계산
            ydat->F += 0.0;
            // y-축선위에서  y > eta인 경우의 그린함수와 -phi의 적분의 계산
            // ydat->F += (1.0/(xb*xb*xb*xb)*(ym*ym*ym)*1.0/pow(xb*xb+yb*yb,3.0/2.0)*(log(fabs(yb))-log(fabs(ym)))*((xb*xb)*(yb*yb)*9.0-yb*pow(xb*xb+yb*yb,3.0/2.0)*6.0+(xb*xb*xb*xb)*2.0+(yb*yb*yb*yb)*6.0))/yb;
          }
        }
      }
    }
  }
  // 잘못된 참조를 한 경우, 에러메시지를 출력하고 종료
  if (xy != 'X' && xy != 'x' && xy != 'Y' && xy != 'y') PrintError ("ExactIntegration");
}

// 해의 표현식의 계수의 확인하는 모듈
double CheckRepresentationFormula (Point *pt, xData *xdat, yData *ydat, int idx, double xm, double xb, double xp, double ym, double yb, double yp, int bdx, int bdy) {
  double xx = xdat->Cu * u_ftn (xb, yb) + xdat->Eu * dudx_ftn (xp, yb) + xdat->Wu * u_ftn (xm, yb) + xdat->Cphi * phi_ftn (xb, yb) + xdat->Ephi * phi_ftn (xp, yb) + xdat->Wphi * phi_ftn (xm, yb) - xdat->F;
  double yy = ydat->Cu * u_ftn (xb, yb) + ydat->Nu * u_ftn (xb, yp) + ydat->Su * u_ftn (xb, ym) + ydat->Cphi * phi_ftn (xb, yb) + ydat->Nphi * phi_ftn (xb, yp) + ydat->Sphi * phi_ftn (xb, ym) - ydat->F;
  double zz = 0.0;

  if (xx > yy) zz = fabs(xx);
  else         zz = fabs(yy);

  // if (zz > 1e-10) {
  if (pt->Index () == idx) {
    printf ("\n");
    printf ("%-8s%23.16e\n", "XM   = ", xm);
    printf ("%-8s%23.16e\n", "XB   = ", xb);
    printf ("%-8s%23.16e\n", "XP   = ", xp);
    printf ("%-8s%23.16e\n", "YM   = ", ym);
    printf ("%-8s%23.16e\n", "YB   = ", yb);
    printf ("%-8s%23.16e\n", "YP   = ", yp);
    printf ("\n");
    printf ("%-8s%c\n"     , "cond = ", pt->Condition ());
    printf ("%-8s%d\n"     , "bdx  = ", bdx);
    printf ("%-8s%d\n"     , "bdy  = ", bdy);
    printf ("\n");
    printf ("%-8s%23.16e\n", "Cu   = ", xdat->Cu  );
    printf ("%-8s%23.16e\n", "Eu   = ", xdat->Eu  );
    printf ("%-8s%23.16e\n", "Wu   = ", xdat->Wu  );
    printf ("%-8s%23.16e\n", "Cphi = ", xdat->Cphi);
    printf ("%-8s%23.16e\n", "Ephi = ", xdat->Ephi);
    printf ("%-8s%23.16e\n", "Wphi = ", xdat->Wphi);
    printf ("%-8s%23.16e\n", "F    = ", xdat->F   );
    printf ("%-8s%23.16e\n", "EF   = ", xdat->Ef   );
    printf ("%-8s%23.16e\n", "CF   = ", xdat->Cf   );
    printf ("%-8s%23.16e\n", "WF   = ", xdat->Wf   );
    printf ("\n");
    printf ("%-8s%23.16e\n", "Cu   = ", ydat->Cu  );
    printf ("%-8s%23.16e\n", "Nu   = ", ydat->Nu  );
    printf ("%-8s%23.16e\n", "Su   = ", ydat->Su  );
    printf ("%-8s%23.16e\n", "Cphi = ", ydat->Cphi);
    printf ("%-8s%23.16e\n", "Nphi = ", ydat->Nphi);
    printf ("%-8s%23.16e\n", "Sphi = ", ydat->Sphi);
    printf ("%-8s%23.16e\n", "F    = ", ydat->F   );
    printf ("%-8s%23.16e\n", "Nf    = ", ydat->Nf   );
    printf ("%-8s%23.16e\n", "Cf    = ", ydat->Cf   );
    printf ("%-8s%23.16e\n", "Sf    = ", ydat->Sf   );
    printf ("\n");
    printf ("%-8s%23.16e\n", "u (xb, yb) = ", u_ftn (xb, yb));
    printf ("%-8s%23.16e\n", "u (xm, yb) = ", u_ftn (xm, yb));
    printf ("%-8s%23.16e\n", "u (xp, yb) = ", u_ftn (xp, yb));
    printf ("%-8s%23.16e\n", "u (xb, ym) = ", u_ftn (xb, ym));
    printf ("%-8s%23.16e\n", "u (xb, yp) = ", u_ftn (xb, yp));
    printf ("%-8s%23.16e\n", "phi (xb, yb) = ", phi_ftn (xb, yb));
    printf ("%-8s%23.16e\n", "phi (xm, yb) = ", phi_ftn (xm, yb));
    printf ("%-8s%23.16e\n", "phi (xp, yb) = ", phi_ftn (xp, yb));
    printf ("%-8s%23.16e\n", "phi (xb, ym) = ", phi_ftn (xb, ym));
    printf ("%-8s%23.16e\n", "phi (xb, yp) = ", phi_ftn (xb, yp));
    printf ("\n");
    printf ("%-35s%23.16e\n", "Representation formula on x = ", (xdat->Cu * u_ftn (xb, yb) + xdat->Eu * dudx_ftn (xp, yb) + xdat->Wu * u_ftn (xm, yb) + xdat->Cphi * phi_ftn (xb, yb) + xdat->Ephi * phi_ftn (xp, yb) + xdat->Wphi * phi_ftn (xm, yb) - xdat->F));
    printf ("%-35s%23.16e\n", "Representation formula on y = ", (ydat->Cu * u_ftn (xb, yb) + ydat->Nu * u_ftn (xb, yp) + ydat->Su * u_ftn (xb, ym) + ydat->Cphi * phi_ftn (xb, yb) + ydat->Nphi * phi_ftn (xb, yp) + ydat->Sphi * phi_ftn (xb, ym) - ydat->F));

    return zz;
  }
  // }
  exit(1);
}

/* y-축선에서의 u의 계수를 계산 */
double CalcVerticalUCoefficient (char NS, double xb, double yb, double yp, double ym, int bdy, double mp) {
  // 위쪽점을 참조하는 경우, y-축선의 위쪽점에서의 u의 계수를 return
  if (NS == 'N') return greens_coefficient_t (yp, ym, yb, yp, xb, yb, 2, bdy, mp, mp);
  // 아래쪽점을 참조하는 경우, y-축선의 아래쪽점에서의 u의 계수를 return
  if (NS == 'S') return greens_coefficient_t (ym, ym, yb, yp, xb, yb, 2, bdy, mp, mp);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcVerticalUCoefficient");
  exit (1);
}

/* x-축선에서의 u의 계수를 계산 */
double CalcHorizontalUCoefficient (char EW, double xb, double yb, double xp, double xm, int bdx, double mp) {
  // 오른쪽점을 참조하는 경우, x-축선의 오른쪽점에서의 u의 계수를 return
  if (EW == 'E') return greens_coefficient_t (xp, xm, xb, xp, xb, yb, 1, bdx, mp, mp);
  // 왼쪽점을 참조하는 경우, x-축선의 왼쪽점에서의 u의 계수를 return
  if (EW == 'W') return greens_coefficient_t (xm, xm, xb, xp, xb, yb, 1, bdx, mp, mp);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcHorizontalUCoefficient");
  exit (1);
}

/* y-축선에서의 phi의 계수를 계산 */
double CalcVerticalPHICoefficient (char NS, double xb, double yb, double yp, double ym, int bdy, double mp) {
  // 위쪽점을 참조하는 경우, y-축선의 위쪽점에서의 phi의 계수를 return
  if (NS == 'N') return -greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mp, mp) - (greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) * (yb - ym) / (yp - ym);
  // 아래쪽점을 참조하는 경우, y-축선의 아래쪽점에서의 phi의 계수를 return
  if (NS == 'S') return -greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mp, mp) - (greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) * (yp - yb) / (yp - ym);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcVerticalPHICoefficient");
  exit (1);
}

/* x-축선에서의 phi의 계수를 계산 */
double CalcHorizontalPHICoefficient (char EW, double xb, double yb, double xp, double xm, int bdx, double mp) {
  // 오른쪽점을 참조하는 경우, x-축선의 오른쪽점에서의 phi의 계수를 return
  if (EW == 'E') return greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + (greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) * (xb - xm) / (xp - xm);
  // 왼쪽점을 참조하는 경우, x-축선의 왼쪽점에서의 phi의 계수를 return
  if (EW == 'W') return greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + (greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) * (xp - xb) / (xp - xm);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcHorizontalPHICoefficient");
  exit (1);
}

/* y-축선에서의 f의 값을 계산 */
double CalcVerticalFCoefficient (char NS, double xb, double yb, double yp, double ym, int bdy, double mp) {
  // 위쪽점을 참조하는 경우, y-축선의 위쪽점에서의 f의 계수를 return
  if (NS == 'N') return greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + (greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) * (yb - ym) / (yp - ym);
  // 아래쪽점을 참조하는 경우, y-축선의 아래쪽점에서의 f의 계수를 return
  if (NS == 'S') return greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + (greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) * (yp - yb) / (yp - ym);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcVerticalFCoefficient");
  exit (1);
}

/* x-축선에서의 f의 값을 계산 */
double CalcHorizontalFCoefficient (char EW, double xb, double yb, double xp, double xm, int bdx, double mp) {
  // 오른쪽점을 참조하는 경우, x-축선의 오른쪽점에서의 phi의 계수를 return
  if (EW == 'E') return greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + (greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) * (xb - xm) / (xp - xm);
  // 왼쪽점을 참조하는 경우, x-축선의 왼쪽점에서의 phi의 계수를 return
  if (EW == 'W') return greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + (greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) * (xp - xb) / (xp - xm);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcHorizontalFCoefficient");
  exit (1);
}

/* y-축선에서의 u의 계수를 계산 */
double CalcVerticalDiffUCoefficient (char NS, double xb, double yb, double yp, double ym, int bdy, double mp) {
  // 위쪽점을 참조하는 경우, y-축선의 위쪽점에서의 u의 계수를 return
  if (NS == 'N') return greens_coefficient_ttau (yp, ym, yb, yp, xb, yb, 2, bdy, mp, mp);
  // 아래쪽점을 참조하는 경우, y-축선의 아래쪽점에서의 u의 계수를 return
  if (NS == 'S') return greens_coefficient_ttau (ym, ym, yb, yp, xb, yb, 2, bdy, mp, mp);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcVerticalDiffUCoefficient");
  exit (1);
}

/* x-축선에서의 u의 계수를 계산 */
double CalcHorizontalDiffUCoefficient (char EW, double xb, double yb, double xp, double xm, int bdx, double mp) {
  // 오른쪽점을 참조하는 경우, x-축선의 오른쪽점에서의 u의 계수를 return
  if (EW == 'E') return greens_coefficient_ttau (xp, xm, xb, xp, xb, yb, 1, bdx, mp, mp);
  // 왼쪽점을 참조하는 경우, x-축선의 왼쪽점에서의 u의 계수를 return
  if (EW == 'W') return greens_coefficient_ttau (xm, xm, xb, xp, xb, yb, 1, bdx, mp, mp);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcHorizontalDiffUCoefficient");
  exit (1);
}

/* y-축선에서의 phi의 계수를 계산 */
double CalcVerticalDiffPHICoefficient (char NS, double xb, double yb, double yp, double ym, int bdy, double mp) {
  // 위쪽점을 참조하는 경우, y-축선의 위쪽점에서의 phi의 계수를 return
  if (NS == 'N') return - greens_integral_tau (4, ym, yb, yp, xb, yb, 2, bdy, mp, mp) - (greens_integral_tau (2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral_tau (3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) * (yb - ym) / (yp - ym);
  // 아래쪽점을 참조하는 경우, y-축선의 아래쪽점에서의 phi의 계수를 return
  if (NS == 'S') return - greens_integral_tau (1, ym, yb, yp, xb, yb, 2, bdy, mp, mp) - (greens_integral_tau (2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral_tau (3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) * (yp - yb) / (yp - ym);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcVerticalDiffPHICoefficient");
  exit (1);
}

/* x-축선에서의 phi의 계수를 계산 */
double CalcHorizontalDiffPHICoefficient (char EW, double xb, double yb, double xp, double xm, int bdx, double mp) {
  // 오른쪽점을 참조하는 경우, x-축선의 오른쪽점에서의 phi의 계수를 return
  if (EW == 'E') return greens_integral_tau (4, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + (greens_integral_tau (2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral_tau (3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) * (xb - xm) / (xp - xm);
  // 왼쪽점을 참조하는 경우, x-축선의 왼쪽점에서의 phi의 계수를 return
  if (EW == 'W') return greens_integral_tau (1, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + (greens_integral_tau (2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral_tau (3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) * (xp - xb) / (xp - xm);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcHorizontalDiffPHICoefficient");
  exit (1);
}

/* y-축선에서의 f의 값을 계산 */
double CalcVerticalDiffFCoefficient (char NS, double xb, double yb, double yp, double ym, int bdy, double mp) {
  // 위쪽점을 참조하는 경우, y-축선의 위쪽점에서의 f의 계수를 return
  if (NS == 'N') return greens_integral_tau (4, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + (greens_integral_tau (2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral_tau (3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) * (yb - ym) / (yp - ym);
  // 아래쪽점을 참조하는 경우, y-축선의 아래쪽점에서의 f의 계수를 return
  if (NS == 'S') return greens_integral_tau (1, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + (greens_integral_tau (2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral_tau (3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) * (yp - yb) / (yp - ym);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcVerticalFDiffCoefficient");
  exit (1);
}

/* x-축선에서의 f의 값을 계산 */
double CalcHorizontalDiffFCoefficient (char EW, double xb, double yb, double xp, double xm, int bdx, double mp) {
  // 오른쪽점을 참조하는 경우, x-축선의 오른쪽점에서의 phi의 계수를 return
  if (EW == 'E') return greens_integral_tau (4, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + (greens_integral_tau (2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral_tau (3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) * (xb - xm) / (xp - xm);
  // 왼쪽점을 참조하는 경우, x-축선의 왼쪽점에서의 phi의 계수를 return
  if (EW == 'W') return greens_integral_tau (1, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + (greens_integral_tau (2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral_tau (3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) * (xp - xb) / (xp - xm);
  // 잘못된 참조의 경우, 에러메시지를 출력하고 종료
  PrintError ("CalcHorizontalFDiffCoefficient");
  exit (1);
}

double PhiInterp (Point *pt) {
  // 국소 x-축선의 시작점과 끝점의 x-좌표와 현재점의 x-좌표
  double xm = pt->MinMaxCoordinate ('x', 'm'), xb = pt->Coord ().Value ('x'), xp = pt->MinMaxCoordinate ('x', 'p');
  // 국소 y-축선의 시작점과 끝점의 y-좌표와 현재점의 y-좌표
  double ym = pt->MinMaxCoordinate ('y', 'm'), yb = pt->Coord ().Value ('y'), yp = pt->MinMaxCoordinate ('y', 'p');
  // 오른쪽 phi의 값을 저장하기 위한 변수
  double phiE = ZeroValue;
  // 왼쪽 phi의 값을 저장하기 위한 변수
  double phiW = ZeroValue;
  // 위쪽 phi의 값을 저장하기 위한 변수
  double phiN = ZeroValue;
  // 아래쪽 phi의 값을 저장하기 위한 변수
  double phiS = ZeroValue;
  // x-축선 위에서 phi의 값을 근사하기 위한 변수
  double phiX = ZeroValue;
  // y-축선 위에서 phi의 값을 근사하기 위한 변수
  double phiY = ZeroValue;

  // 오른쪽점이 존재하는 경우
  if (pt->EWNS ('E', 'E') != NULL) {
    phiE = pt->EWNS ('E', 'E')->Phi ()->Value ();
  }
  // 오른쪽 위의 점과 오른쪽 아래의 점이 존재하는 경우
  else if (pt->EWNS ('E', 'N') != NULL && pt->EWNS ('E', 'S') != NULL) {
    // 오른쪽 위의 점의 y-좌표를 저장
    yp = pt->EWNS ('E', 'N')->Coord ().Value ('y');
    // 오른쪽 아래의 점의 y-좌표를 저장
    ym = pt->EWNS ('E', 'S')->Coord ().Value ('y');
    // 오른쪽점에서의 phi의 값을 근사
    phiE = pt->EWNS ('E', 'N')->Phi ()->Value () * (yb - ym) / (yp - ym) + pt->EWNS ('E', 'S')->Phi ()->Value () * (yp - yb) / (yp - ym);
  }
  // 둘 다 아닌 경우 에러메시지를 출력
  else PrintError ("PhiInterp");

  // 왼쪽점이 존재하는 경우
  if (pt->EWNS ('W', 'W') != NULL) {
    phiW = pt->EWNS ('W', 'W')->Phi ()->Value ();
  }
  // 왼쪽 위의 점과 왼쪽 아래의 점이 존재하는 경우
  else if (pt->EWNS ('W', 'N') != NULL && pt->EWNS ('W', 'S') != NULL) {
    // 왼쪽 위의 점의 y-좌표를 저장
    yp = pt->EWNS ('W', 'N')->Coord ().Value ('y');
    // 왼쪽 아래의 점의 y-좌표를 저장
    ym = pt->EWNS ('W', 'S')->Coord ().Value ('y');
    // 왼쪽점에서의 phi의 값을 근사
    phiW = pt->EWNS ('W', 'N')->Phi ()->Value () * (yb - ym) / (yp - ym) + pt->EWNS ('W', 'S')->Phi ()->Value () * (yp - yb) / (yp - ym);
  }
  // 둘 다 아닌 경우 에러메시지를 출력
  else PrintError ("PhiInterp");

  //  위쪽점이 존재하는 경우
  if (pt->EWNS ('N', 'N') != NULL) {
    phiN = pt->EWNS ('N', 'N')->Phi ()->Value ();
  }
  //  위쪽의 오른쪽 점과  위쪽의 왼쪽점이 존재하는 경우
  else if (pt->EWNS ('N', 'E') != NULL && pt->EWNS ('N', 'W') != NULL) {
    //  위쪽의 오른쪽 점의 x-좌표를 저장
    xp = pt->EWNS ('N', 'E')->Coord ().Value ('x');
    //  위쪽의 왼쪽 점의 x-좌표를 저장
    xm = pt->EWNS ('N', 'W')->Coord ().Value ('x');
    //  위쪽점에서의 phi의 값을 근사
    phiN = pt->EWNS ('N', 'E')->Phi ()->Value () * (xb - xm) / (xp - xm) + pt->EWNS ('N', 'W')->Phi ()->Value () * (xp - xb) / (xp - xm);
  }
  // 둘 다 아닌 경우 에러메시지를 출력
  else PrintError ("PhiInterp");

  // 아래쪽점이 존재하는 경우
  if (pt->EWNS ('S', 'S') != NULL) {
    phiS = pt->EWNS ('S', 'S')->Phi ()->Value ();
  }
  // 아래쪽의 오른쪽 점과 아래쪽의 왼쪽점이 존재하는 경우
  else if (pt->EWNS ('S', 'E') != NULL && pt->EWNS ('S', 'W') != NULL) {
    // 아래쪽의 오른쪽 점의 x-좌표를 저장
    xp = pt->EWNS ('S', 'E')->Coord ().Value ('x');
    // 아래쪽의 왼쪽 점의 x-좌표를 저장
    xm = pt->EWNS ('S', 'W')->Coord ().Value ('x');
    // 아래쪽점에서의 phi의 값을 근사
    phiS = pt->EWNS ('S', 'E')->Phi ()->Value () * (xb - xm) / (xp - xm) + pt->EWNS ('S', 'W')->Phi ()->Value () * (xp - xb) / (xp - xm);
  }
  // 둘 다 아닌 경우 에러메시지를 출력
  else PrintError ("PhiInterp");

  // 국소 x-축선의 시작점과 끝점의 x-좌표와 현재점의 x-좌표
  xm = pt->MinMaxCoordinate ('x', 'm'), xb = pt->Coord ().Value ('x'), xp = pt->MinMaxCoordinate ('x', 'p');
  // 국소 y-축선의 시작점과 끝점의 y-좌표와 현재점의 y-좌표
  ym = pt->MinMaxCoordinate ('y', 'm'), yb = pt->Coord ().Value ('y'), yp = pt->MinMaxCoordinate ('y', 'p');

  // x-축선 위에서 phi의 값의 근사
  phiX = phiE * (xb - xm) / (xp - xm) + phiW * (xp - xb) / (xp - xm);
  // y-축선 위에서 phi의 값의 근사
  phiY = phiN * (yb - ym) / (yp - ym) + phiS * (yp - yb) / (yp - ym);

  // return 5.0E-1 * (phiX + phiY);
  return phiY;
}

void SettingMaterialProperties (Point *pt, double *mpx1, double *mpx2, double *mpy1, double *mpy2) {
  char errorMassage[1024];

  if (pt->Condition () != 'I') {
    sprintf (errorMassage, "Error occurs in SettingMaterialProperties, pt->Condition () == %c, it should be 'I'", pt->Condition ());
    PrintError (errorMassage);
    exit (1);
  }

  // 오른쪽점이 존재한는 경우, 오른쪽점의 conductivity를 국소 x-축선의 끝점의 conductivity에 저장
  if (pt->EWNS ('E', 'E') != NULL) *mpx2 = pt->EWNS ('E', 'E')->MaterialProperty ();
  // 오른쪽의 위쪽점이 존재한는 경우, 오른쪽의 위쪽점의 conductivity를 국소 x-축선의 끝점의 conductivity에 저장
  else if (pt->EWNS ('E', 'N') != NULL) *mpx2 = pt->EWNS ('E', 'N')->MaterialProperty ();
  // 오른쪽의 아래쪽점이 존재한는 경우, 오른쪽의 아래쪽점의 conductivity를 국소 x-축선의 끝점의 conductivity에 저장
  else if (pt->EWNS ('E', 'S') != NULL) *mpx2 = pt->EWNS ('E', 'S')->MaterialProperty ();
  // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
  else PrintError ("Interface has no East Points!!");
  // 왼쪽점이 존재한는 경우, 왼쪽점의 conductivity를 국소 x-축선의 시작점의 conductivity에 저장
  if (pt->EWNS ('W', 'W') != NULL) *mpx1 = pt->EWNS ('W', 'W')->MaterialProperty ();
  // 왼쪽의 위쪽점이 존재한는 경우, 왼쪽의 위쪽점의 conductivity를 국소 x-축선의 시작점의 conductivity에 저장
  else if (pt->EWNS ('W', 'N') != NULL) *mpx1 = pt->EWNS ('W', 'N')->MaterialProperty ();
  // 왼쪽의 아래쪽점이 존재한는 경우, 왼쪽의 아래쪽점의 conductivity를 국소 x-축선의 시작점의 conductivity에 저장
  else if (pt->EWNS ('W', 'S') != NULL) *mpx1 = pt->EWNS ('W', 'S')->MaterialProperty ();
  // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
  else PrintError ("Interface has no West Points!!");
  // 위쪽점이 존재한는 경우, 위쪽점의 conductivity를 국소 y-축선의 끝점의 conductivity에 저장
  if (pt->EWNS ('N', 'N') != NULL) *mpy2 = pt->EWNS ('N', 'N')->MaterialProperty ();
  // 위쪽의 오른쪽점이 존재한는 경우, 위쪽의 오른쪽점의 conductivity를 국소 y-축선의 끝점의 conductivity에 저장
  else if (pt->EWNS ('N', 'E') != NULL) *mpy2 = pt->EWNS ('N', 'E')->MaterialProperty ();
  // 위쪽의 왼쪽점이 존재한는 경우, 위쪽의 왼쪽점의 conductivity를 국소 y-축선의 끝점의 conductivity에 저장
  else if (pt->EWNS ('N', 'W') != NULL) *mpy2 = pt->EWNS ('N', 'W')->MaterialProperty ();
  // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
  else PrintError ("Interface has no North Points!!");
  // 아래쪽점이 존재한는 경우, 아래쪽점의 conductivity를 국소 y-축선의 시작점의 conductivity에 저장
  if (pt->EWNS ('S', 'S') != NULL) *mpy1 = pt->EWNS ('S', 'S')->MaterialProperty ();
  // 아래쪽의 오른쪽점이 존재한는 경우, 아래쪽의 오른쪽점의 conductivity를 국소 y-축선의 시작점의 conductivity에 저장
  else if (pt->EWNS ('S', 'E') != NULL) *mpy1 = pt->EWNS ('S', 'E')->MaterialProperty ();
  // 아래쪽의 왼쪽점이 존재한는 경우, 아래쪽의 왼쪽점의 conductivity를 국소 y-축선의 시작점의 conductivity에 저장
  else if (pt->EWNS ('S', 'W') != NULL) *mpy1 = pt->EWNS ('S', 'W')->MaterialProperty ();
  // 모두 존재하지 않는 경우, 에러메시지를 출력하고 종료
  else PrintError ("Interface has no South Points!!");
}

void InitializationCoef (xData *xdat, yData *ydat, string xy) {
  if ((xy.find("x") != string::npos) || (xy.find("X") != string::npos)) {
    // x-축선에서의 해의 표현식의 계수들을 0으로 초기화
    xdat->F  = ZeroValue;
    xdat->Cu = ZeroValue, xdat->Cphi = ZeroValue; xdat->Cf = ZeroValue, xdat->Cdu = ZeroValue;
    xdat->Eu = ZeroValue, xdat->ENu = ZeroValue, xdat->ESu = ZeroValue;
    xdat->Wu = ZeroValue, xdat->WNu = ZeroValue, xdat->WSu = ZeroValue;
    xdat->Edu = ZeroValue, xdat->ENdu = ZeroValue, xdat->ESdu = ZeroValue;
    xdat->Wdu = ZeroValue, xdat->WNdu = ZeroValue, xdat->WSdu = ZeroValue;
    xdat->Ephi = ZeroValue, xdat->ENphi = ZeroValue, xdat->ESphi = ZeroValue;
    xdat->Wphi = ZeroValue, xdat->WNphi = ZeroValue, xdat->WSphi = ZeroValue;
    xdat->Ef = ZeroValue, xdat->ENf = ZeroValue, xdat->ESf = ZeroValue;
    xdat->Wf = ZeroValue, xdat->WNf = ZeroValue, xdat->WSf = ZeroValue;
  }

  if ((xy.find("y") != string::npos) || (xy.find("Y") != string::npos)) {
    // y-축선에서의 해의 표현식의 계수들을 0으로 초기화
    ydat->F = ZeroValue;
    ydat->Cu = ZeroValue, ydat->Cphi = ZeroValue; ydat->Cf = ZeroValue, ydat->Cdu = ZeroValue;
    ydat->Nu = ZeroValue, ydat->NEu = ZeroValue, ydat->NWu = ZeroValue;
    ydat->Su = ZeroValue, ydat->SEu = ZeroValue, ydat->SWu = ZeroValue;
    ydat->Ndu = ZeroValue, ydat->NEdu = ZeroValue, ydat->NWdu = ZeroValue;
    ydat->Sdu = ZeroValue, ydat->SEdu = ZeroValue, ydat->SWdu = ZeroValue;
    ydat->Nphi = ZeroValue, ydat->NEphi = ZeroValue, ydat->NWphi = ZeroValue;
    ydat->Sphi = ZeroValue, ydat->SEphi = ZeroValue, ydat->SWphi = ZeroValue;
    ydat->Nf = ZeroValue, ydat->NEf = ZeroValue, ydat->NWf = ZeroValue;
    ydat->Sf = ZeroValue, ydat->SEf = ZeroValue, ydat->SWf = ZeroValue;
  }
}

void SettingGreenShape (Point *pt, int *bdx, int *bdy) {
  char azimuthNX[2][3] = {"E2", "W1"};
  char azimuthNY[2][3] = {"N2", "S1"};
  // char azimuthF[4][3] = {"E4", "W3", "N4", "S3"};
  if (pt->Condition () == 'N') {
    for (const auto &i : azimuthNX) if (pt->IsBoundary (i[0])) *bdx = int (i[1]) - '0';
    for (const auto &i : azimuthNY) if (pt->IsBoundary (i[0])) *bdy = int (i[1]) - '0';
  }
  // for (const auto &i : azimuthF) if (pt->EWNS (i[0], i[0])) if (pt->EWNS (i[0], i[0])->Condition () == 'F') *bdx = int (i[1]) - '0';
}

void SettingCoefficient (Point *pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2,
  const int nD, const bool is_sol) {
    double  sign = is_sol ? 1.0E0 : -1.0E0;
    // 국소 x-축선에서 u의 현재점에서의 계수
    xdat->Cu = - UnitValue;
    // 국소 x-축선에서 u의 국소 x-축선에서 끝점에서의 계수
    xdat->Eu = greens_coefficient_t (xp, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 u의 국소 x-축선에서 시작점에서의 계수
    xdat->Wu = greens_coefficient_t (xm, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 phi의 현재점에서의 계수
    xdat->Cphi = greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) + greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 phi의 국소 x-축선에서 끝점에서의 계수
    xdat->Ephi = greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 phi의 국소 x-축선에서 시작점에서의 계수
    xdat->Wphi = greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 f의 현재점에서의 계수
    xdat->Cf = greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) + greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 f의 국소 x-축선에서 끝점에서의 계수
    xdat->Ef = greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 f의 국소 x-축선에서 시작점에서의 계수
    xdat->Wf = greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // x-축선의 양끝점중 한 개의 점의 경계조건이 Infinity 혹은 Singularity일 때, 정확한 적분을 계산
    ExactIntegration ('x', xdat, ydat, nD, bdx, bdy, xm, xb, xp, ym, yb, yp);
    // 국소 y-축선에서 u의 현재점에서의 계수
    ydat->Cu = - sign * UnitValue;
    // 국소 y-축선에서 u의 국소 y-축선에서 끝점에서의 계수
    ydat->Nu = sign * greens_coefficient_t (yp, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 u의 국소 y-축선에서 시작점에서의 계수
    ydat->Su = sign * greens_coefficient_t (ym, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 phi의 현재점에서의 계수
    ydat->Cphi = - sign * greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) - sign * greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 phi의 국소 y-축선에서 끝점에서의 계수
    ydat->Nphi = - sign * greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 phi의 국소 y-축선에서 시작점에서의 계수
    ydat->Sphi = - sign * greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 f의 현재점에서의 계수
    ydat->Cf = sign * greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) + sign * greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 f의 국소 y-축선에서 끝점에서의 계수
    ydat->Nf = sign * greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 f의 국소 y-축선에서 시작점에서의 계수
    ydat->Sf = sign * greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // y-축선의 양끝점중 한 개의 점의 경계조건이 Infinity 혹은 Singularity일 때, 정확한 적분을 계산
    ExactIntegration ('y', xdat, ydat, nD, bdx, bdy, xm, xb, xp, ym, yb, yp);

    // For Heat Equation
    if (::SolverType) SettingTimeIntegration (pt, xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1, mpy2, is_sol);
  }

void SettingNeumannCoefficient (Point *pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2,
  const bool is_sol) {
  if (is_sol) {
    // 국소 x-축선에서 u의 국소 x-축선에서 끝점에서의 계수
    xdat->Eu += greens_coefficient_ttau (xp, xm, xb, xp, xb, yb, 1, 0, mpx1, mpx2);
    // 국소 x-축선에서 u의 국소 x-축선에서 시작점에서의 계수
    xdat->Wu += greens_coefficient_ttau (xm, xm, xb, xp, xb, yb, 1, 0, mpx1, mpx2);
    // 국소 x-축선에서 phi의 현재점에서의 계수
    if (!IsEqualDouble (xb, xm)) xdat->Cphi += greens_integral_tau (2, xm, xb, xp, xb, yb, 1, 0, mpx1, mpx2);
    if (!IsEqualDouble (xb, xp)) xdat->Cphi += greens_integral_tau (3, xm, xb, xp, xb, yb, 1, 0, mpx1, mpx2);
    // 국소 x-축선에서 phi의 국소 x-축선에서 끝점에서의 계수
    if (!IsEqualDouble (xb, xp)) xdat->Ephi += greens_integral_tau (4, xm, xb, xp, xb, yb, 1, 0, mpx1, mpx2);
    // 국소 x-축선에서 phi의 국소 x-축선에서 시작점에서의 계수
    if (!IsEqualDouble (xb, xm)) xdat->Wphi += greens_integral_tau (1, xm, xb, xp, xb, yb, 1, 0, mpx1, mpx2);
    // 국소 x-축선에서 f의 현재점에서의 계수
    if (!IsEqualDouble (xb, xm)) xdat->Cf += greens_integral_tau (2, xm, xb, xp, xb, yb, 1, 0, mpx1, mpx2);
    if (!IsEqualDouble (xb, xp)) xdat->Cf += greens_integral_tau (3, xm, xb, xp, xb, yb, 1, 0, mpx1, mpx2);
    // 국소 x-축선에서 f의 국소 x-축선에서 끝점에서의 계수
    if (!IsEqualDouble (xb, xp)) xdat->Ef += greens_integral_tau (4, xm, xb, xp, xb, yb, 1, 0, mpx1, mpx2);
    // 국소 x-축선에서 f의 국소 x-축선에서 시작점에서의 계수
    if (!IsEqualDouble (xb, xm)) xdat->Wf += greens_integral_tau (1, xm, xb, xp, xb, yb, 1, 0, mpx1, mpx2);
    // 국소 y-축선에서 u의 국소 y-축선에서 끝점에서의 계수
    ydat->Nu += greens_coefficient_ttau (yp, ym, yb, yp, xb, yb, 2, 0, mpy1, mpy2);
    // 국소 y-축선에서 u의 국소 y-축선에서 시작점에서의 계수
    ydat->Su += greens_coefficient_ttau (ym, ym, yb, yp, xb, yb, 2, 0, mpy1, mpy2);
    // 국소 y-축선에서 phi의 현재점에서의 계수
    if (!IsEqualDouble (yb, ym)) ydat->Cphi -= greens_integral_tau (2, ym, yb, yp, xb, yb, 2, 0, mpy1, mpy2);
    if (!IsEqualDouble (yb, yp)) ydat->Cphi -= greens_integral_tau (3, ym, yb, yp, xb, yb, 2, 0, mpy1, mpy2);
    // 국소 y-축선에서 phi의 국소 y-축선에서 끝점에서의 계수
    if (!IsEqualDouble (yb, yp)) ydat->Nphi -= greens_integral_tau (4, ym, yb, yp, xb, yb, 2, 0, mpy1, mpy2);
    // 국소 y-축선에서 phi의 국소 y-축선에서 시작점에서의 계수
    if (!IsEqualDouble (yb, ym)) ydat->Sphi -= greens_integral_tau (1, ym, yb, yp, xb, yb, 2, 0, mpy1, mpy2);
    // 국소 y-축선에서 f의 현재점에서의 계수
    if (!IsEqualDouble (yb, ym)) ydat->Cf += greens_integral_tau (2, ym, yb, yp, xb, yb, 2, 0, mpy1, mpy2);
    if (!IsEqualDouble (yb, yp)) ydat->Cf += greens_integral_tau (3, ym, yb, yp, xb, yb, 2, 0, mpy1, mpy2);
    // 국소 y-축선에서 f의 국소 y-축선에서 끝점에서의 계수
    if (!IsEqualDouble (yb, yp)) ydat->Nf += greens_integral_tau (4, ym, yb, yp, xb, yb, 2, 0, mpy1, mpy2);
    // 국소 y-축선에서 f의 국소 y-축선에서 시작점에서의 계수
    if (!IsEqualDouble (yb, ym)) ydat->Sf += greens_integral_tau (1, ym, yb, yp, xb, yb, 2, 0, mpy1, mpy2);

    // For Heat Equation
    if (::SolverType) SettingNeumannTimeIntegration (pt, xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1, mpy2, is_sol);

    for (const auto &i : {&xdat->Cu, &xdat->Eu, &xdat->Wu, &xdat->Cphi, &xdat->Ephi, &xdat->Wphi, &xdat->Cf, &xdat->Ef, &xdat->Wf}) if (IsEqualDouble (pt->Normal ().Value ('x'), ZeroValue)) *i = ZeroValue; else *i *= pt->Normal ().Value ('x');
    for (const auto &i : {&ydat->Cu, &ydat->Nu, &ydat->Su, &ydat->Cphi, &ydat->Nphi, &ydat->Sphi, &ydat->Cf, &ydat->Nf, &ydat->Sf}) if (IsEqualDouble (pt->Normal ().Value ('y'), ZeroValue)) *i = ZeroValue; else *i *= pt->Normal ().Value ('y');
  }

  else {
    PrintError ("SettingNeumannCoefficient");


    double XM = xm, XB = xb, XP = xp;
    double YM = ym, YB = yb, YP = yp;

    unordered_map<char, double*> m;
    unordered_map<char, char> opposite_azimuth;

    opposite_azimuth['E'] = 'W', opposite_azimuth['W'] = 'E', opposite_azimuth['N'] = 'S', opposite_azimuth['S'] = 'N';
    m['E'] = &XP, m['W'] = &XM, m['N'] = &YP, m['S'] = &YM;

    bool sh = true;

    for (const auto &i : {'N', 'S'})
    if (pt->IsBoundary (opposite_azimuth[i]))
    for (const auto &j : {i, 'E', 'W'})
    if (pt->EWNS2nd (i, j)) *m[i] = pt->EWNS2nd (i, j)->Coord ().Value ('y'), sh = false;

    if (IsEqualDouble (XP - XM, YP - YM) || sh)
    for (const auto &i : {'E', 'W'})
    if (pt->IsBoundary (opposite_azimuth[i]))
    for (const auto &j : {i, 'N', 'S'})
    if (pt->EWNS2nd (i, j)) *m[i] = pt->EWNS2nd (i, j)->Coord ().Value ('x');

    SettingNeumannCoefficient (pt->Pressure (), xdat, ydat, XM, XB, XP, YM, YB, YP, bdx, bdy, mpx1, mpx2, mpy1, mpy2, true);
  }
}

void SettingDirichletCoefficient (Point *pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2,
  const int nD, const bool is_sol) {

    if (is_sol) ydat->Cu = 1.0E0, ydat->F = pt->Boundaryvalue ();
    else {
      double XM = xm, XB = xb, XP = xp;
      double YM = ym, YB = yb, YP = yp;

      char errorMassage[256];

      Point *calc_pt = NULL;

      unordered_map<char, double*> m;
      unordered_map<char, char> opposite_azimuth;
      unordered_map<char, char> coord;
      unordered_map<char, string> axis;

      m['E'] = &XM, m['W'] = &XP, m['N'] = &YM, m['S'] = &YP;
      opposite_azimuth['E'] = 'W', opposite_azimuth['W'] = 'E', opposite_azimuth['N'] = 'S', opposite_azimuth['S'] = 'N';
      coord['E'] = 'x', coord['W'] = 'x', coord['N'] = 'y', coord['S'] = 'y';
      axis['x'] = "Y", axis['y'] = "X";

      for (const auto &i : {'E', 'W', 'N', 'S'})
      if (pt->Axis (coord[i]) > -1)
      if (pt->IsBoundary (opposite_azimuth[i]))
      if (pt->EWNS (i, i))
      if (pt->EWNS (i, i)->EWNS (i, i))
      if (pt->EWNS (i, i)->EWNS (i, i)->Condition () == 'C') {
        calc_pt = pt->EWNS (i, i)->EWNS (i, i);
        XM = calc_pt->Pressure ()->MinMaxCoordinate ('x', 'm'), XB = calc_pt->Pressure ()->Coord ().Value ('x'), XP = calc_pt->Pressure ()->MinMaxCoordinate ('x', 'p');
        YM = calc_pt->Pressure ()->MinMaxCoordinate ('y', 'm'), YB = calc_pt->Pressure ()->Coord ().Value ('y'), YP = calc_pt->Pressure ()->MinMaxCoordinate ('y', 'p');
        *m[i] = pt->Pressure ()->Coord ().Value (coord[i]);
        break;
      }

      if (!calc_pt) sprintf (errorMassage, "SettingDirichletCoefficient, calc_pt = NULL point"), pt->Velocity ('u')->PrintDebuggingData ("11110", xdat, ydat, true), PrintError (errorMassage);

      SettingCoefficient (calc_pt->Pressure (), xdat, ydat, XM, XB, XP, YM, YB, YP, 0, 0, mpx1, mpx2, mpy1, mpy2, nD, true);

      for (const auto &i : {'E', 'W', 'N', 'S'})
      if (pt->Axis (coord[i]) > -1)
      if (pt->IsBoundary (opposite_azimuth[i]))
      if (pt->EWNS (i, i))
      if (pt->EWNS (i, i)->EWNS (i, i))
      if (pt->EWNS (i, i)->EWNS (i, i)->Condition () == 'C') {InitializationCoef (xdat, ydat, axis[coord[i]]); break;}

      RightHandSideDirichlet (pt, calc_pt, xdat, ydat);

      // For Heat Equation
      if (::SolverType) SettingDirichletTimeIntegrationRightHand (pt, calc_pt, xdat, ydat, XM, XB, XP, YM, YB, YP, bdx, bdy, mpx1, mpx2, mpy1, mpy2, is_sol);

      if (pt->Pressure ()->Axis ('x') < 0 && pt->Pressure ()->Axis ('y') < 0) sprintf (errorMassage, "SettingDirichletCoefficient, all axis < 0"), PrintError (errorMassage);
    }
  }

void CalcInterfaceCoef (Point *pt, xData *xdat, yData *ydat, const char *azimuthHV, const double xb, const double yb, const double tm, const double tp, const double mp) {
  double (*uCoef) (char, double, double, double, double, int, double);
  double (*phiCoef) (char, double, double, double, double, int, double);
  double (*fCoef) (char, double, double, double, double, int, double);

  double t[2];
  double tt;
  double *u[3], *phi[3], *f[4];
  int bd = 0;

  char c;
  char sub_azimuth[2];
  char errorMassage[1024];

  char azimuth = azimuthHV[0];
  char HV      = azimuthHV[1];

  if (azimuth == 'E' || azimuth == 'W') {
    tt = yb, c = 'y', sub_azimuth[0] = 'N', sub_azimuth[1] = 'S', f[3] = &xdat->F;
    if (azimuth == 'E') {
      u[0] = &xdat->Eu, u[1] = &xdat->ENu, u[2] = &xdat->ESu;
      phi[0] = &xdat->Ephi, phi[1] = &xdat->ENphi, phi[2] = &xdat->ESphi;
      f[0] = &xdat->Ef, f[1] = &xdat->ENf, f[2] = &xdat->ESf;
    }
    if (azimuth == 'W') {
      u[0] = &xdat->Wu, u[1] = &xdat->WNu, u[2] = &xdat->WSu;
      phi[0] = &xdat->Wphi, phi[1] = &xdat->WNphi, phi[2] = &xdat->WSphi;
      f[0] = &xdat->Wf, f[1] = &xdat->WNf, f[2] = &xdat->WSf;
    }
  }
  else if (azimuth == 'N' || azimuth == 'S') {
    tt = xb, c = 'x', sub_azimuth[0] = 'E', sub_azimuth[1] = 'W', f[3] = &ydat->F;
    if (azimuth == 'N') {
      u[0] = &ydat->Nu, u[1] = &ydat->NEu, u[2] = &ydat->NWu;
      phi[0] = &ydat->Nphi, phi[1] = &ydat->NEphi, phi[2] = &ydat->NWphi;
      f[0] = &ydat->Nf, f[1] = &ydat->NEf, f[2] = &ydat->NWf;
    }
    if (azimuth == 'S') {
      u[0] = &ydat->Su, u[1] = &ydat->SEu, u[2] = &ydat->SWu;
      phi[0] = &ydat->Sphi, phi[1] = &ydat->SEphi, phi[2] = &ydat->SWphi;
      f[0] = &ydat->Sf, f[1] = &ydat->SEf, f[2] = &ydat->SWf;
    }
  }
  else {
    sprintf (errorMassage, "CalcInterfaceCoef, azimuth = %c, please check azimuth", azimuth);
    PrintError (errorMassage);
  }

  if (HV == 'V' || HV == 'v') {
    uCoef = CalcVerticalUCoefficient;
    phiCoef = CalcVerticalPHICoefficient;
    fCoef = CalcVerticalFCoefficient;
  }

  else if (HV == 'H' || HV == 'h') {
    uCoef = CalcHorizontalUCoefficient;
    phiCoef = CalcHorizontalPHICoefficient;
    fCoef = CalcHorizontalFCoefficient;
  }

  else {
    sprintf (errorMassage, "CalcInterfaceCoef, HV = %c, please check HV", HV);
    PrintError (errorMassage);
  }

  if (pt->EWNS (azimuth, azimuth) == NULL) {
    for (size_t i = 0; i < 2; i++) {
      if (pt->EWNS (azimuth, sub_azimuth[i]) == NULL) {
        sprintf (errorMassage, "CalcInterfaceCoef, (%c, %c) point is NULL point but (%c, %c) point is also NULL point", azimuth, azimuth, azimuth, sub_azimuth[i]);
        PrintError (errorMassage);
      }
    }

    for (size_t i = 0; i < 2; i++) if (pt->EWNS (azimuth, sub_azimuth[i])->Condition () == 'F') bd = 4 - i;
    for (size_t i = 0; i < 2; i++) t[i] = pt->EWNS (azimuth, sub_azimuth[i])->Coord ().Value (c);
    for (size_t i = 1; i < 3; i++) {
      *u[i]   =   uCoef (sub_azimuth[i - 1], xb, yb, t[0], t[1], bd, mp) * (*u[0]);
      *phi[i] = phiCoef (sub_azimuth[i - 1], xb, yb, t[0], t[1], bd, mp) * (*u[0]) + (*phi[0]) * fabs (tt - t[i % 2]) / (t[0] - t[1]);
      *f[i]   =   fCoef (sub_azimuth[i - 1], xb, yb, t[0], t[1], bd, mp) * (*u[0]) +   (*f[0]) * fabs (tt - t[i % 2]) / (t[0] - t[1]);

      if (::SolverType) *u[i] -= 5.0E-1 * fCoef (sub_azimuth[i - 1], xb, yb, t[0], t[1], bd, mp) * (*u[0]) / pt->Dt ();
      if (::SolverType) *f[3] +=          fCoef (sub_azimuth[i - 1], xb, yb, t[0], t[1], bd, mp) * (*u[0]) / pt->Dt () * pt->EWNS (azimuth, sub_azimuth[i - 1])->Pre ()->Value ();
    }
  }
}

void ApproximateInterfacePt (Point *pt, xData *xdat, yData *ydat, double xm, double xb, double xp, double ym, double yb, double yp, double mpx1, double mpx2, double mpy1, double mpy2, bool is_sol) {
  char azimuthHV[4][3] = {"EV", "WV", "NH", "SH"};
  double mp[4] = {mpx2, mpx1, mpy2, mpy1};
  unordered_map<char, double> tx, ty;
  tx['E'] = xp, tx['W'] = xm, tx['N'] = xb, tx['S'] = xb;
  ty['E'] = yb, ty['W'] = yb, ty['N'] = yp, ty['S'] = ym;

  for (size_t i = 0; i < 4; i++) CalcInterfaceCoef (pt, xdat, ydat, azimuthHV[i], tx[azimuthHV[i][0]], ty[azimuthHV[i][0]], ym, yp, mp[i]);
}
void CalcDiffDiffCoef (Point *pt, xData *xdat, yData *ydat, const char *azimuthHV, const double xb, const double yb, const double tm, const double tp, const double mp) {
  double (*uCoef) (char, double, double, double, double, int, double);
  double (*phiCoef) (char, double, double, double, double, int, double);
  double (*fCoef) (char, double, double, double, double, int, double);

  double t[2];
  double tt;
  double *u[3], *du[3], *phi[3], *f[4];
  int bd = 0;

  char c;
  char sub_azimuth[2];
  char errorMassage[1024];

  char azimuth = azimuthHV[0];
  char HV      = azimuthHV[1];

  if (azimuth == 'E' || azimuth == 'W') {
    tt = yb, c = 'y', sub_azimuth[0] = 'N', sub_azimuth[1] = 'S', f[3] = &xdat->F;
    if (azimuth == 'E') {
      u[0] = &xdat->Eu, u[1] = &xdat->ENu, u[2] = &xdat->ESu;
      du[0] = &xdat->Edu, du[1] = &xdat->ENdu, du[2] = &xdat->ESdu;
      phi[0] = &xdat->Ephi, phi[1] = &xdat->ENphi, phi[2] = &xdat->ESphi;
      f[0] = &xdat->Ef, f[1] = &xdat->ENf, f[2] = &xdat->ESf;
    }
    if (azimuth == 'W') {
      u[0] = &xdat->Wu, u[1] = &xdat->WNu, u[2] = &xdat->WSu;
      du[0] = &xdat->Wdu, du[1] = &xdat->WNdu, du[2] = &xdat->WSdu;
      phi[0] = &xdat->Wphi, phi[1] = &xdat->WNphi, phi[2] = &xdat->WSphi;
      f[0] = &xdat->Wf, f[1] = &xdat->WNf, f[2] = &xdat->WSf;
    }
  }
  else if (azimuth == 'N' || azimuth == 'S') {
    tt = xb, c = 'x', sub_azimuth[0] = 'E', sub_azimuth[1] = 'W', f[3] = &ydat->F;
    if (azimuth == 'N') {
      u[0] = &ydat->Nu, u[1] = &ydat->NEu, u[2] = &ydat->NWu;
      du[0] = &ydat->Ndu, du[1] = &ydat->NEdu, du[2] = &ydat->NWdu;
      phi[0] = &ydat->Nphi, phi[1] = &ydat->NEphi, phi[2] = &ydat->NWphi;
      f[0] = &ydat->Nf, f[1] = &ydat->NEf, f[2] = &ydat->NWf;
    }
    if (azimuth == 'S') {
      u[0] = &ydat->Su, u[1] = &ydat->SEu, u[2] = &ydat->SWu;
      du[0] = &ydat->Sdu, du[1] = &ydat->SEdu, du[2] = &ydat->SWdu;
      phi[0] = &ydat->Sphi, phi[1] = &ydat->SEphi, phi[2] = &ydat->SWphi;
      f[0] = &ydat->Sf, f[1] = &ydat->SEf, f[2] = &ydat->SWf;
    }
  }
  else {
    sprintf (errorMassage, "CalcDiffDiffCoef, azimuth = %c, please check azimuth", azimuth);
    PrintError (errorMassage);
  }

  if (HV == 'V' || HV == 'v') {
    uCoef = CalcVerticalDiffUCoefficient;
    phiCoef = CalcVerticalDiffPHICoefficient;
    fCoef = CalcVerticalDiffFCoefficient;
  }

  else if (HV == 'H' || HV == 'h') {
    uCoef = CalcHorizontalDiffUCoefficient;
    phiCoef = CalcHorizontalDiffPHICoefficient;
    fCoef = CalcHorizontalDiffFCoefficient;
  }

  else {
    sprintf (errorMassage, "CalcDiffDiffCoef, HV = %c, please check HV", HV);
    PrintError (errorMassage);
  }

  if (pt->EWNS (azimuth, azimuth) == NULL) {
    for (size_t i = 0; i < 2; i++) {
      if (pt->EWNS (azimuth, sub_azimuth[i]) == NULL) {
        sprintf (errorMassage, "CalcDiffDiffCoef, (%c, %c) point is NULL point but (%c, %c) point is also NULL point", azimuth, azimuth, azimuth, sub_azimuth[i]);
        PrintError (errorMassage);
      }
    }

    for (size_t i = 0; i < 2; i++) if (pt->EWNS (azimuth, sub_azimuth[i])->Condition () == 'F') bd = 4 - i;
    for (size_t i = 0; i < 2; i++) t[i] = pt->EWNS (azimuth, sub_azimuth[i])->Coord ().Value (c);
    for (size_t i = 1; i < 3; i++) {
      *du[i] = (*du[0]) * fabs (tt - t[i % 2]) / (t[0] - t[1]);
      // *u[i]   +=   uCoef (sub_azimuth[i - 1], xb, yb, t[0], t[1], bd, mp) * (*du[0]);
      // *phi[i] += phiCoef (sub_azimuth[i - 1], xb, yb, t[0], t[1], bd, mp) * (*du[0]);
      // *f[i]   +=   fCoef (sub_azimuth[i - 1], xb, yb, t[0], t[1], bd, mp) * (*du[0]);
      //
      // if (::SolverType) *u[i] -= 5.0E-1 * fCoef (sub_azimuth[i - 1], xb, yb, t[0], t[1], bd, mp) * (*du[0]) / pt->Dt ();
      // if (::SolverType) *f[3] +=          fCoef (sub_azimuth[i - 1], xb, yb, t[0], t[1], bd, mp) * (*du[0]) / pt->Dt () * pt->EWNS (azimuth, sub_azimuth[i - 1])->Pre ()->Value ();
    }
  }
}


void ApproximateDiffDiffPt (Point *pt, xData *xdat, yData *ydat, double xm, double xb, double xp, double ym, double yb, double yp, double mpx1, double mpx2, double mpy1, double mpy2, bool is_sol) {
  char azimuthHV[4][3] = {"EV", "WV", "NH", "SH"};
  double mp[4] = {mpx2, mpx1, mpy2, mpy1};
  unordered_map<char, double> tx, ty;
  tx['E'] = xp, tx['W'] = xm, tx['N'] = xb, tx['S'] = xb;
  ty['E'] = yb, ty['W'] = yb, ty['N'] = yp, ty['S'] = ym;

  for (size_t i = 0; i < 4; i++) CalcInterfaceCoef (pt, xdat, ydat, azimuthHV[i], tx[azimuthHV[i][0]], ty[azimuthHV[i][0]], ym, yp, mp[i]),
                                 CalcDiffDiffCoef  (pt, xdat, ydat, azimuthHV[i], tx[azimuthHV[i][0]], ty[azimuthHV[i][0]], ym, yp, mp[i]);

  // if (*pt == 942) PrintError ("-----");
}

void TransposeNeumannBoundaryData (Point *pt, xData *xdat, yData *ydat) {
  // double *phi[2];
  // char errorMassage[1024];
  // char azimuth[4] = {'E', 'W', 'N', 'S'};

  ydat->F += pt->Boundaryvalue ();

  // for (const auto &j : azimuth) {
  //   if      (j == 'E') phi[0] = &xdat->Cphi, phi[1] = &xdat->Ephi;
  //   else if (j == 'W') phi[0] = &xdat->Cphi, phi[1] = &xdat->Wphi;
  //   else if (j == 'N') phi[0] = &ydat->Cphi, phi[1] = &ydat->Nphi;
  //   else if (j == 'S') phi[0] = &ydat->Cphi, phi[1] = &ydat->Sphi;
  //   else sprintf (errorMassage, "TransposeNeumannBoundaryData, azimuth = %c, please check azimuth\n", j), PrintError (errorMassage), exit (1);
  //
  //   if (pt->IsBoundary (j)) *phi[0] += *phi[1], *phi[1] = ZeroValue;
  // }
}

void TransposeOtherBoundaryData (Point *pt, xData *xdat, yData *ydat) {
  double *phi[2];
  double *psi[2];
  double *u = NULL;
  double *v[2];
  double *f = NULL;
  char sub_azimuth[2];
  char errorMassage[1024];
  char azimuth[4] = {'E', 'W', 'N', 'S'};

  for (const auto &j : azimuth) {
    if      (j == 'E') phi[0] = &xdat->Cphi, phi[1] = &xdat->Ephi, u = &xdat->Eu, f = &xdat->F, sub_azimuth[0] = 'N', sub_azimuth[1] = 'S', v[0] = &xdat->ENu, v[1] = &xdat->ESu, psi[0] = &xdat->ENphi, psi[1] = &xdat->ESphi;
    else if (j == 'W') phi[0] = &xdat->Cphi, phi[1] = &xdat->Wphi, u = &xdat->Wu, f = &xdat->F, sub_azimuth[0] = 'N', sub_azimuth[1] = 'S', v[0] = &xdat->WNu, v[1] = &xdat->WSu, psi[0] = &xdat->WNphi, psi[1] = &xdat->WSphi;
    else if (j == 'N') phi[0] = &ydat->Cphi, phi[1] = &ydat->Nphi, u = &ydat->Nu, f = &ydat->F, sub_azimuth[0] = 'E', sub_azimuth[1] = 'W', v[0] = &ydat->NEu, v[1] = &ydat->NWu, psi[0] = &ydat->NEphi, psi[1] = &ydat->NWphi;
    else if (j == 'S') phi[0] = &ydat->Cphi, phi[1] = &ydat->Sphi, u = &ydat->Su, f = &ydat->F, sub_azimuth[0] = 'E', sub_azimuth[1] = 'W', v[0] = &ydat->SEu, v[1] = &ydat->SWu, psi[0] = &ydat->SEphi, psi[1] = &ydat->SWphi;
    else sprintf (errorMassage, "TransposeOtherBoundaryData, azimuth = %c, please check azimuth\n", j), PrintError (errorMassage), exit (1);

    if (pt->EWNS (j, j) != NULL) {
      // if      (pt->EWNS (j, j)->Condition () == 'D') *f -= *u * pt->EWNS (j, j)->Boundaryvalue (), *u = ZeroValue, *phi[0] += *phi[1], *phi[1] = ZeroValue;
      // if      (pt->EWNS (j, j)->Condition () == 'D')
      if      (pt->EWNS (j, j)->Condition () == 'I') *phi[0] += *phi[1], *phi[1] = ZeroValue;
      else if (pt->EWNS (j, j)->Condition () == 'F') *u = ZeroValue, *phi[1] = ZeroValue;
    }
    else {
      for (size_t i = 0; i < 2; i++) if (pt->EWNS (j, sub_azimuth[i]) == NULL) sprintf (errorMassage, "TransposeOtherBoundaryData, (%c %c) point does not exist!!\n", j, sub_azimuth[i]);
      for (size_t i = 0; i < 2; i++) {
        // if (pt->EWNS (j, sub_azimuth[i])->Condition () == 'D') *f -= *v[i] * pt->EWNS (j, sub_azimuth[i])->Boundaryvalue (), *v[i] = ZeroValue, *psi[(i + 1) % 2] += *psi[i], *psi[i] = ZeroValue;
        if (pt->EWNS (j, sub_azimuth[i])->Condition () == 'I') *psi[(i + 1) % 2] += *psi[i], *psi[i] = ZeroValue;
        if (pt->EWNS (j, sub_azimuth[i])->Condition () == 'F') *v[i] = ZeroValue, *psi[i] = ZeroValue;
      }
    }
  }
}

void SettingTimeIntegration (Point *pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2,
  const bool is_sol) {
    double  sign = is_sol ? 1.0E0 : -1.0E0;

    xdat->Cu -= 5.0E-1 *        greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Cu -= 5.0E-1 *        greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Eu -= 5.0E-1 *        greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Wu -= 5.0E-1 *        greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();

    ydat->Cu -= 5.0E-1 * sign * greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Cu -= 5.0E-1 * sign * greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Nu -= 5.0E-1 * sign * greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Su -= 5.0E-1 * sign * greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
}

void SettingDiffTimeIntegation (Point *pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2) {

    xdat->Cu -= 5.0E-1 * greens_integral_tau (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Cu -= 5.0E-1 * greens_integral_tau (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Eu -= 5.0E-1 * greens_integral_tau (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Wu -= 5.0E-1 * greens_integral_tau (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();

    ydat->Cu -= 5.0E-1 * greens_integral_tau (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Cu -= 5.0E-1 * greens_integral_tau (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Nu -= 5.0E-1 * greens_integral_tau (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Su -= 5.0E-1 * greens_integral_tau (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
}

void SettingDiffDiffTimeIntegation (Point *pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2) {

    xdat->Cu += 5.0E-1 * greens_integral_ttau (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Cu += 5.0E-1 * greens_integral_ttau (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Eu += 5.0E-1 * greens_integral_ttau (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Wu += 5.0E-1 * greens_integral_ttau (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();

    ydat->Cu += 5.0E-1 * greens_integral_ttau (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Cu += 5.0E-1 * greens_integral_ttau (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Nu += 5.0E-1 * greens_integral_ttau (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Su += 5.0E-1 * greens_integral_ttau (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
}

void SettingNeumannTimeIntegration (Point *pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2,
  const bool is_sol) {

    xdat->Cu -= 5.0E-1 * greens_integral_tau (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Cu -= 5.0E-1 * greens_integral_tau (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Eu -= 5.0E-1 * greens_integral_tau (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    xdat->Wu -= 5.0E-1 * greens_integral_tau (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();

    ydat->Cu -= 5.0E-1 * greens_integral_tau (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Cu -= 5.0E-1 * greens_integral_tau (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Nu -= 5.0E-1 * greens_integral_tau (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    ydat->Su -= 5.0E-1 * greens_integral_tau (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
}

void SettingTimeIntegrationRightHand (Point *pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2,
  const bool is_sol) {
    double  sign = is_sol ? 1.0E0 : -1.0E0;

    unordered_map <char, double> fe;
    unordered_map <char, double> fw;
    unordered_map <char, double> fn;
    unordered_map <char, double> fs;
    unordered_map <char, unordered_map <char, double>> f;

    char errorMassage[256];
    unordered_map <char, double*> ff;

    f['E'] = fe, f['W'] = fw, f['N'] = fn, f['S'] = fs;
    f['E']['E'] = - 5.0E-1 *        greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    f['W']['W'] = - 5.0E-1 *        greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    f['N']['N'] = - 5.0E-1 * sign * greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    f['S']['S'] = - 5.0E-1 * sign * greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();

    if (pt->EWNS ('E', 'N') && pt->EWNS ('E', 'S'))
    f['E']['N'] = f['E']['E'] * (yb - pt->EWNS ('E', 'S')->Coord ().Value ('y')) / (pt->EWNS ('E', 'N')->Coord ().Value ('y') - pt->EWNS ('E', 'S')->Coord ().Value ('y')),
    f['E']['S'] = f['E']['E'] * (pt->EWNS ('E', 'N')->Coord ().Value ('y') - yb) / (pt->EWNS ('E', 'N')->Coord ().Value ('y') - pt->EWNS ('E', 'S')->Coord ().Value ('y'));

    if (pt->EWNS ('W', 'N') && pt->EWNS ('W', 'S'))
    f['W']['N'] = f['W']['W'] * (yb - pt->EWNS ('W', 'S')->Coord ().Value ('y')) / (pt->EWNS ('W', 'N')->Coord ().Value ('y') - pt->EWNS ('W', 'S')->Coord ().Value ('y')),
    f['W']['S'] = f['W']['W'] * (pt->EWNS ('W', 'N')->Coord ().Value ('y') - yb) / (pt->EWNS ('W', 'N')->Coord ().Value ('y') - pt->EWNS ('W', 'S')->Coord ().Value ('y'));

    if (pt->EWNS ('N', 'E') && pt->EWNS ('N', 'W'))
    f['N']['E'] = f['N']['N'] * (xb - pt->EWNS ('N', 'W')->Coord ().Value ('x')) / (pt->EWNS ('N', 'E')->Coord ().Value ('x') - pt->EWNS ('N', 'W')->Coord ().Value ('x')),
    f['N']['W'] = f['N']['N'] * (pt->EWNS ('N', 'E')->Coord ().Value ('x') - xb) / (pt->EWNS ('N', 'E')->Coord ().Value ('x') - pt->EWNS ('N', 'W')->Coord ().Value ('x'));

    if (pt->EWNS ('S', 'E') && pt->EWNS ('S', 'W'))
    f['S']['E'] = f['S']['S'] * (xb - pt->EWNS ('S', 'W')->Coord ().Value ('x')) / (pt->EWNS ('S', 'E')->Coord ().Value ('x') - pt->EWNS ('S', 'W')->Coord ().Value ('x')),
    f['S']['W'] = f['S']['S'] * (pt->EWNS ('S', 'E')->Coord ().Value ('x') - xb) / (pt->EWNS ('S', 'E')->Coord ().Value ('x') - pt->EWNS ('S', 'W')->Coord ().Value ('x'));

    xdat->F -= 5.0E-1 *        greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt () * pt->Pre ()->Value ();
    xdat->F -= 5.0E-1 *        greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt () * pt->Pre ()->Value ();

    ydat->F -= 5.0E-1 * sign * greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt () * pt->Pre ()->Value ();
    ydat->F -= 5.0E-1 * sign * greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt () * pt->Pre ()->Value ();

    for (const auto &i : {'E', 'W'})
    if (pt->EWNS (i, i)) xdat->F += f[i][i] * pt->EWNS (i, i)->Pre ()->Value ();
    else
    for (const auto &j : {'N', 'S'})
    if (pt->EWNS (i, j)) xdat->F += f[i][j] * pt->EWNS (i, j)->Pre ()->Value ();
    else sprintf (errorMassage, "SettingTimeIntegrationRightHand, pt->EWNS (%c, %c) does not exist, please check azimuth", i, j), PrintError (errorMassage);

    for (const auto &i : {'N', 'S'})
    if (pt->EWNS (i, i)) ydat->F += f[i][i] * pt->EWNS (i, i)->Pre ()->Value ();
    else
    for (const auto &j : {'E', 'W'})
    if (pt->EWNS (i, j)) ydat->F += f[i][j] * pt->EWNS (i, j)->Pre ()->Value ();
    else sprintf (errorMassage, "SettingTimeIntegrationRightHand, pt->EWNS (%c, %c) does not exist, please check azimuth", i, j), PrintError (errorMassage);
}

void SettingDiffTimeIntegrationRightHand (Point *pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2) {

    unordered_map <char, double> fe;
    unordered_map <char, double> fw;
    unordered_map <char, double> fn;
    unordered_map <char, double> fs;
    unordered_map <char, unordered_map <char, double>> f;
    char errorMassage[256];

    f['E'] = fe, f['W'] = fw, f['N'] = fn, f['S'] = fs;
    f['E']['E'] = - 5.0E-1 * greens_integral_tau (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    f['W']['W'] = - 5.0E-1 * greens_integral_tau (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    f['N']['N'] = - 5.0E-1 * greens_integral_tau (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    f['S']['S'] = - 5.0E-1 * greens_integral_tau (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();

    if (pt->EWNS ('E', 'N') && pt->EWNS ('E', 'S'))
    f['E']['N'] = f['E']['E'] * (yb - pt->EWNS ('E', 'S')->Coord ().Value ('y')) / (pt->EWNS ('E', 'N')->Coord ().Value ('y') - pt->EWNS ('E', 'S')->Coord ().Value ('y')),
    f['E']['S'] = f['E']['E'] * (pt->EWNS ('E', 'N')->Coord ().Value ('y') - yb) / (pt->EWNS ('E', 'N')->Coord ().Value ('y') - pt->EWNS ('E', 'S')->Coord ().Value ('y'));

    if (pt->EWNS ('W', 'N') && pt->EWNS ('W', 'S'))
    f['W']['N'] = f['W']['W'] * (yb - pt->EWNS ('W', 'S')->Coord ().Value ('y')) / (pt->EWNS ('W', 'N')->Coord ().Value ('y') - pt->EWNS ('W', 'S')->Coord ().Value ('y')),
    f['W']['S'] = f['W']['W'] * (pt->EWNS ('W', 'N')->Coord ().Value ('y') - yb) / (pt->EWNS ('W', 'N')->Coord ().Value ('y') - pt->EWNS ('W', 'S')->Coord ().Value ('y'));

    if (pt->EWNS ('N', 'E') && pt->EWNS ('N', 'W'))
    f['N']['E'] = f['N']['N'] * (xb - pt->EWNS ('N', 'W')->Coord ().Value ('x')) / (pt->EWNS ('N', 'E')->Coord ().Value ('x') - pt->EWNS ('N', 'W')->Coord ().Value ('x')),
    f['N']['W'] = f['N']['N'] * (pt->EWNS ('N', 'E')->Coord ().Value ('x') - xb) / (pt->EWNS ('N', 'E')->Coord ().Value ('x') - pt->EWNS ('N', 'W')->Coord ().Value ('x'));

    if (pt->EWNS ('S', 'E') && pt->EWNS ('S', 'W'))
    f['S']['E'] = f['S']['S'] * (xb - pt->EWNS ('S', 'W')->Coord ().Value ('x')) / (pt->EWNS ('S', 'E')->Coord ().Value ('x') - pt->EWNS ('S', 'W')->Coord ().Value ('x')),
    f['S']['W'] = f['S']['S'] * (pt->EWNS ('S', 'E')->Coord ().Value ('x') - xb) / (pt->EWNS ('S', 'E')->Coord ().Value ('x') - pt->EWNS ('S', 'W')->Coord ().Value ('x'));

    xdat->F -= 5.0E-1 * greens_integral_tau (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt () * pt->Pre ()->Value ();
    xdat->F -= 5.0E-1 * greens_integral_tau (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt () * pt->Pre ()->Value ();

    ydat->F -= 5.0E-1 * greens_integral_tau (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt () * pt->Pre ()->Value ();
    ydat->F -= 5.0E-1 * greens_integral_tau (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt () * pt->Pre ()->Value ();

    for (const auto &i : {'E', 'W'})
    if (pt->EWNS (i, i)) xdat->F += f[i][i] * pt->EWNS (i, i)->Pre ()->Value ();
    else
    for (const auto &j : {'N', 'S'})
    if (pt->EWNS (i, j)) xdat->F += f[i][j] * pt->EWNS (i, j)->Pre ()->Value ();
    else sprintf (errorMassage, "SettingDiffTimeIntegrationRightHand, pt->EWNS (%c, %c) does not exist, please check azimuth", i, j), PrintError (errorMassage);

    for (const auto &i : {'N', 'S'})
    if (pt->EWNS (i, i)) ydat->F += f[i][i] * pt->EWNS (i, i)->Pre ()->Value ();
    else
    for (const auto &j : {'E', 'W'})
    if (pt->EWNS (i, j)) ydat->F += f[i][j] * pt->EWNS (i, j)->Pre ()->Value ();
    else sprintf (errorMassage, "SettingDiffTimeIntegrationRightHand, pt->EWNS (%c, %c) does not exist, please check azimuth", i, j), PrintError (errorMassage);
}

void SettingDiffDiffTimeIntegrationRightHand (Point *pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2) {

  unordered_map <char, double> fe;
  unordered_map <char, double> fw;
  unordered_map <char, double> fn;
  unordered_map <char, double> fs;
  unordered_map <char, unordered_map <char, double>> f;
  char errorMassage[256];

  f['E'] = fe, f['W'] = fw, f['N'] = fn, f['S'] = fs;
  f['E']['E'] = 5.0E-1 * greens_integral_ttau (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
  f['W']['W'] = 5.0E-1 * greens_integral_ttau (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
  f['N']['N'] = 5.0E-1 * greens_integral_ttau (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
  f['S']['S'] = 5.0E-1 * greens_integral_ttau (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();

  if (pt->EWNS ('E', 'N') && pt->EWNS ('E', 'S'))
  f['E']['N'] = f['E']['E'] * (yb - pt->EWNS ('E', 'S')->Coord ().Value ('y')) / (pt->EWNS ('E', 'N')->Coord ().Value ('y') - pt->EWNS ('E', 'S')->Coord ().Value ('y')),
  f['E']['S'] = f['E']['E'] * (pt->EWNS ('E', 'N')->Coord ().Value ('y') - yb) / (pt->EWNS ('E', 'N')->Coord ().Value ('y') - pt->EWNS ('E', 'S')->Coord ().Value ('y'));

  if (pt->EWNS ('W', 'N') && pt->EWNS ('W', 'S'))
  f['W']['N'] = f['W']['W'] * (yb - pt->EWNS ('W', 'S')->Coord ().Value ('y')) / (pt->EWNS ('W', 'N')->Coord ().Value ('y') - pt->EWNS ('W', 'S')->Coord ().Value ('y')),
  f['W']['S'] = f['W']['W'] * (pt->EWNS ('W', 'N')->Coord ().Value ('y') - yb) / (pt->EWNS ('W', 'N')->Coord ().Value ('y') - pt->EWNS ('W', 'S')->Coord ().Value ('y'));

  if (pt->EWNS ('N', 'E') && pt->EWNS ('N', 'W'))
  f['N']['E'] = f['N']['N'] * (xb - pt->EWNS ('N', 'W')->Coord ().Value ('x')) / (pt->EWNS ('N', 'E')->Coord ().Value ('x') - pt->EWNS ('N', 'W')->Coord ().Value ('x')),
  f['N']['W'] = f['N']['N'] * (pt->EWNS ('N', 'E')->Coord ().Value ('x') - xb) / (pt->EWNS ('N', 'E')->Coord ().Value ('x') - pt->EWNS ('N', 'W')->Coord ().Value ('x'));

  if (pt->EWNS ('S', 'E') && pt->EWNS ('S', 'W'))
  f['S']['E'] = f['S']['S'] * (xb - pt->EWNS ('S', 'W')->Coord ().Value ('x')) / (pt->EWNS ('S', 'E')->Coord ().Value ('x') - pt->EWNS ('S', 'W')->Coord ().Value ('x')),
  f['S']['W'] = f['S']['S'] * (pt->EWNS ('S', 'E')->Coord ().Value ('x') - xb) / (pt->EWNS ('S', 'E')->Coord ().Value ('x') - pt->EWNS ('S', 'W')->Coord ().Value ('x'));

  xdat->F += 5.0E-1 * greens_integral_ttau (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt () * pt->Pre ()->Value ();
  xdat->F += 5.0E-1 * greens_integral_ttau (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt () * pt->Pre ()->Value ();

  ydat->F += 5.0E-1 * greens_integral_ttau (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt () * pt->Pre ()->Value ();
  ydat->F += 5.0E-1 * greens_integral_ttau (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt () * pt->Pre ()->Value ();

  for (const auto &i : {'E', 'W'})
  if (pt->EWNS (i, i)) xdat->F += f[i][i] * pt->EWNS (i, i)->Pre ()->Value ();
  else
  for (const auto &j : {'N', 'S'})
  if (pt->EWNS (i, j)) xdat->F += f[i][j] * pt->EWNS (i, j)->Pre ()->Value ();
  else sprintf (errorMassage, "SettingDiffDiffTimeIntegrationRightHand, pt->EWNS (%c, %c) does not exist, please check azimuth", i, j), PrintError (errorMassage);

  for (const auto &i : {'N', 'S'})
  if (pt->EWNS (i, i)) ydat->F += f[i][i] * pt->EWNS (i, i)->Pre ()->Value ();
  else
  for (const auto &j : {'E', 'W'})
  if (pt->EWNS (i, j)) ydat->F += f[i][j] * pt->EWNS (i, j)->Pre ()->Value ();
  else sprintf (errorMassage, "SettingDiffDiffTimeIntegrationRightHand, pt->EWNS (%c, %c) does not exist, please check azimuth", i, j), PrintError (errorMassage);
}

void SettingNeumannTimeIntegrationRightHand (Point *pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2,
  const bool is_sol) {

    unordered_map <char, double> fe;
    unordered_map <char, double> fw;
    unordered_map <char, double> fn;
    unordered_map <char, double> fs;
    unordered_map <char, unordered_map <char, double>> f;

    char errorMassage[256];

    unordered_map <char, char(*)[2]> s;
    unordered_map <char, double*> ff;

    f['E'] = fe, f['W'] = fw, f['N'] = fn, f['S'] = fs;
    f['E']['E'] = - 5.0E-1 * greens_integral_tau (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    f['W']['W'] = - 5.0E-1 * greens_integral_tau (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt ();
    f['N']['N'] = - 5.0E-1 * greens_integral_tau (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();
    f['S']['S'] = - 5.0E-1 * greens_integral_tau (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt ();

    if (pt->EWNS ('E', 'N') && pt->EWNS ('E', 'S'))
    f['E']['N'] = f['E']['E'] * (yb - pt->EWNS ('E', 'S')->Coord ().Value ('y')) / (pt->EWNS ('E', 'N')->Coord ().Value ('y') - pt->EWNS ('E', 'S')->Coord ().Value ('y')),
    f['E']['S'] = f['E']['E'] * (pt->EWNS ('E', 'N')->Coord ().Value ('y') - yb) / (pt->EWNS ('E', 'N')->Coord ().Value ('y') - pt->EWNS ('E', 'S')->Coord ().Value ('y'));

    if (pt->EWNS ('W', 'N') && pt->EWNS ('W', 'S'))
    f['W']['N'] = f['W']['W'] * (yb - pt->EWNS ('W', 'S')->Coord ().Value ('y')) / (pt->EWNS ('W', 'N')->Coord ().Value ('y') - pt->EWNS ('W', 'S')->Coord ().Value ('y')),
    f['W']['S'] = f['W']['W'] * (pt->EWNS ('W', 'N')->Coord ().Value ('y') - yb) / (pt->EWNS ('W', 'N')->Coord ().Value ('y') - pt->EWNS ('W', 'S')->Coord ().Value ('y'));

    if (pt->EWNS ('N', 'E') && pt->EWNS ('N', 'W'))
    f['N']['E'] = f['N']['N'] * (xb - pt->EWNS ('N', 'W')->Coord ().Value ('x')) / (pt->EWNS ('N', 'E')->Coord ().Value ('x') - pt->EWNS ('N', 'W')->Coord ().Value ('x')),
    f['N']['W'] = f['N']['N'] * (pt->EWNS ('N', 'E')->Coord ().Value ('x') - xb) / (pt->EWNS ('N', 'E')->Coord ().Value ('x') - pt->EWNS ('N', 'W')->Coord ().Value ('x'));

    if (pt->EWNS ('S', 'E') && pt->EWNS ('S', 'W'))
    f['S']['E'] = f['S']['S'] * (xb - pt->EWNS ('S', 'W')->Coord ().Value ('x')) / (pt->EWNS ('S', 'E')->Coord ().Value ('x') - pt->EWNS ('S', 'W')->Coord ().Value ('x')),
    f['S']['W'] = f['S']['S'] * (pt->EWNS ('S', 'E')->Coord ().Value ('x') - xb) / (pt->EWNS ('S', 'E')->Coord ().Value ('x') - pt->EWNS ('S', 'W')->Coord ().Value ('x'));

    if (!IsEqualDouble (pt->Normal ().Value ('x'), ZeroValue)) xdat->F -= 5.0E-1 * greens_integral_tau (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt () * pt->Pre ()->Value () * pt->Normal ().Value ('x');
    if (!IsEqualDouble (pt->Normal ().Value ('x'), ZeroValue)) xdat->F -= 5.0E-1 * greens_integral_tau (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / pt->Dt () * pt->Pre ()->Value () * pt->Normal ().Value ('x');

    if (!IsEqualDouble (pt->Normal ().Value ('y'), ZeroValue)) ydat->F -= 5.0E-1 * greens_integral_tau (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt () * pt->Pre ()->Value () * pt->Normal ().Value ('y');
    if (!IsEqualDouble (pt->Normal ().Value ('y'), ZeroValue)) ydat->F -= 5.0E-1 * greens_integral_tau (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / pt->Dt () * pt->Pre ()->Value () * pt->Normal ().Value ('y');

    for (const auto &i : {'E', 'W'})
    if (pt->EWNS (i, i)) {if (!IsEqualDouble (pt->Normal ().Value ('x'), ZeroValue)) xdat->F += f[i][i] * pt->EWNS (i, i)->Pre ()->Value () * pt->Normal ().Value ('x');}
    else
    for (const auto &j : {'N', 'S'})
    if (pt->EWNS (i, j)) {if (!IsEqualDouble (pt->Normal ().Value ('x'), ZeroValue)) xdat->F += f[i][j] * pt->EWNS (i, j)->Pre ()->Value () * pt->Normal ().Value ('x');}
    else sprintf (errorMassage, "SettingNeumannTimeIntegrationRightHand, pt->EWNS (%c, %c) does not exist, please check azimuth", i, j), PrintError (errorMassage);

    for (const auto &i : {'N', 'S'})
    if (pt->EWNS (i, i)) {if (!IsEqualDouble (pt->Normal ().Value ('y'), ZeroValue)) ydat->F += f[i][i] * pt->EWNS (i, i)->Pre ()->Value () * pt->Normal ().Value ('y');}
    else
    for (const auto &j : {'E', 'W'})
    if (pt->EWNS (i, j)) {if (!IsEqualDouble (pt->Normal ().Value ('y'), ZeroValue)) ydat->F += f[i][j] * pt->EWNS (i, j)->Pre ()->Value () * pt->Normal ().Value ('y');}
    else sprintf (errorMassage, "SettingNeumannTimeIntegrationRightHand, pt->EWNS (%c, %c) does not exist, please check azimuth", i, j), PrintError (errorMassage);
}

void SettingDirichletTimeIntegrationRightHand (Point *pt, Point *calc_pt, xData *xdat, yData *ydat,
  const double xm, const double xb, const double xp,
  const double ym, const double yb, const double yp,
  const int bdx, const int bdy,
  const double mpx1, const double mpx2,
  const double mpy1, const double mpy2,
  const bool is_sol) {
  // double  sign = is_sol ? 1.0E0 : -1.0E0;

  unordered_map <char, double> fe;
  unordered_map <char, double> fw;
  unordered_map <char, double> fn;
  unordered_map <char, double> fs;
  unordered_map <char, unordered_map <char, double>> f;

  f['E'] = fe, f['W'] = fw, f['N'] = fn, f['S'] = fs;
  f['E']['E'] = - 5.0E-1 *  greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / calc_pt->Dt ();
  f['W']['W'] = - 5.0E-1 *  greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / calc_pt->Dt ();
  f['N']['N'] = - 5.0E-1 *  greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / calc_pt->Dt ();
  f['S']['S'] = - 5.0E-1 *  greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / calc_pt->Dt ();

  xdat->F -= 5.0E-1 *  greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / calc_pt->Dt () * calc_pt->Pre ()->Value ();
  xdat->F -= 5.0E-1 *  greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) / calc_pt->Dt () * calc_pt->Pre ()->Value ();

  ydat->F -= 5.0E-1 *  greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / calc_pt->Dt () * calc_pt->Pre ()->Value ();
  ydat->F -= 5.0E-1 *  greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) / calc_pt->Dt () * calc_pt->Pre ()->Value ();

  unordered_map<char, char> opposite_azimuth;
  opposite_azimuth['E'] = 'W', opposite_azimuth['W'] = 'E', opposite_azimuth['N'] = 'S', opposite_azimuth['S'] = 'N';

  for (const auto &i : {'E', 'W'})
  if (calc_pt->EWNS (i, i))
  if (calc_pt->EWNS (i, i)->EWNS (i, i))
  if (calc_pt->EWNS (i, i)->EWNS (i, i) == pt) {
    xdat->F += f[i][i] * pt->Pre ()->Value ();
    if (calc_pt->EWNS (opposite_azimuth[i], opposite_azimuth[i])) xdat->F += f[opposite_azimuth[i]][opposite_azimuth[i]] * calc_pt->EWNS (opposite_azimuth[i], opposite_azimuth[i])->Pre ()->Value ();
    else PrintError ("SettingDirichletTimeIntegrationRightHand");
  }

  for (const auto &i : {'N', 'S'})
  if (calc_pt->EWNS (i, i))
  if (calc_pt->EWNS (i, i)->EWNS (i, i))
  if (calc_pt->EWNS (i, i)->EWNS (i, i) == pt) {
    ydat->F += f[i][i] * pt->Pre ()->Value ();
    if (calc_pt->EWNS (opposite_azimuth[i], opposite_azimuth[i])) ydat->F += f[opposite_azimuth[i]][opposite_azimuth[i]] * calc_pt->EWNS (opposite_azimuth[i], opposite_azimuth[i])->Pre ()->Value ();
    else PrintError ("SettingDirichletTimeIntegrationRightHand");
  }

  unordered_map<char, char> coord;
  coord['E'] = 'x', coord['W'] = 'x', coord['N'] = 'y', coord['S'] = 'y';

  unordered_map<char, string> axis;
  axis['x'] = "Y", axis['y'] = "X";
  for (const auto &i : {'E', 'W', 'N', 'S'})
  if (pt->Axis (coord[i]) > -1)
  if (pt->IsBoundary (opposite_azimuth[i]))
  if (pt->EWNS (i, i))
  if (pt->EWNS (i, i)->EWNS (i, i))
  if (pt->EWNS (i, i)->EWNS (i, i)->Condition () == 'C') {InitializationCoef (xdat, ydat, axis[coord[i]]); break;}
}

#endif