#ifndef UTIL_H
#define UTIL_H

#include "Header.hpp"
#include "Greenfunctions.hpp"
#include "ftns.hpp"

void ProgressBar (int presentValue, int maximumValue) {

  int barWidth = 70;
  double progress = float (presentValue) / float (maximumValue);
  int pos = barWidth * progress;

  if (presentValue < maximumValue) {
    printf ("[");
    for (size_t i = 0; i < barWidth; ++i) {
      if (i < pos) printf ("=");
      else if (i == pos) printf (">");
      else printf (".");
    }
    printf ("] %.2f%% (%d/%d)\r", progress * 1.0E2, presentValue, maximumValue);
  } else {
    printf ("[");
    for (size_t i = 0; i < barWidth; ++i) {
      if (i < barWidth) printf ("=");
      else if (i == barWidth) printf (">");
    }
    printf ("] %.2f%% (%d/%d)\r",1.0E2, presentValue, maximumValue);
    printf ("\n");
  }

}

void SettingSolver (string Solvertype) {
  /* Solver Type               */
  /* Elliptic solver     : 0   */
  /* Heat transfer solver: 1   */
  char errorMassage[256];

  for (auto &i : Solvertype) i = toupper (i);

  if      (!Solvertype.compare ("ELLIPTIC")) ::SolverType = 0;
  else if (!Solvertype.compare ("HEAT"))     ::SolverType = 1;
  else sprintf (errorMassage, "SettingSolver, Solvertype = \"%s\" is wrong", Solvertype.c_str ()), PrintError (errorMassage);

}

/* double형의 두 실수가 같은지 비교 (NearZero값으로 비교) */
bool IsEqualDouble (double v1, double v2) {

  // v1과 v2의 차이가 NearZero보다 작으면 true
    return fabs(v1 - v2) < NearZero;

  // 크면 false
}

/* double형의 두 실수가 같은지 비교 (tolerance가 주어진 비교) */
bool IsEqualDouble (double v1, double v2, double tol) {

  // v1과 v2의 차이가 주어진 tolerance 작으면 true
    return fabs(v1 - v2) < tol;

  // 크면 false
}

int Digits (int integer) {
  if (integer < 0) return 1;

  int digits = 0;
  while (integer) {
    integer /= 10;
    digits++;
  }

  return digits;
}

/* integer set인 memberSet의 0 ~ range안에 target이 있으면 있으면 위치를 return 없으면 -1을 return */
int is_member (int target, int *memberSet, int range) {

  // range가 0보다 작거나 같으면 -1을 return
  if (range <= 0) return -1;

  // 0부터 range까지 memberSet을 검색하면서 target과 같으면 그 위치를 return
  for (size_t i = 0; i < range; i++) if (target == memberSet[i]) return i;

  // memberSet안에 target과 같은 member가 없으면 -1을 return
  return -1;
}

/* double set인 memberSet의 0 ~ range안에 target이 있으면 있으면 위치를 return 없으면 -1을 return */
int is_member (double target, double *memberSet, int range) {

  // range가 0보다 작거나 같으면 -1을 return
  if (range <= 0) return -1;

  // 0부터 range까지 memberSet을 검색하면서 target과 같으면 그 위치를 return
  for (size_t i = 0; i < range; i++) if (IsEqualDouble (target, memberSet[i])) return i;

  // memberSet안에 target과 같은 member가 없으면 -1을 return
  return -1;
}

/* 두 점의 거리를 측정 */
double point_distance (double x1, double y1, double x2, double y2) {
  // 두 점 사이의 거리를 return
  return sqrt ((x1 - x2) *(x1 - x2) +(y1 - y2) *(y1 - y2));

}

/* 반올림을 하는 모듈 */
double round (double value, int pos) {
  // 임시로 갓을 저장하기 위한 변수
  double temp;
  // 원하는 소수점 자리수만큼 10의 누승을 함
  temp = value * pow (10, pos);
  // 0.5를 더한후 버림하면 반올림이 됨
  temp = floor (temp + 0.5);
  // 다시 원래 소수점 자리수로
  temp *= pow (10, -pos);
  // 값을 return
  return temp;
}


/* Interface위의 점 또는 neumann boundary condition을 갖는 경계점의 개수를 세는 모듈 */
int CountInterface (AxialData *adat, Point *pt) {
  // Interface위의 점의 개수를 세는 정수
  int j = 0;
  // 모든점에 대해서 반복
  for (size_t i = 0; i < adat->Pts_Num (); i++) if (pt[i].Condition () == 'I') j++;
  // 센 개수를 return
  return j;
}

/* 점들을 모든점과 내부점, phi의 값을 계산하는 점들로 분류하는 모듈 */
void sortPts (AxialData *adat, Point *pt) {
  // 점들을 검색할 때 이용하는 정수
  int j = 0;
  // 모든점에 대해서 반복
  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    // 내부점인 경우
    if (pt[i].Condition () == 'C') {
      // phi의 값을 계산하는 점의 index에서 모든점의 index
      adat->SetPtsTOpts ('H', j, i);
      // 모든점의 index에서 phi의 값을 계산하는 점의 index
      adat->SetPtsTOpts ('T', i, j);
      // 찾은점의 개수에 1을 더한다
      j += 1;
    }
  }
  // 모든점에 대해서 반복
  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    // 주위의 conductivity가 같은 Interface위의 점인 경우
    if (pt[i].Condition () == 'M' || pt[i].Condition () == 'N' || pt[i].Condition () == 'D') {
      // phi의 값을 계산하는 점의 index에서 모든점의 index
      adat->SetPtsTOpts ('H', j, i);
      // 모든점의 index에서 phi의 값을 계산하는 점의 index
      adat->SetPtsTOpts ('T', i, j);
      // 찾은점의 개수에 1을 더한다
      j += 1;
    }
  }
  // j를 다시 0으로 초기화 (처음부터 다시 셈)
  j = 0;
  // 모든점에 대해서 반복
  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    // 내부점인 경우
    if (pt[i].Condition () == 'C') {
      // 내부점의 index에서 모든점의 index
      adat->SetPtsTOpts ('I', j, i);
      // 모든점의 index에서 내부점이 index
      adat->SetPtsTOpts ('P', i, j);
      // 찾은점의 개수에 1을 더한다
      j += 1;
    }
  }
  // 모든점에 대해서 반복
  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    // Interface위의 점인 경우
    if (pt[i].Condition () == 'I' || pt[i].Condition () == 'M' || pt[i].Condition () == 'N' || pt[i].Condition () == 'D') {
      // 내부점의 index에서 모든점의 index
      adat->SetPtsTOpts ('I', j, i);
      // 모든점의 index에서 내부점이 index
      adat->SetPtsTOpts ('P', i, j);
      // 찾은점의 개수에 1을 더한다
      j += 1;
    }
  }
}

/* 4pt Gaussian quadrature rule integration */
double gauss_quadrature (std::function<double (double)> fp, double alpha, double beta) {
  /*
  fp: 적분을 할 함수 (input: double형, output: double형)
  alpha: 적분 시작점
  beta: 적분 끝점
  */
  // N: Gaussian quadrature rule에서 사용 할 point의 개수
  int N = 4;
  double x[N], w[N], c[2];
  double xx = 0.0;
  double gq = 0.0;

  // x: [-1, 1]에서의 N개의 점의 좌표
  x[0] =  sqrt (3.0/7.0 - 2.0/7.0 * sqrt (6.0/5.0));
  x[1] = -sqrt (3.0/7.0 - 2.0/7.0 * sqrt (6.0/5.0));
  x[2] =  sqrt (3.0/7.0 + 2.0/7.0 * sqrt (6.0/5.0));
  x[3] = -sqrt (3.0/7.0 + 2.0/7.0 * sqrt (6.0/5.0));

  // w: x에 대응하는 weight
  w[0] = (18.0 + sqrt (30.0))/36.0;
  w[1] = (18.0 + sqrt (30.0))/36.0;
  w[2] = (18.0 - sqrt (30.0))/36.0;
  w[3] = (18.0 - sqrt (30.0))/36.0;

  c[0] = (beta - alpha)/2.0;
  c[1] = (beta + alpha)/2.0;

  // [alpha, beta]에서의 N개의 함수의 값에 weight를 곱해서 더한다.
  for (size_t i = 0; i < N; i++) {
    // [alpha, beta]에서의 N개의 점의 좌표
    xx  = c[0] * x[i] + c[1];
    // 함수의 값에 weight를 곱해서 더한다.
    gq += w[i] * fp (xx);
  }
  return c[0] * gq;
}

/* 해의 relative L2-error를 계산해 주는 모듈 */
double Calc_Error (AxialData *adat, Point *pts) {
  // 현재점의 x-좌표와 y-좌표
  double xb, yb;
  // 현재점을 포함하는 국소 x-축선의 왼쪽 끝점과 오른쪽 끝점의 x-좌표
  double xm, xp;
  // 현재점을 포함하는 국소 y-축선의 아래쪽 끝점과 위쪽 끝점의 y-좌표
  double ym, yp;
  // exact solution의 L2-norm을 저장하기 위한 변수
  double error1 = 0.0;
  // error의 L2-norm을 저정하기 위한 변수
  double error2 = 0.0;
  // exact solution의 infinite-norm을 저장하기 위한 변수
  double norm1 = 0.0;
  // error의 infinite-norm을 저정하기 위한 변수
  double norm2 = 0.0;
  // error가 최대가 되는 점의 index를 찾기위한 변수
  int MaxErrorIndex = -1;
  // 모든점에 대해서 반복
  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    // 참조하는 점이 Dirichlet 경계점이거나 Neumann 경계점이거나 Infinity 점이거나 Singularity점인 경우 pass함
    if (HeadVelocity (&pts[i])->Condition () == 'F' || HeadVelocity (&pts[i])->Condition () == 'S') continue;
    // 참조하는 점의 x-좌표
    xb = HeadVelocity (&pts[i])->Coord ().Value ('x');
    // 참조하는 점의 y-좌표
    yb = HeadVelocity (&pts[i])->Coord ().Value ('y');
    // 현재점을 포함하는 국소 x-축선의 왼쪽 끝점의 x-좌표
    xm = HeadVelocity (&pts[i])->MinMaxCoordinate ('x', 'm');
    // 현재점을 포함하는 국소 x-축선의 오른쪽 끝점의 x-좌표
    xp = HeadVelocity (&pts[i])->MinMaxCoordinate ('x', 'p');
    // 현재점을 포함하는 국소 y-축선의 아래쪽 끝점의 y-좌표
    ym = HeadVelocity (&pts[i])->MinMaxCoordinate ('y', 'm');
    // 현재점을 포함하는 국소 y-축선의 위쪽 끝점의 y-좌표
    yp = HeadVelocity (&pts[i])->MinMaxCoordinate ('y', 'p');

    // 참조하는 점에서의 exact solution의 값이 저장한 exact solution의 최대값보다 더 큰 경우, 대체해서 저장
    if (fabs (u_ftn_Dirichlet (&pts[i])) > norm1) norm1 = fabs (u_ftn_Dirichlet (&pts[i]));
    // 참조하는 점에서의 error가 저장한 error의 최대값보다 더 큰 경우, 대체해서 저장
    if (fabs (pts[i].Value () - u_ftn_Dirichlet (&pts[i])) > norm2) {norm2 = fabs (pts[i].Value () - u_ftn_Dirichlet (&pts[i])); MaxErrorIndex = pts[i].Index ();}
    // 참조하는 점에서의 exact solution의 L2-norm을 저장한 L2-norm에 더함
    error1 += 0.25 * u_ftn_Dirichlet (&pts[i]) * u_ftn_Dirichlet (&pts[i]) * (xp - xm) * (yp - ym);
    // 참조하는 점에서의 error의 L2-norm을 저장한 L2-norm에 더함
    error2 += 0.25 *(pts[i].Value () - u_ftn_Dirichlet (&pts[i])) *(pts[i].Value () - u_ftn_Dirichlet (&pts[i])) * (xp - xm) * (yp - ym);
  }
  // 줄바꿈
  printf ("\n");
  printf ("Error calculation of %s\n", pts[0].Mark ().c_str ());
  // exact solution의 infinite-norm을 출력
  printf ("%-50s%23.16e\n", "Infinite-norm of the solution = ", norm1);
  // error의 infinite-norm을 출력
  printf ("%-50s%23.16e\n", "Infinite-norm of the error    = ", norm2);
  // relative infinite-norm을 출력
  printf ("%-50s%23.16e\n", "Relative Infinite-error       = ", norm2 / norm1);
  // exact solution의 L2-norm을 출력
  printf ("%-50s%23.16e\n", "L2-norm of the solution       = ", sqrt (error1));
  // error의 L2-norm을 출력
  printf ("%-50s%23.16e\n", "L2-norm of the error          = ", sqrt (error2));
  // relative L2-norm을 출력
  printf ("%-50s%23.16e\n", "Relative L2-error             = ", sqrt (error2) / sqrt (error1));
  // error가 최대가 되는 점의 index를 출력
  printf ("%-50s%d\n", "Maximum error occurs at ", MaxErrorIndex);
  // error가 최대가 되는 점의 index를 출력
  // printf ("%-50s%23.16e\n", "Undetermined constant beta in unbounded domain = ", Calc_Undetermined_Coefficient_in_Unbounded_Domian (adat, pts));
  // error의 L2-norm을 return
  return sqrt (error2) / sqrt (error1);
}

/* 미분의 relative L2-error를 계산해 주는 모듈 */
double Calc_Derivative_Error (AxialData *adat, Point *pts) {
  // 현재점의 x-좌표와 y-좌표
  double xb, yb;
  // 현재점을 포함하는 국소 x-축선의 왼쪽 끝점과 오른쪽 끝점의 x-좌표
  double xm, xp;
  // 현재점을 포함하는 국소 y-축선의 아래쪽 끝점과 위쪽 끝점의 y-좌표
  double ym, yp;
  // exact solution의 L2-norm을 저장하기 위한 변수
  double error1 = 0.0;
  // error의 L2-norm을 저정하기 위한 변수
  double error2 = 0.0;
  // exact solution의 infinite-norm을 저장하기 위한 변수
  double norm1 = 0.0;
  // error의 infinite-norm을 저정하기 위한 변수
  double norm2 = 0.0;
  // error가 최대가 되는 점의 index를 찾기위한 변수
  int MaxErrorIndex = -1;
  // exact 미분의 크기를 저장하기 위한 변수
  double Du_ftn = 0.0;
  // 미분의 크기를 저장하기 위한 변수
  double Du = 0.0;
  // 모든점에 대해서 반복
  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    // 참조하는 점이 Dirichlet 경계점이거나 Neumann 경계점이거나 Infinity 점이거나 Singularity점인 경우 pass함
    if (pts[i].Condition () == 'D' || pts[i].Condition () == 'N' || pts[i].Condition () == 'F' || pts[i].Condition () == 'S') continue;
    // 참조하는 점의 x-좌표
    xb = HeadVelocity (&pts[i])->Coord ().Value ('x');
    // 참조하는 점의 y-좌표
    yb = HeadVelocity (&pts[i])->Coord ().Value ('y');
    // 현재점을 포함하는 국소 x-축선의 왼쪽 끝점의 x-좌표
    xm = HeadVelocity (&pts[i])->MinMaxCoordinate ('x', 'm');
    // 현재점을 포함하는 국소 x-축선의 오른쪽 끝점의 x-좌표
    xp = HeadVelocity (&pts[i])->MinMaxCoordinate ('x', 'p');
    // 현재점을 포함하는 국소 y-축선의 아래쪽 끝점의 y-좌표
    ym = HeadVelocity (&pts[i])->MinMaxCoordinate ('y', 'm');
    // 현재점을 포함하는 국소 y-축선의 위쪽 끝점의 y-좌표
    yp = HeadVelocity (&pts[i])->MinMaxCoordinate ('y', 'p');
    // exact 미분의 크기를 저장
    Du_ftn = sqrt (dudx_ftn (xb, yb) * dudx_ftn (xb, yb) + dudy_ftn (xb, yb) * dudy_ftn (xb, yb));
    // 미분의 크기를 저장
    Du = sqrt (pts[i].Diff ('x')->Value () * pts[i].Diff ('x')->Value () + pts[i].Diff ('y')->Value () * pts[i].Diff ('y')->Value ());
    // 참조하는 점에서의 exact solution의 값이 저장한 exact solution의 최대값보다 더 큰 경우, 대체해서 저장
    if (fabs (Du_ftn) > norm1) norm1 = fabs (Du_ftn);
    // 참조하는 점에서의 error가 저장한 error의 최대값보다 더 큰 경우, 대체해서 저장
    if (fabs (Du - Du_ftn) > norm2) {norm2 = fabs (Du - Du_ftn); MaxErrorIndex = pts[i].Index ();}
    // 참조하는 점에서의 exact solution의 L2-norm을 저장한 L2-norm에 더함
    error1 += 0.25 * Du_ftn * Du_ftn *(xp - xm) *(yp - ym);
    // 참조하는 점에서의 error의 L2-norm을 저장한 L2-norm에 더함
    error2 += 0.25 *(Du - Du_ftn) *(Du - Du_ftn) *(xp - xm) *(yp - ym);
  }
  // 줄바꿈
  printf ("\n");
  // exact solution의 infinite-norm을 출력
  printf ("%-50s%23.16e\n", "Infinite-norm of the derivative = ", norm1);
  // error의 infinite-norm을 출력
  printf ("%-50s%23.16e\n", "Infinite-norm of the error      = ", norm2);
  // relative infinite-norm을 출력
  printf ("%-50s%23.16e\n", "Relative Infinite-error         = ", norm2 / norm1);
  // exact solution의 L2-norm을 출력
  printf ("%-50s%23.16e\n", "L2-norm of the derivative       = ", sqrt (error1));
  // error의 L2-norm을 출력
  printf ("%-50s%23.16e\n", "L2-norm of the error            = ", sqrt (error2));
  // relative L2-norm을 출력
  printf ("%-50s%23.16e\n", "Relative L2-error               = ", sqrt (error2) / sqrt (error1));
  // error가 최대가 되는 점의 index를 출력
  printf ("%-50s%d\n", "Maximum error occurs at ", MaxErrorIndex);
  // error의 L2-norm을 return
  return sqrt (error2) / sqrt (error1);
}

/* Unbounded domain에서 undetermined coefficient를 계산하는 모듈 */
double Calc_Undetermined_Coefficient_in_Unbounded_Domian (AxialData *adat, Point* pt) {
  // undetermined coefficient의 평균을 계산하기 위해서 개수를 세는 변수
  int N_undeterminedCoefficient = 0;
  // 2차원에서 점근적함수가 u (x, y) ~ U * y + c * 1인 경우의 변수 U
  // double U = 1.0;
  // undetermined coefficient의 평균을 계산하기 위한 변수
  double undeterminedCoefficient = 0.0;
  // 점근적함수에 input의 x-좌표로 사용하기 위한 변수
  double x = 0.0;
  // 점근적함수에 input의 y-좌표로 사용하기 위한 변수
  double y = 0.0;
  // 점근적 함수 u_infinity = u1_infiity + beta * u2_infitiy에서의 u1_infiity
  std::function<double(double, double)> u1_infinity = [&] (double x, double y) {return 4.0E0 * log (1.0E0 / sqrt (x * x + y * y));};
  // 점근적 함수 u_infinity = u1_infiity + beta * u2_infitiy에서의 u2_infiity
  std::function<double(double, double)> u2_infinity = [&] (double x, double y) {return 1.0E0;};
  // 모든점에 대해서 반복
  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    // 현재점이 경계점인 경우, pass
    if (pt[i].Condition () == 'D' || pt[i].Condition () == 'N' || pt[i].Condition () == 'F') continue;
    // 오른쪽점이 존재하고 오른쪽점의 경계조건이 Infinite조건인 경우,
    if (pt[i].EWNS ('E', 'E') != NULL) if (pt[i].EWNS ('E', 'E')->Condition () == 'F') {
      // 왼쪽점의 x-좌표를 저장
      x = pt[i].EWNS ('W', 'W')->Coord ().Value ('x');
      // 왼쪽점의 y-좌표를 저장
      y = pt[i].EWNS ('W', 'W')->Coord ().Value ('y');
      // undeterminedCoefficient에 축선의 왼쪽끝에서 undetermined coefficient를 계산해서 더한다
      undeterminedCoefficient += (pt[i].EWNS ('W', 'W')->Value () - u1_infinity (x, y)) / u2_infinity (x, y);
      // undetermined coefficient에 더한 개수에 1을 더한다
      N_undeterminedCoefficient += 1;
    }
    // 왼쪽점이 존재하고 왼쪽점의 경계조건이 Infinite조건인 경우,
    if (pt[i].EWNS ('W', 'W') != NULL) if (pt[i].EWNS ('W', 'W')->Condition () == 'F') {
      // 오른쪽점의 x-좌표를 저장
      x = pt[i].EWNS ('E', 'E')->Coord ().Value ('x');
      // 오른쪽점의 y-좌표를 저장
      y = pt[i].EWNS ('E', 'E')->Coord ().Value ('y');
      // undeterminedCoefficient에 축선의 오른쪽끝에서 undetermined coefficient를 계산해서 더한다
      undeterminedCoefficient += (pt[i].EWNS ('E', 'E')->Value () - u1_infinity (x, y)) / u2_infinity (x, y);
      // undetermined coefficient에 더한 개수에 1을 더한다
      N_undeterminedCoefficient += 1;
    }
    // 위쪽점이 존재하고 위쪽점의 경계조건이 Infinite조건인 경우,
    if (pt[i].EWNS ('N', 'N') != NULL) if (pt[i].EWNS ('N', 'N')->Condition () == 'F') {
      // 아래쪽점의 x-좌표를 저장
      x = pt[i].EWNS ('S', 'S')->Coord ().Value ('x');
      // 아래쪽점의 y-좌표를 저장
      y = pt[i].EWNS ('S', 'S')->Coord ().Value ('y');
      // undeterminedCoefficient에 축선의 아래쪽끝에서 undetermined coefficient를 계산해서 더한다
      undeterminedCoefficient += (pt[i].EWNS ('S', 'S')->Value () - u1_infinity (x, y)) / u2_infinity (x, y);
      // undetermined coefficient에 더한 개수에 1을 더한다
      N_undeterminedCoefficient += 1;
    }
    // 아래쪽점이 존재하고 아래쪽점의 경계조건이 Infinite조건인 경우,
    if (pt[i].EWNS ('S', 'S') != NULL) if (pt[i].EWNS ('S', 'S')->Condition () == 'F') {
      // 위쪽점의 x-좌표를 저장
      x = pt[i].EWNS ('N', 'N')->Coord ().Value ('x');
      // 위쪽점의 y-좌표를 저장
      y = pt[i].EWNS ('N', 'N')->Coord ().Value ('y');
      // undeterminedCoefficient에 축선의 위쪽끝에서 undetermined coefficient를 계산해서 더한다
      undeterminedCoefficient += (pt[i].EWNS ('N', 'N')->Value () - u1_infinity (x, y)) / u2_infinity (x, y);
      // undetermined coefficient에 더한 개수에 1을 더한다
      N_undeterminedCoefficient += 1;
    }

  }
  // undetermined coefficientd의 개수가 0보다 크다면 undeterminedCoefficient를 더한 개수로 나누어서 평균을 return
  if (N_undeterminedCoefficient > 0) return undeterminedCoefficient / N_undeterminedCoefficient;
  // undetermined coefficientd의 개수가 0보다 작다면 0을 return
  else                               return 0.0;
}

/* 에러메시지 출력 */
void PrintError (const char* massage) {
  // massage: 덧붙일 에러 메시지
  printf ("\n");
  for (size_t i = 0; i < 28; i++) printf ("%c", '='); printf ("%s", "Some error has occured"); for (size_t i = 0; i < 28; i++) printf ("%c", '='); printf ("\n");
  printf ("\t%s%s\n", "Error in ", massage);
  for (size_t i = 0; i < 78; i++) printf ("%c", '='); printf ("\n"); printf ("\n");

  exit (1);
}

/* 점의 정보를 출력 */
void PrintPtsInfo (Point* pt) {
  // 줄을 띄운다.
  printf ("\n");
  // 점의 index를 출력
  printf ("index = %d\n", pt->Index ());
  // 점의 경계조견을 출력
  printf ("Boundary condition = %c\n", pt->Condition ());
  // 오른쪽점이 존재하는 경우, 오른쪽점의 경계조건을 출력
  if (pt->EWNS ('E', 'E') != NULL) printf ("Boundary condition of E = %c\n", pt->EWNS ('E', 'E')->Condition ());
  // 왼쪽점이 존재하는 경우, 왼쪽점의 경계조건을 출력
  if (pt->EWNS ('W', 'W') != NULL) printf ("Boundary condition of W = %c\n", pt->EWNS ('W', 'W')->Condition ());
  // 위쪽점이 존재하는 경우, 위쪽점의 경계조건을 출력
  if (pt->EWNS ('N', 'N') != NULL) printf ("Boundary condition of N = %c\n", pt->EWNS ('N', 'N')->Condition ());
  // 아래쪽점이 존재하는 경우, 아래쪽점의 경계조건을 출력
  if (pt->EWNS ('S', 'S') != NULL) printf ("Boundary condition of S = %c\n", pt->EWNS ('S', 'S')->Condition ());
  // 오른쪽의 위쪽점이 존재하는 경우, 오른쪽의 위쪽점의 경계조건을 출력
  if (pt->EWNS ('E', 'N') != NULL) printf ("Boundary condition of EN = %c\n", pt->EWNS ('E', 'N')->Condition ());
  // 오른쪽의 아래쪽점이 존재하는 경우, 오른쪽의 아래쪽점의 경계조건을 출력
  if (pt->EWNS ('E', 'S') != NULL) printf ("Boundary condition of ES = %c\n", pt->EWNS ('E', 'S')->Condition ());
  // 왼쪽의 위쪽점이 존재하는 경우, 왼쪽의 위쪽점의 경계조건을 출력
  if (pt->EWNS ('W', 'N') != NULL) printf ("Boundary condition of WN = %c\n", pt->EWNS ('W', 'N')->Condition ());
  // 왼쪽의 아래쪽점이 존재하는 경우, 왼쪽의 아래쪽점의 경계조건을 출력
  if (pt->EWNS ('W', 'S') != NULL) printf ("Boundary condition of WS = %c\n", pt->EWNS ('W', 'S')->Condition ());
  // 위쪽의 오른쪽점이 존재하는 경우, 위쪽의 오른쪽점의 경계조건을 출력
  if (pt->EWNS ('N', 'E') != NULL) printf ("Boundary condition of NE = %c\n", pt->EWNS ('N', 'E')->Condition ());
  // 위쪽의 왼쪽점이 존재하는 경우, 위쪽의 왼쪽점의 경계조건을 출력
  if (pt->EWNS ('N', 'W') != NULL) printf ("Boundary condition of NW = %c\n", pt->EWNS ('N', 'W')->Condition ());
  // 아래쪽의 오른쪽점이 존재하는 경우, 아래의 오른쪽쪽점의 경계조건을 출력
  if (pt->EWNS ('S', 'E') != NULL) printf ("Boundary condition of SE = %c\n", pt->EWNS ('S', 'E')->Condition ());
  // 아래쪽의 왼쪽점이 존재하는 경우, 아래쪽의 왼쪽점의 경계조건을 출력
  if (pt->EWNS ('S', 'W') != NULL) printf ("Boundary condition of SW = %c\n", pt->EWNS ('S', 'W')->Condition ());

  printf ("\n\n[");
  // 오른쪽점이 존재하는 경우, 오른쪽점의 index를 출력
  if (pt->EWNS ('E', 'E') != NULL) printf ("%d, ", pt->EWNS ('E', 'E')->Index ());
  // 왼쪽점이 존재하는 경우, 왼쪽점의 index를 출력
  if (pt->EWNS ('W', 'W') != NULL) printf ("%d, ", pt->EWNS ('W', 'W')->Index ());
  // 위쪽점이 존재하는 경우, 위쪽점의 index를 출력
  if (pt->EWNS ('N', 'N') != NULL) printf ("%d, ", pt->EWNS ('N', 'N')->Index ());
  // 아래쪽점이 존재하는 경우, 아래쪽점의 index를 출력
  if (pt->EWNS ('S', 'S') != NULL) printf ("%d, ", pt->EWNS ('S', 'S')->Index ());
  // 오른쪽의 위쪽점이 존재하는 경우, 오른쪽의 위쪽점의 index를 출력
  if (pt->EWNS ('E', 'N') != NULL) printf ("%d, ", pt->EWNS ('E', 'N')->Index ());
  // 오른쪽의 아래쪽점이 존재하는 경우, 오른쪽의 아래쪽점의 index를 출력
  if (pt->EWNS ('E', 'S') != NULL) printf ("%d, ", pt->EWNS ('E', 'S')->Index ());
  // 왼쪽의 위쪽점이 존재하는 경우, 왼쪽의 위쪽점의 index를 출력
  if (pt->EWNS ('W', 'N') != NULL) printf ("%d, ", pt->EWNS ('W', 'N')->Index ());
  // 왼쪽의 아래쪽점이 존재하는 경우, 왼쪽의 아래쪽점의 index를 출력
  if (pt->EWNS ('W', 'S') != NULL) printf ("%d, ", pt->EWNS ('W', 'S')->Index ());
  // 위쪽의 오른쪽점이 존재하는 경우, 위쪽의 오른쪽점의 index를 출력
  if (pt->EWNS ('N', 'E') != NULL) printf ("%d, ", pt->EWNS ('N', 'E')->Index ());
  // 위쪽의 왼쪽점이 존재하는 경우, 위쪽의 왼쪽점의 index를 출력
  if (pt->EWNS ('N', 'W') != NULL) printf ("%d, ", pt->EWNS ('N', 'W')->Index ());
  // 아래쪽의 오른쪽점이 존재하는 경우, 아래의 오른쪽쪽점의 index를 출력
  if (pt->EWNS ('S', 'E') != NULL) printf ("%d, ", pt->EWNS ('S', 'E')->Index ());
  // 아래쪽의 왼쪽점이 존재하는 경우, 아래쪽의 왼쪽점의 index를 출력
  if (pt->EWNS ('S', 'W') != NULL) printf ("%d, ", pt->EWNS ('S', 'W')->Index ());

  printf ("]\n\n");

  // 점의 좌표를 출력
  printf ("XB = %23.16e;\n", pt->Coord ().Value ('x'));
  printf ("YB = %23.16e;\n", pt->Coord ().Value ('y'));
  printf ("XM = %23.16e;\n", pt->MinMaxCoordinate ('x', 'm'));
  printf ("YM = %23.16e;\n", pt->MinMaxCoordinate ('y', 'm'));
  printf ("XP = %23.16e;\n", pt->MinMaxCoordinate ('x', 'p'));
  printf ("YP = %23.16e;\n", pt->MinMaxCoordinate ('y', 'p'));
  // 점의 좌표를 출력
  printf ("(xb , yb) = (%23.16e, %23.16e)\n", pt->Coord ().Value ('x'), pt->Coord ().Value ('y'));
  // (xm, ym)을 출력
  printf ("(xm , ym) = (%23.16e, %23.16e)\n", pt->MinMaxCoordinate ('x', 'm'), pt->MinMaxCoordinate ('y', 'm'));
  // (xp, yp)을 출력
  printf ("(xp , yp) = (%23.16e, %23.16e)\n", pt->MinMaxCoordinate ('x', 'p'), pt->MinMaxCoordinate ('y', 'p'));
  // 줄을 띄운다.
  printf ("\n");
}

/* xData의 정보를 출력 */
void PrintXdata (xData* xdat) {
  // 줄을 띄운다.
  printf ("\n");
  // 각 정보를 출력
  printf ("   Cu = %23.16e\n", xdat->Cu);
  printf ("   Eu = %23.16e\n", xdat->Eu);
  printf ("   Wu = %23.16e\n", xdat->Wu);
  printf ("  ENu = %23.16e\n", xdat->ENu);
  printf ("  ESu = %23.16e\n", xdat->ESu);
  printf ("  WNu = %23.16e\n", xdat->WNu);
  printf ("  WSu = %23.16e\n", xdat->WSu);
  printf (" Cphi = %23.16e\n", xdat->Cphi);
  printf (" Ephi = %23.16e\n", xdat->Ephi);
  printf (" Wphi = %23.16e\n", xdat->Wphi);
  printf ("ENphi = %23.16e\n", xdat->ENphi);
  printf ("ESphi = %23.16e\n", xdat->ESphi);
  printf ("WNphi = %23.16e\n", xdat->WNphi);
  printf ("WSphi = %23.16e\n", xdat->WSphi);
  printf ("   Cf = %23.16e\n", xdat->Cf);
  printf ("   Ef = %23.16e\n", xdat->Ef);
  printf ("   Wf = %23.16e\n", xdat->Wf);
  printf ("  ENf = %23.16e\n", xdat->ENf);
  printf ("  ESf = %23.16e\n", xdat->ESf);
  printf ("  WNf = %23.16e\n", xdat->WNf);
  printf ("  WSf = %23.16e\n", xdat->WSf);
  printf ("    F = %23.16e\n", xdat->F);
  // 줄을 띄운다.
  printf ("\n");
}

/* xData의 정보를 출력 */
void PrintYdata (yData* ydat) {
  // 줄을 띄운다.
  printf ("\n");
  // 각 정보를 출력
  printf ("   Cu = %23.16e\n", ydat->Cu);
  printf ("   Nu = %23.16e\n", ydat->Nu);
  printf ("   Su = %23.16e\n", ydat->Su);
  printf ("  NEu = %23.16e\n", ydat->NEu);
  printf ("  NWu = %23.16e\n", ydat->NWu);
  printf ("  SEu = %23.16e\n", ydat->SEu);
  printf ("  SWu = %23.16e\n", ydat->SWu);
  printf (" Cphi = %23.16e\n", ydat->Cphi);
  printf (" Nphi = %23.16e\n", ydat->Nphi);
  printf (" Sphi = %23.16e\n", ydat->Sphi);
  printf ("NEphi = %23.16e\n", ydat->NEphi);
  printf ("NWphi = %23.16e\n", ydat->NWphi);
  printf ("SEphi = %23.16e\n", ydat->SEphi);
  printf ("SWphi = %23.16e\n", ydat->SWphi);
  printf ("   Cf = %23.16e\n", ydat->Cf);
  printf ("   Nf = %23.16e\n", ydat->Nf);
  printf ("   Sf = %23.16e\n", ydat->Sf);
  printf ("  NEf = %23.16e\n", ydat->NEf);
  printf ("  NWf = %23.16e\n", ydat->NWf);
  printf ("  SEf = %23.16e\n", ydat->SEf);
  printf ("  SWf = %23.16e\n", ydat->SWf);
  printf ("    F = %23.16e\n", ydat->F);
  // 줄을 띄운다.
  printf ("\n");
}

// 점의 오른쪽, 왼쪽, 위쪽, 아래쪽을 대입하는 모듈
void SetAllEWNS (AxialData *adat, Point *src, Point *dest, Point *dests) {
  // 오른쪽점이 존재하는 경우, 오른쪽점을 대입한다.
  if (src->EWNS ('E', 'E') != NULL)  dest->SetEWNS ('E', 'E', &dests[src->EWNS ('E', 'E')->Index ()]);
  // 왼쪽점이 존재하는 경우, 왼쪽점을 대입한다.
  if (src->EWNS ('W', 'W') != NULL)  dest->SetEWNS ('W', 'W', &dests[src->EWNS ('W', 'W')->Index ()]);
  // 위쪽점이 존재하는 경우, 위쪽점을 대입한다.
  if (src->EWNS ('N', 'N') != NULL)  dest->SetEWNS ('N', 'N', &dests[src->EWNS ('N', 'N')->Index ()]);
  // 아래쪽점이 존재하는 경우, 아래쪽점을 대입한다.
  if (src->EWNS ('S', 'S') != NULL)  dest->SetEWNS ('S', 'S', &dests[src->EWNS ('S', 'S')->Index ()]);
  // 오른쪽의 위쪽점이 존재하는 경우, 오른쪽의 위쪽점을 대입한다
  if (src->EWNS ('E', 'N') != NULL)  dest->SetEWNS ('E', 'N', &dests[src->EWNS ('E', 'N')->Index ()]);
  // 오른쪽의 아래쪽점이 존재하는 경우, 오른쪽의 아래쪽점을 대입한다
  if (src->EWNS ('E', 'S') != NULL)  dest->SetEWNS ('E', 'S', &dests[src->EWNS ('E', 'S')->Index ()]);
  // 왼쪽의 위쪽점이 존재하는 경우, 왼쪽의 위쪽점을 대입한다
  if (src->EWNS ('W', 'N') != NULL)  dest->SetEWNS ('W', 'N', &dests[src->EWNS ('W', 'N')->Index ()]);
  // 왼쪽의 아래쪽점이 존재하는 경우, 왼쪽의 아래쪽점을 대입한다
  if (src->EWNS ('W', 'S') != NULL)  dest->SetEWNS ('W', 'S', &dests[src->EWNS ('W', 'S')->Index ()]);
  // 위쪽의 오른쪽점이 존재하는 경우, 오른쪽의 위쪽점을 대입한다
  if (src->EWNS ('N', 'E') != NULL)  dest->SetEWNS ('N', 'E', &dests[src->EWNS ('N', 'E')->Index ()]);
  // 위쪽의 왼쪽점이 존재하는 경우, 위쪽의 왼쪽점을 대입한다
  if (src->EWNS ('N', 'W') != NULL)  dest->SetEWNS ('N', 'W', &dests[src->EWNS ('N', 'W')->Index ()]);
  // 아래쪽의 오른쪽점이 존재하는 경우, 아래쪽의 오른쪽점을 대입한다
  if (src->EWNS ('S', 'E') != NULL)  dest->SetEWNS ('S', 'E', &dests[src->EWNS ('S', 'E')->Index ()]);
  // 아래쪽의 왼쪽점이 존재하는 경우, 아래쪽의 왼쪽점을 대입한다
  if (src->EWNS ('S', 'W') != NULL)  dest->SetEWNS ('S', 'W', &dests[src->EWNS ('S', 'W')->Index ()]);
}

// 점의 두번째 오른쪽, 왼쪽, 위쪽, 아래쪽을 대입하는 모듈
void SetAllEWNS2nd (AxialData *adat, Point *src, Point *dest, Point *dests) {
  // 오른쪽점이 존재하는 경우, 오른쪽점을 대입한다.
  if (src->EWNS2nd ('E', 'E') != NULL)  dest->SetEWNS2nd ('E', 'E', &dests[src->EWNS2nd ('E', 'E')->Index ()]);
  // 왼쪽점이 존재하는 경우, 왼쪽점을 대입한다.
  if (src->EWNS2nd ('W', 'W') != NULL)  dest->SetEWNS2nd ('W', 'W', &dests[src->EWNS2nd ('W', 'W')->Index ()]);
  // 위쪽점이 존재하는 경우, 위쪽점을 대입한다.
  if (src->EWNS2nd ('N', 'N') != NULL)  dest->SetEWNS2nd ('N', 'N', &dests[src->EWNS2nd ('N', 'N')->Index ()]);
  // 아래쪽점이 존재하는 경우, 아래쪽점을 대입한다.
  if (src->EWNS2nd ('S', 'S') != NULL)  dest->SetEWNS2nd ('S', 'S', &dests[src->EWNS2nd ('S', 'S')->Index ()]);
  // 오른쪽의 위쪽점이 존재하는 경우, 오른쪽의 위쪽점을 대입한다
  if (src->EWNS2nd ('E', 'N') != NULL)  dest->SetEWNS2nd ('E', 'N', &dests[src->EWNS2nd ('E', 'N')->Index ()]);
  // 오른쪽의 아래쪽점이 존재하는 경우, 오른쪽의 아래쪽점을 대입한다
  if (src->EWNS2nd ('E', 'S') != NULL)  dest->SetEWNS2nd ('E', 'S', &dests[src->EWNS2nd ('E', 'S')->Index ()]);
  // 왼쪽의 위쪽점이 존재하는 경우, 왼쪽의 위쪽점을 대입한다
  if (src->EWNS2nd ('W', 'N') != NULL)  dest->SetEWNS2nd ('W', 'N', &dests[src->EWNS2nd ('W', 'N')->Index ()]);
  // 왼쪽의 아래쪽점이 존재하는 경우, 왼쪽의 아래쪽점을 대입한다
  if (src->EWNS2nd ('W', 'S') != NULL)  dest->SetEWNS2nd ('W', 'S', &dests[src->EWNS2nd ('W', 'S')->Index ()]);
  // 위쪽의 오른쪽점이 존재하는 경우, 오른쪽의 위쪽점을 대입한다
  if (src->EWNS2nd ('N', 'E') != NULL)  dest->SetEWNS2nd ('N', 'E', &dests[src->EWNS2nd ('N', 'E')->Index ()]);
  // 위쪽의 왼쪽점이 존재하는 경우, 위쪽의 왼쪽점을 대입한다
  if (src->EWNS2nd ('N', 'W') != NULL)  dest->SetEWNS2nd ('N', 'W', &dests[src->EWNS2nd ('N', 'W')->Index ()]);
  // 아래쪽의 오른쪽점이 존재하는 경우, 아래쪽의 오른쪽점을 대입한다
  if (src->EWNS2nd ('S', 'E') != NULL)  dest->SetEWNS2nd ('S', 'E', &dests[src->EWNS2nd ('S', 'E')->Index ()]);
  // 아래쪽의 왼쪽점이 존재하는 경우, 아래쪽의 왼쪽점을 대입한다
  if (src->EWNS2nd ('S', 'W') != NULL)  dest->SetEWNS2nd ('S', 'W', &dests[src->EWNS2nd ('S', 'W')->Index ()]);
}

Point * HeadVelocity (Point *pt) {
  char errorMassage[256];

  if (pt->Mark ().length () == 1) return pt;
  else if (pt->Mark ().length () < 4) {
    if      (pt->Mark ()[0] == 'u') return pt->Velocity ('u');
    else if (pt->Mark ()[0] == 'v') return pt->Velocity ('v');
    else if (!pt->Mark ().compare ("phi")) return pt->Velocity ('u');
    else if (!pt->Mark ().compare ("psi")) return pt->Velocity ('v');
  }
  else if (pt->Mark ().length () == 4) {
    if      (pt->Mark ()[3] == 'u') return pt->Velocity ('u');
    else if (pt->Mark ()[3] == 'v') return pt->Velocity ('v');
  }
  else if (pt->Mark ().length () == 6) {
    if      (pt->Mark ()[5] == 'u') return pt->Velocity ('u');
    else if (pt->Mark ()[5] == 'v') return pt->Velocity ('v');
  }
  else if (pt->Mark ().length () == 11) {
    if      (pt->Mark ()[10] == 'u') return pt->Velocity ('u');
    else if (pt->Mark ()[10] == 'v') return pt->Velocity ('v');
  }

  sprintf (errorMassage, "HeadVelocity, pt->Mark () == %s, please check point's mark", pt->Mark ().c_str ());
  PrintError (errorMassage);
  exit (1);
}

void SettingPoints (AxialData *adat, Point *ptU, Point *ptV, Point *ptP,
  Point *pt_diff_u_x, Point *pt_diff_u_y, Point *pt_diff_v_x, Point *pt_diff_v_y,
  Point *pt_diff_u_xx, Point *pt_diff_u_yy, Point *pt_diff_v_xx, Point *pt_diff_v_yy,
  Point *pt_hat_u, Point *pt_hat_v, Point *pt_pre_u, Point* pt_pre_v, Point *pt_old_u, Point *pt_old_v,
  Point *pt_phi, Point *pt_psi) {

  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    // u의 값을 저장하는 점의 표시
    ptU[i].SetMark ("u");
    // v의 값을 저장하는 점의 표시
    ptV[i].SetMark ("v");
    // p의 값을 저장하는 점의 표시
    ptP[i].SetMark ("p");
    // u의 x성분으로의 미분값을 저장하는 점의 표시
    pt_diff_u_x[i].SetMark ("ux");
    // u의 y성분으로의 미분값을 저장하는 점의 표시
    pt_diff_u_y[i].SetMark ("uy");
    // v의 x성분으로의 미분값을 저장하는 점의 표시
    pt_diff_v_x[i].SetMark ("vx");
    // v의 y성분으로의 미분값을 저장하는 점의 표시
    pt_diff_v_y[i].SetMark ("vy");
    // u의 x성분으로의 미분의 x성분으로의 미분값을 저장하는 점의 표시
    pt_diff_u_xx[i].SetMark ("uxx");
    // u의 y성분으로의 미분의 y성분으로의 미분값을 저장하는 점의 표시
    pt_diff_u_yy[i].SetMark ("uyy");
    // v의 x성분으로의 미분의 x성분으로의 미분값을 저장하는 점의 표시
    pt_diff_v_xx[i].SetMark ("vxx");
    // v의 y성분으로의 미분의 y성분으로의 미분값을 저장하는 점의 표시
    pt_diff_v_yy[i].SetMark ("vyy");
    // 중간 속도 u를 저장하는 점의 표시
    pt_hat_u[i].SetMark ("hatu");
    // 중간 속도 v를 저장하는 점의 표시
    pt_hat_v[i].SetMark ("hatv");
    // 이전의 속도 u를 저장하는 점의 표시
    pt_pre_u[i].SetMark ("preu");
    // 이전의 속도 v를 저장하는 점의 표시
    pt_pre_v[i].SetMark ("prev");
    // 전전의 속도 u를 저장하는 점의 표시
    pt_old_u[i].SetMark ("oldu");
    // 전전의 속도 v를 저장하는 점의 표시
    pt_old_v[i].SetMark ("oldv");
    // phi를 저장하는 점의 표시
    pt_phi[i].SetMark ("phi");
    // psi를 저장하는 점의 표시
    pt_psi[i].SetMark ("psi");

    // velocityd에 자신의 주소를 저장
    ptU[i].SetVelocity ('u', &ptU[i]);
    ptV[i].SetVelocity ('v', &ptV[i]);
    ptP[i].SetVelocity ('u', &ptU[i]);
    ptP[i].SetVelocity ('v', &ptV[i]);
    pt_diff_u_x[i].SetVelocity ('u', &ptU[i]);
    pt_diff_u_y[i].SetVelocity ('u', &ptU[i]);
    pt_diff_u_xx[i].SetVelocity ('u', &ptU[i]);
    pt_diff_u_yy[i].SetVelocity ('u', &ptU[i]);
    pt_diff_v_x[i].SetVelocity ('v', &ptV[i]);
    pt_diff_v_y[i].SetVelocity ('v', &ptV[i]);
    pt_diff_v_xx[i].SetVelocity ('v', &ptV[i]);
    pt_diff_v_yy[i].SetVelocity ('v', &ptV[i]);
    pt_hat_u[i].SetVelocity ('u', &ptU[i]);
    pt_hat_v[i].SetVelocity ('v', &ptV[i]);
    pt_pre_u[i].SetVelocity ('u', &ptU[i]);
    pt_pre_v[i].SetVelocity ('v', &ptV[i]);
    pt_old_u[i].SetVelocity ('u', &ptU[i]);
    pt_old_v[i].SetVelocity ('v', &ptV[i]);
    pt_phi[i].SetVelocity ('u', &ptU[i]);
    pt_psi[i].SetVelocity ('v', &ptV[i]);

    // opposite velcotiy를 저장
    ptU[i].SetVelocity ('v', &ptV[i]);
    ptV[i].SetVelocity ('u', &ptU[i]);
    pt_diff_u_x[i].SetVelocity ('v', &ptV[i]);
    pt_diff_u_y[i].SetVelocity ('v', &ptV[i]);
    pt_diff_u_xx[i].SetVelocity ('v', &ptV[i]);
    pt_diff_u_yy[i].SetVelocity ('v', &ptV[i]);
    pt_diff_v_x[i].SetVelocity ('u', &ptU[i]);
    pt_diff_v_y[i].SetVelocity ('u', &ptU[i]);
    pt_diff_v_xx[i].SetVelocity ('u', &ptU[i]);
    pt_diff_v_yy[i].SetVelocity ('u', &ptU[i]);
    pt_hat_u[i].SetVelocity ('v', &ptV[i]);
    pt_hat_v[i].SetVelocity ('u', &ptU[i]);
    pt_pre_u[i].SetVelocity ('v', &ptV[i]);
    pt_pre_v[i].SetVelocity ('u', &ptU[i]);
    pt_old_u[i].SetVelocity ('v', &ptV[i]);
    pt_old_v[i].SetVelocity ('u', &ptU[i]);
    pt_phi[i].SetVelocity ('v', &ptU[i]);
    pt_psi[i].SetVelocity ('u', &ptV[i]);

    // 압력을 저장
    ptU[i].SetPressure (&ptP[i]);
    ptV[i].SetPressure (&ptP[i]);
    ptP[i].SetPressure (&ptP[i]);
    pt_diff_u_x[i].SetPressure (&ptP[i]);
    pt_diff_u_y[i].SetPressure (&ptP[i]);
    pt_diff_u_xx[i].SetPressure (&ptP[i]);
    pt_diff_u_yy[i].SetPressure (&ptP[i]);
    pt_diff_v_x[i].SetPressure (&ptP[i]);
    pt_diff_v_y[i].SetPressure (&ptP[i]);
    pt_diff_v_xx[i].SetPressure (&ptP[i]);
    pt_diff_v_yy[i].SetPressure (&ptP[i]);
    pt_hat_u[i].SetPressure (&ptP[i]);
    pt_hat_v[i].SetPressure (&ptP[i]);
    pt_pre_u[i].SetPressure (&ptP[i]);
    pt_pre_v[i].SetPressure (&ptP[i]);
    pt_old_u[i].SetPressure (&ptP[i]);
    pt_old_v[i].SetPressure (&ptP[i]);
    pt_phi[i].SetPressure (&ptP[i]);
    pt_psi[i].SetPressure (&ptP[i]);

    // 미분을 저장
    ptU[i].SetDiff ('x', &pt_diff_u_x[i]);
    ptU[i].SetDiff ('y', &pt_diff_u_y[i]);
    ptV[i].SetDiff ('x', &pt_diff_v_x[i]);
    ptV[i].SetDiff ('y', &pt_diff_v_y[i]);

    pt_diff_u_x[i].SetDiff ('x', &pt_diff_u_xx[i]);
    pt_diff_u_y[i].SetDiff ('y', &pt_diff_u_yy[i]);
    pt_diff_v_x[i].SetDiff ('x', &pt_diff_v_xx[i]);
    pt_diff_v_y[i].SetDiff ('y', &pt_diff_v_yy[i]);

    // 미분하기 전의 점을 저장
    pt_diff_u_x[i].SetIntg (&ptU[i]);
    pt_diff_u_y[i].SetIntg (&ptU[i]);
    pt_diff_v_x[i].SetIntg (&ptV[i]);
    pt_diff_v_y[i].SetIntg (&ptV[i]);

    pt_diff_u_xx[i].SetIntg (&pt_diff_u_x[i]);
    pt_diff_u_yy[i].SetIntg (&pt_diff_u_y[i]);
    pt_diff_v_xx[i].SetIntg (&pt_diff_v_x[i]);
    pt_diff_v_yy[i].SetIntg (&pt_diff_v_y[i]);

    // 중간의 속도의 점을 저장
    ptU[i].SetHat (&pt_hat_u[i]);
    ptV[i].SetHat (&pt_hat_v[i]);

    // 이전의 속도의 점을 저장
    ptU[i].SetPre (&pt_pre_u[i]);
    ptV[i].SetPre (&pt_pre_v[i]);

    // 전전의 속도의 점을 저장
    ptU[i].SetOld (&pt_old_u[i]);
    ptV[i].SetOld (&pt_old_v[i]);

    // phi의 점을 저장
    ptU[i].SetPhi (&pt_phi[i]);
    ptV[i].SetPhi (&pt_psi[i]);
  }

}

#endif
