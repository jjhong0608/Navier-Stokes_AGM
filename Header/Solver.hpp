#ifndef ELLIPTICSOLVER_HPP
#define ELLIPTICSOLVER_HPP

#include "CalcDiff.hpp"

void Solver (ControlData *cdat, AxialData *adat, Point *pt) {
  // 행렬의 크기를 측정함을 알림
  printf ("\n%s\n", "Measure the size of a matrix:");
  // x-축선위에서의 해의 표현식의 계수들의 struct 선언
  xData xdat;

  // y-축선위에서의 해의 표현식의 계수들의 struct 선언
  yData ydat;

  // 행렬의 정보를 가지는 class의 선언
  MatrixProcess matps (adat);

  // 행렬을 만드는데 포함되는 점을 세기위한 정수의 선언
  int j = 0;

  //
  bool is_sol = true;

  bool ND = false;

  // 행렬의 크기가 2n x 2n이므로 모든점에 대해서 두번씩 반복
  for (size_t i = 0; i < 2 * adat->Pts_Num (); i++) {

    // 양쪽의 conductivity가 다른 interface위의 점인 경우 넘어간다
    if (i >= adat->Pts_Num ()) if (pt[i % adat->Pts_Num ()] == 'I') continue;

    // 내부점이거나, interface위의 점 또는 neumann 경계점인 경우
    if (pt[i % adat->Pts_Num ()] == 'C' || pt[i % adat->Pts_Num ()] == 'I' || pt[i % adat->Pts_Num ()] == 'M' || pt[i % adat->Pts_Num ()] == 'N' || pt[i % adat->Pts_Num ()] == 'D') {

      is_sol = i < adat->Pts_Num();

      // Neumann condition이 주어진 점의 두번째 식을 Dirichlet condition이 주어진 점과 같은 방식으로 계산
      if (pt[i % adat->Pts_Num ()] == 'N' && !is_sol) ND = true, pt[i % adat->Pts_Num ()].SetCondition ('D');

      // 내부점이거나, Interface인 경우에 해의 표현식의 계수들을 계산
      CalcRepresenCoef (&pt[i % adat->Pts_Num ()], &xdat, &ydat, is_sol);

      // 경계 근처에서의 해의 표현식을 처리
      TransposeBoundaryData (&pt[i % adat->Pts_Num ()], &xdat, &ydat, is_sol);

      // 행렬의 0이 아닌 원소의 개수를 셈
      matps.countEnt_num (j, adat, &pt[i % adat->Pts_Num ()], &xdat, &ydat);

      if (ND) pt[i % adat->Pts_Num ()].SetCondition ('N'), ND = false;

      // 행렬을 만드는데 포함되는 점에 1을 더함
      j += 1;

    }

    // 모든점에 대한 현재 진행중인 점의 index를 알림
    ProgressBar (i + 1, 2 * adat->Pts_Num ());
  }

  // 행렬의 크기를 측정했음을 알림
  printf ("\n%s\n", "Matrix size measurement done");

  // 행렬을 구성함을 알림
  printf ("\n%s\n", "Make the matrix:");

  // 행렬 정보의 초기화
  matps.initialization (adat, j);

  // 행렬을 만드는데 포함되는 점을 세기위한 정수의 초기화
  j = 0;

  // 행렬의 크기가 2n x 2n이므로 모든점에 대해서 두번씩 반복
  for (size_t i = 0; i < 2 * adat->Pts_Num (); i++) {

    // 양쪽의 conductivity가 다른 interface위의 점인 경우 넘어간다
    if (i >= adat->Pts_Num ()) if (pt[i % adat->Pts_Num ()] == 'I') continue;

    // 내부점이거나, interface위의 점인 경우
    if (pt[i % adat->Pts_Num ()] == 'C' || pt[i % adat->Pts_Num ()] == 'I' || pt[i % adat->Pts_Num ()] == 'M' || pt[i % adat->Pts_Num ()] == 'N' || pt[i % adat->Pts_Num ()].Condition () == 'D') {

      is_sol = i < adat->Pts_Num();

      // Neumann condition이 주어진 점의 두번째 식을 Dirichlet condition이 주어진 점과 같은 방식으로 계산
      if (pt[i % adat->Pts_Num ()] == 'N' && !is_sol) ND = true, pt[i % adat->Pts_Num ()].SetCondition ('D');

      // 내부점이거나, Interface인 경우에 해의 표현식의 계수들을 계산
      CalcRepresenCoef (&pt[i % adat->Pts_Num ()], &xdat, &ydat, is_sol);

      // 경계 근처에서의 해의 표현식을 처리
      TransposeBoundaryData (&pt[i % adat->Pts_Num ()], &xdat, &ydat, is_sol);

      // 행렬을 구성
      matps.MakeMatrixSystem (j, adat, &pt[i % adat->Pts_Num ()], &xdat, &ydat);

      if (ND) pt[i % adat->Pts_Num ()].SetCondition ('N'), ND = false;

      // 행렬을 만드는데 포함되는 점에 1을 더함
      j += 1;

    }

    // 모든점에 대한 현재 진행중인 점의 index를 알림
    // printf ("%s%lu%s%d", "\r", i + 1, "/", 2 * adat->Pts_Num ());
    ProgressBar (i + 1, 2 * adat->Pts_Num ());
  }
  // 행렬의 구성을 내보내기
  matps.ExportMatrixData (adat);
  // // 주위점의 정보를 출력
  // FILE *EWNS_output = fopen ("EWNS.dat", "w");
  // for (size_t i = 0; i < adat->Pts_Num (); i++) {
  //   pt[i].ExportEWNS (EWNS_output);
  // }
  // fclose (EWNS_output);
  // // 두번째 주위점의 정보를 출력
  // FILE *EWNS2nd_output = fopen ("EWNS2nd.dat", "w");
  // for (size_t i = 0; i < adat->Pts_Num (); i++) {
  //   pt[i].ExportEWNS2nd (EWNS2nd_output);
  // }
  // fclose (EWNS2nd_output);

  // 행렬의 구성이 끝났음을 알림
  printf ("\n%s\n", "Matrix made");

  // 행렬식의 계산이 시작됨을 알림
  printf ("\n%s\n", "Calcuate the matrix:");

  // 행렬식의 계산
  matps.calcMatrix (cdat, adat, pt);

  // 행렬식의 계산이 끝났음을 알림
  printf ("\n%s\n", "Matrix caculated");

  // 경계에서의 phi값을 부여
  // AssignPhivalue (adat, pt);

  // 미분의 값을 계산을 시작함을 알림
  printf ("\n%s\n", "Calcuate the differentiation:");

  for (size_t i = 0; i < adat->Pts_Num (); i++) {
      pt[i].Diff ('x')->SetValue (CalcDiff ('x', pt[i].Diff ('x'), &xdat, &ydat));
      pt[i].Diff ('y')->SetValue (CalcDiff ('y', pt[i].Diff ('y'), &xdat, &ydat));
      // printf ("%s%lu%s%d", "\r", i + 1, "/", adat->Pts_Num ());
      ProgressBar (i + 1, adat->Pts_Num ());
  }

  // 미분의 값의 계산이 끝났을 알림
  printf ("\n%s\n", "The differentiation is calculated");


  // 미분의 값을 계산을 시작함을 알림
  printf ("\n%s\n", "Calcuate the twice differentiation:");

  for (size_t i = 0; i < adat->Pts_Num (); i++) {
      pt[i].Diff ('x')->Diff ('x')->SetValue (CalcDiffDiff ('x', pt[i].Diff ('x')->Diff ('x'), &xdat, &ydat));
      pt[i].Diff ('y')->Diff ('y')->SetValue (CalcDiffDiff ('y', pt[i].Diff ('y')->Diff ('y'), &xdat, &ydat));
      // printf ("%s%lu%s%d", "\r", i + 1, "/", adat->Pts_Num ());
      ProgressBar (i + 1, adat->Pts_Num ());
  }

  // 미분의 값의 계산이 끝났을 알림
  printf ("\n%s\n", "The twice differentiation is calculated");

//  double result = ZeroValue;
//  double maxResult = ZeroValue;
//  int    maxIdx = -1;
//
//  for (size_t i = 0; i < adat->Pts_Num (); i++) {
//    // if (pt[i] == 'N') pt[i].SetCondition ('D'), ND = true;
//    matps.PrintDebuggingData (adat, &pt[i], &xdat, &ydat, &result);
//    if (maxResult < fabs (result)) maxResult = fabs (result), maxIdx = i;
//    // if (ND) pt[i].SetCondition ('N'), ND = false;
//  }
//  // if (pt[maxIdx] == 'N') pt[maxIdx].SetCondition ('D'), ND = true;
//  pt[maxIdx].Diff ('x')->Diff ('x')->PrintDebuggingData ("11111", &xdat, &ydat, true);
//  printf ("maxResult = %23.16e\n", maxResult);
//  // if (ND) pt[maxIdx].SetCondition ('N'), ND = false;

}

#endif
