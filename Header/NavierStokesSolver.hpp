#ifndef NAVIERSTOKESSOLVER
#define NAVIERSTOKESSOLVER

#include "Solver.hpp"

void NavierStokesSolver (ControlData *cdat, AxialData *adat, Point *ptU, Point *ptV, Point *ptP) {

  // 현재시각이 종료시각보다 작은 경우 반복
  while (ptP[0].NextTime ()) {

    // 현재시각을 update
    for (size_t i = 0; i < adat->Pts_Num (); i++) ptU[i].NextTime ();
    for (size_t i = 0; i < adat->Pts_Num (); i++) ptV[i].NextTime ();
    for (size_t i = 1; i < adat->Pts_Num (); i++) ptP[i].NextTime ();

    // u-velocity의 elliptic equation을 계산
    Solver (cdat, adat, ptU);
    // v-velocity의 elliptic equation을 계산
    Solver (cdat, adat, ptV);

  }

}

#endif
