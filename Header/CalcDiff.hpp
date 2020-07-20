#ifndef CALCDIFF_H
#define CALCDIFF_H

#include "CalcAtNeumannPt.hpp"

/* 해의 미분의 표현식의 계수를 계산 */
void SettingDiffCoefficient(Point *pt, xData *xdat, yData *ydat,
                            const double xm, const double xb, const double xp,
                            const double ym, const double yb, const double yp,
                            const int bdx, const int bdy,
                            const double mpx1, const double mpx2,
                            const double mpy1, const double mpy2) {

    // 국소 x-축선에서 u의 현재점에서의 계수
    // xdat->Cu = - 1.0E0;
    // 국소 x-축선에서 u의 국소 x-축선에서 끝점에서의 계수
    xdat->Eu = greens_coefficient_ttau(xp, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 u의 국소 x-축선에서 시작점에서의 계수
    xdat->Wu = greens_coefficient_ttau(xm, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 phi의 현재점에서의 계수
    if (!pt->IsBoundary('W')) xdat->Cphi += greens_integral_tau(2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    if (!pt->IsBoundary('E')) xdat->Cphi += greens_integral_tau(3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 phi의 국소 x-축선에서 끝점에서의 계수
    if (!pt->IsBoundary('E')) xdat->Ephi = greens_integral_tau(4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 phi의 국소 x-축선에서 시작점에서의 계수
    if (!pt->IsBoundary('W')) xdat->Wphi = greens_integral_tau(1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 f의 현재점에서의 계수
    if (!pt->IsBoundary('W')) xdat->Cf += greens_integral_tau(2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    if (!pt->IsBoundary('E')) xdat->Cf += greens_integral_tau(3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 f의 국소 x-축선에서 끝점에서의 계수
    if (!pt->IsBoundary('E')) xdat->Ef = greens_integral_tau(4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 f의 국소 x-축선에서 시작점에서의 계수
    if (!pt->IsBoundary('W')) xdat->Wf = greens_integral_tau(1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 y-축선에서 u의 현재점에서의 계수
    // ydat->Cu = - 1.0E0;
    // 국소 y-축선에서 u의 국소 y-축선에서 끝점에서의 계수
    ydat->Nu = greens_coefficient_ttau(yp, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 u의 국소 y-축선에서 시작점에서의 계수
    ydat->Su = greens_coefficient_ttau(ym, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 phi의 현재점에서의 계수
    if (!pt->IsBoundary('S')) ydat->Cphi += -greens_integral_tau(2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    if (!pt->IsBoundary('N')) ydat->Cphi += -greens_integral_tau(3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 phi의 국소 y-축선에서 끝점에서의 계수
    if (!pt->IsBoundary('N')) ydat->Nphi = -greens_integral_tau(4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 phi의 국소 y-축선에서 시작점에서의 계수
    if (!pt->IsBoundary('S')) ydat->Sphi = -greens_integral_tau(1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 f의 현재점에서의 계수
    if (!pt->IsBoundary('S')) ydat->Cf += greens_integral_tau(2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    if (!pt->IsBoundary('N')) ydat->Cf += greens_integral_tau(3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 f의 국소 y-축선에서 끝점에서의 계수
    if (!pt->IsBoundary('N')) ydat->Nf = greens_integral_tau(4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 f의 국소 y-축선에서 시작점에서의 계수
    if (!pt->IsBoundary('S')) ydat->Sf = greens_integral_tau(1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);

    if (::SolverType)
        SettingDiffTimeIntegation(HeadVelocity(pt), xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1,
                                  mpy2);
}

void SettingDiffDiffCoefficient(Point *pt, xData *xdat, yData *ydat,
                                const double xm, const double xb, const double xp,
                                const double ym, const double yb, const double yp,
                                const int bdx, const int bdy,
                                const double mpx1, const double mpx2,
                                const double mpy1, const double mpy2) {

    // 국소 x-축선에서 u의 현재점에서의 계수
    // xdat->Cu = - 1.0E0;
    // 국소 x-축선에서 u의 국소 x-축선에서 끝점에서의 계수
    xdat->Edu = greens_coefficient_ttau(xp, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 u의 국소 x-축선에서 시작점에서의 계수
    xdat->Wdu = greens_coefficient_ttau(xm, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 phi의 현재점에서의 계수
    if (!pt->IsBoundary('W')) xdat->Cphi += -greens_integral_ttau(2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    if (!pt->IsBoundary('E')) xdat->Cphi += -greens_integral_ttau(3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 phi의 국소 x-축선에서 끝점에서의 계수
    if (!pt->IsBoundary('E')) xdat->Ephi = -greens_integral_ttau(4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 phi의 국소 x-축선에서 시작점에서의 계수
    if (!pt->IsBoundary('W')) xdat->Wphi = -greens_integral_ttau(1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 f의 현재점에서의 계수
    if (!pt->IsBoundary('W')) xdat->Cf += -greens_integral_ttau(2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    if (!pt->IsBoundary('E')) xdat->Cf += -greens_integral_ttau(3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 f의 국소 x-축선에서 끝점에서의 계수
    if (!pt->IsBoundary('E')) xdat->Ef = -greens_integral_ttau(4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 x-축선에서 f의 국소 x-축선에서 시작점에서의 계수
    if (!pt->IsBoundary('W')) xdat->Wf = -greens_integral_ttau(1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
    // 국소 y-축선에서 u의 현재점에서의 계수
    // ydat->Cu = - 1.0E0;
    // 국소 y-축선에서 u의 국소 y-축선에서 끝점에서의 계수
    ydat->Ndu = greens_coefficient_ttau(yp, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 u의 국소 y-축선에서 시작점에서의 계수
    ydat->Sdu = greens_coefficient_ttau(ym, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 phi의 현재점에서의 계수
    if (!pt->IsBoundary('S')) ydat->Cphi += greens_integral_ttau(2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    if (!pt->IsBoundary('N')) ydat->Cphi += greens_integral_ttau(3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 phi의 국소 y-축선에서 끝점에서의 계수
    if (!pt->IsBoundary('N')) ydat->Nphi = greens_integral_ttau(4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 phi의 국소 y-축선에서 시작점에서의 계수
    if (!pt->IsBoundary('S')) ydat->Sphi = greens_integral_ttau(1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 f의 현재점에서의 계수
    if (!pt->IsBoundary('S')) ydat->Cf += -greens_integral_ttau(2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    if (!pt->IsBoundary('N')) ydat->Cf += -greens_integral_ttau(3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 f의 국소 y-축선에서 끝점에서의 계수
    if (!pt->IsBoundary('N')) ydat->Nf = -greens_integral_ttau(4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
    // 국소 y-축선에서 f의 국소 y-축선에서 시작점에서의 계수
    if (!pt->IsBoundary('S')) ydat->Sf = -greens_integral_ttau(1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);

    if (::SolverType)
        SettingDiffDiffTimeIntegation(HeadVelocity(pt), xdat, ydat, xm, xb, xp, ym, yb, yp, bdx, bdy, mpx1, mpx2, mpy1,
                                      mpy2);
}

/* 해의 미분을 계산 */
double CalcDiff(char xy, Point *pt, xData *xdat, yData *ydat) {
    double arrEnt[26] = {ZeroValue,}, arrVal[26] = {ZeroValue,}, returnval = ZeroValue;
    char azimuth[12][3] = {"EE", "WW", "NN", "SS", "EN", "ES", "WN", "WS", "NE", "NW", "SE", "SW"};
    int idx = 0;
    // char errorMassage[256];

    // 해의 미분의 표현식의 계수를 계산
    CalcRepresenCoef(pt, xdat, ydat, true);
    // 현재점에서의 해의 값
    arrVal[0] = HeadVelocity(pt)->Value(), arrVal[13] = HeadVelocity(pt)->Phi()->Value();
    for (const auto &i : azimuth) {
        idx++;
        if (HeadVelocity(pt)->EWNS(i[0], i[1]))
            arrVal[idx] = HeadVelocity(pt)->EWNS(i[0], i[1])->Value(), arrVal[idx + 13] = HeadVelocity(pt)->EWNS(i[0],
                                                                                                                 i[1])->Phi()->Value();
    }

    // x-성분의 미분을 계산하는 경우
    if (xy == 'X' || xy == 'x') {
        returnval = -xdat->F,
        arrEnt[0] = xdat->Cu, arrEnt[1] = xdat->Eu, arrEnt[2] = xdat->Wu, arrEnt[5] = xdat->ENu, arrEnt[6] = xdat->ESu, arrEnt[7] = xdat->WNu, arrEnt[8] = xdat->WSu,
        arrEnt[13] = xdat->Cphi, arrEnt[14] = xdat->Ephi, arrEnt[15] = xdat->Wphi, arrEnt[18] = xdat->ENphi, arrEnt[19] = xdat->ESphi, arrEnt[20] = xdat->WNphi, arrEnt[21] = xdat->WSphi;
    }

    // y-성분의 미분을 계산하는 경우
    if (xy == 'Y' || xy == 'y') {
        returnval = -ydat->F,
        arrEnt[0] = ydat->Cu, arrEnt[3] = ydat->Nu, arrEnt[4] = ydat->Su, arrEnt[9] = ydat->NEu, arrEnt[10] = ydat->NWu, arrEnt[11] = ydat->SEu, arrEnt[12] = ydat->SWu,
        arrEnt[13] = ydat->Cphi, arrEnt[16] = ydat->Nphi, arrEnt[17] = ydat->Sphi, arrEnt[22] = ydat->NEphi, arrEnt[23] = ydat->NWphi, arrEnt[24] = ydat->SEphi, arrEnt[25] = ydat->SWphi;
    }

    /* --------------------------- */

//  for (size_t i = 0; i < 26; i++) {
//      if (isnan(arrEnt[i])) {
//          printf ("%s%02lu%s%f\n", "arrEnt[", i, "] = ", arrEnt[i]);
//          printf ("%s%d\n", "index = ", pt->Index ());
//          PrintPtsInfo (pt);
//      }
//  }
    // if (IsEqualDouble (arrEnt[0], 0.0)) {
    //   printf ("%s%02d%s%f\n", "arrEnt[", 0, "] = ", arrEnt[0]);
    //   printf ("%s%d\n", "index = ", pt->Index ());
    //   PrintPtsInfo (pt);
    // }

    /* --------------------------- */

    // 해의 미분의 표현식의 계수와 해의 값을 곲해서 더한다
    for (size_t i = 0; i < 26; i++) returnval += arrEnt[i] * arrVal[i];

    // if (returnval != returnval) for (size_t i = 1; i < 26; i++) printf ("arrEnt[%zu] = %23.16e\tarrVal[%zu] = %23.16e\n", i, arrEnt[i], i, arrVal[i]);
    // if (returnval != returnval) pt->PrintDebuggingData ("11111", xdat, ydat, true), sprintf (errorMassage, "CalcDiff:: returnval = NaN"), PrintError (errorMassage);

    if (isnan(returnval)) {
        for (const auto &i : {'E', 'W', 'N', 'S'}) {
            if (HeadVelocity(pt)->EWNS(i, i))
                if (HeadVelocity(pt)->EWNS(i, i)->Condition() == 'C')
                    return HeadVelocity(pt)->EWNS(i, i)->Diff(xy)->Value();
        }
    }

    // 현재점의 계수로 나누어서 return
    return returnval;
}

double CalcDiffDiff(char xy, Point *pt, xData *xdat, yData *ydat) {
    double arrEnt[39] = {ZeroValue,}, arrVal[39] = {ZeroValue,}, returnval = ZeroValue;
    char azimuth[12][3] = {"EE", "WW", "NN", "SS", "EN", "ES", "WN", "WS", "NE", "NW", "SE", "SW"};
    int idx = 0;
    // char errorMassage[256];

    // 해의 미분의 표현식의 계수를 계산
    CalcRepresenCoef(pt, xdat, ydat, true);
    // 현재점에서의 해의 값
    arrVal[0] = HeadVelocity(pt)->Value(), arrVal[13] = HeadVelocity(pt)->Phi()->Value(),
    arrVal[26] = HeadVelocity(pt)->Diff(xy)->Value();
    for (const auto &i : azimuth) {
        idx++;
        if (HeadVelocity(pt)->EWNS(i[0], i[1])) {
            arrVal[idx] = HeadVelocity(pt)->EWNS(i[0], i[1])->Value(),
            arrVal[idx + 13] = HeadVelocity(pt)->EWNS(i[0], i[1])->Phi()->Value(),
            arrVal[idx + 26] = HeadVelocity(pt)->EWNS(i[0], i[1])->Diff(xy)->Value();
        }
    }

    // x-성분의 미분을 계산하는 경우
    if (xy == 'X' || xy == 'x') {
        returnval = -xdat->F - (5.0E-1 * HeadVelocity(pt)->F() + arrVal[13]) / HeadVelocity(pt)->MaterialProperty();
        if (::SolverType) {
            returnval += (5.0E-1 * (HeadVelocity(pt)->Value() - HeadVelocity(pt)->Pre()->Value()) /
                          HeadVelocity(pt)->Dt()) / HeadVelocity(pt)->MaterialProperty();
            returnval -= HeadVelocity(pt)->Diff('x')->Diff('x')->Pre()->Value();
        }
        arrEnt[0] = xdat->Cu, arrEnt[1] = xdat->Eu, arrEnt[2] = xdat->Wu, arrEnt[5] = xdat->ENu, arrEnt[6] = xdat->ESu, arrEnt[7] = xdat->WNu, arrEnt[8] = xdat->WSu,
        arrEnt[13] = xdat->Cphi, arrEnt[14] = xdat->Ephi, arrEnt[15] = xdat->Wphi, arrEnt[18] = xdat->ENphi, arrEnt[19] = xdat->ESphi, arrEnt[20] = xdat->WNphi, arrEnt[21] = xdat->WSphi,
        arrEnt[26] = xdat->Cdu, arrEnt[27] = xdat->Edu, arrEnt[28] = xdat->Wdu, arrEnt[31] = xdat->ENdu, arrEnt[32] = xdat->ESdu, arrEnt[33] = xdat->WNdu, arrEnt[34] = xdat->WSdu;
    }

    // y-성분의 미분을 계산하는 경우
    if (xy == 'Y' || xy == 'y') {
        returnval = -ydat->F - (5.0E-1 * HeadVelocity(pt)->F() - arrVal[13]) / HeadVelocity(pt)->MaterialProperty();
        if (::SolverType) {
            returnval += (5.0E-1 * (HeadVelocity(pt)->Value() - HeadVelocity(pt)->Pre()->Value()) /
                          HeadVelocity(pt)->Dt()) / HeadVelocity(pt)->MaterialProperty();
            returnval -= HeadVelocity(pt)->Diff('y')->Diff('y')->Pre()->Value();
        }
        arrEnt[0] = ydat->Cu, arrEnt[3] = ydat->Nu, arrEnt[4] = ydat->Su, arrEnt[9] = ydat->NEu, arrEnt[10] = ydat->NWu, arrEnt[11] = ydat->SEu, arrEnt[12] = ydat->SWu,
        arrEnt[13] = ydat->Cphi, arrEnt[16] = ydat->Nphi, arrEnt[17] = ydat->Sphi, arrEnt[22] = ydat->NEphi, arrEnt[23] = ydat->NWphi, arrEnt[24] = ydat->SEphi, arrEnt[25] = ydat->SWphi,
        arrEnt[26] = ydat->Cdu, arrEnt[29] = ydat->Ndu, arrEnt[30] = ydat->Sdu, arrEnt[35] = ydat->NEdu, arrEnt[36] = ydat->NWdu, arrEnt[37] = ydat->SEdu, arrEnt[38] = ydat->SWdu;
    }

    /* --------------------------- */

//  for (size_t i = 0; i < 26; i++) {
//      if (isnan(arrEnt[i])) {
//          printf ("%s%02lu%s%f\n", "arrEnt[", i, "] = ", arrEnt[i]);
//          printf ("%s%d\n", "index = ", HeadVelocity (pt)->Index ());
//          PrintPtsInfo (HeadVelocity (pt));
//      }
//  }
//  for (size_t i = 0; i < 26; i++) {
//    if (isinf (arrEnt[i])) {
//      printf ("%s%02lu%s%f\n", "arrEnt[", i, "] = ", arrEnt[i]);
//      printf ("%s%d\n", "index = ", HeadVelocity (pt)->Index ());
//      PrintPtsInfo (HeadVelocity (pt));
//    }
//  }
    // if (IsEqualDouble (arrEnt[0], ZeroValue)) {
    //   printf ("%s%02d%s%f\n", "arrEnt[", 0, "] = ", arrEnt[0]);
    //   printf ("%s%d\n", "index = ", HeadVelocity (pt)->Index ());
    //   PrintPtsInfo (HeadVelocity (pt));
    // }
//  if (returnval != returnval) {
//    printf ("%s%f\n", "returnval = ", returnval);
//    printf ("%s%d\n", "index = ", HeadVelocity (pt)->Index ());
//    PrintPtsInfo (HeadVelocity (pt));
//  }

    /* --------------------------- */


    // 해의 미분의 표현식의 계수와 해의 값을 곲해서 더한다
    for (size_t i = 0; i < 39; i++) returnval += arrEnt[i] * arrVal[i];

    // if (returnval != returnval) pt->PrintDebuggingData ("11111", xdat, ydat, true), sprintf (errorMassage, "CalcDiffDiff:: returnval = NaN"), PrintError (errorMassage);

    if (isnan(returnval)) {
        for (const auto &i : {'E', 'W', 'N', 'S'}) {
            if (HeadVelocity(pt)->EWNS(i, i))
                if (HeadVelocity(pt)->EWNS(i, i)->Condition() == 'C')
                    return HeadVelocity(pt)->EWNS(i, i)->Diff(xy)->Diff(xy)->Value();
        }
    }

    // 현재점의 계수로 나누어서 return
    return returnval;
}

#endif
