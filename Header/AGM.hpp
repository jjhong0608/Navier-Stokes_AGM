#ifndef AGM_HPP
#define AGM_HPP

#include <utility>

#include "NavierStokesSolver.hpp"

int
AGM(string AGL_output_file, string AGM_output_file, double initialTime, double terminalTime, int timeStep, double dt) {
    // AGM solver의 시작을 알림
    printf("+");
    for (int i = 0; i < 35; i++) printf("-");
    printf("+");
    printf("\n%s\n", "| AGM Navier-Stokes equation solver |");
    printf("+");
    for (int i = 0; i < 35; i++) printf("-");
    printf("+\n");

    SettingSolver("Heat");

    // 제어정보를 가지는 class의 선언
    ControlData cdat;

    // 제어정보를 파일로부터 읽어들임
    cdat.LoadCtrlData(std::move(AGL_output_file), std::move(AGM_output_file), initialTime, terminalTime, timeStep, dt);

    // 읽어들인 제어정보를 출력
    cdat.ShowCtrlData();

    // 점과 축선의 정보를 가지는 class의 선언
    AxialData adat;

    // 축선생성기로부터 만들어진 파일로부터 점과 축선의 정보를 읽어들임
    adat.LoadAxialData(cdat.Axialfile());

    // 각 점의 왼쪽, 오른쪽, 위쪽, 아래쪽의 index를 저장하는 변수의 초기화
    adat.SortEWNS();

    // 속도 u를 저장하는 점의 class의 선언
    auto *ptU = new Point[adat.Pts_Num()];
    // 속도 v를 저장하는 점의 class의 선언
    auto *ptV = new Point[adat.Pts_Num()];
    // 압력 p를 저장하는 점의 class의 선언
    auto *ptP = new Point[adat.Pts_Num()];
    // 속도 u의 x성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_diff_u_x = new Point[adat.Pts_Num()];
    // 속도 u의 y성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_diff_u_y = new Point[adat.Pts_Num()];
    // 속도 v의 x성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_diff_v_x = new Point[adat.Pts_Num()];
    // 속도 v의 y성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_diff_v_y = new Point[adat.Pts_Num()];
    // 속도 u의 x성분으로의 미분의 x성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_diff_u_xx = new Point[adat.Pts_Num()];
    // 속도 u의 y성분으로의 미분의 y성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_diff_u_yy = new Point[adat.Pts_Num()];
    // 속도 v의 x성분으로의 미분의 x성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_diff_v_xx = new Point[adat.Pts_Num()];
    // 속도 v의 y성분으로의 미분의 y성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_diff_v_yy = new Point[adat.Pts_Num()];
    // 중간 속도 u를 저장하는 점의 class의 선언
    auto *pt_hat_u = new Point[adat.Pts_Num()];
    // 중간 속도 v를 저장하는 점의 class의 선언
    auto *pt_hat_v = new Point[adat.Pts_Num()];
    // 이전의 속도 u를 저장하는 점의 class의 선언
    auto *pt_pre_u = new Point[adat.Pts_Num()];
    // 이전의 속도 v를 저장하는 점의 class의 선언
    auto *pt_pre_v = new Point[adat.Pts_Num()];
    // 이전의 속도 u의 x성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_pre_diff_u_x = new Point[adat.Pts_Num()];
    // 이전의 속도 u의 y성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_pre_diff_u_y = new Point[adat.Pts_Num()];
    // 이전의 속도 v의 x성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_pre_diff_v_x = new Point[adat.Pts_Num()];
    // 이전의 속도 v의 y성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_pre_diff_v_y = new Point[adat.Pts_Num()];
    // 이전의 속도 u의 x성분으로의 미분의 x성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_pre_diff_u_xx = new Point[adat.Pts_Num()];
    // 이전의 속도 u의 y성분으로의 미분의 y성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_pre_diff_u_yy = new Point[adat.Pts_Num()];
    // 이전의 속도 v의 x성분으로의 미분의 x성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_pre_diff_v_xx = new Point[adat.Pts_Num()];
    // 이전의 속도 v의 y성분으로의 미분의 y성분으로의 미분을 저장하는 점의 class의 선언
    auto *pt_pre_diff_v_yy = new Point[adat.Pts_Num()];
    // 전전의 속도 u를 저장하는 점의 class의 선언
    auto *pt_old_u = new Point[adat.Pts_Num()];
    // 전전의 속도 v를 저장하는 점의 class의 선언
    auto *pt_old_v = new Point[adat.Pts_Num()];
    // phi를 저장하는 점의 class의 선언
    auto *pt_phi = new Point[adat.Pts_Num()];
    // psi를 저장하는 점의 class의 선언
    auto *pt_psi = new Point[adat.Pts_Num()];
    // 이전의 phi를 저장하는 점의 class의 선언
    auto *pt_pre_phi = new Point[adat.Pts_Num()];
    // 이전의 psi를 저장하는 점의 class의 선언
    auto *pt_pre_psi = new Point[adat.Pts_Num()];

    SettingPoints(&adat, ptU, ptV, ptP, pt_diff_u_x, pt_diff_u_y, pt_diff_v_x, pt_diff_v_y, pt_diff_u_xx, pt_diff_u_yy,
                  pt_diff_v_xx, pt_diff_v_yy, pt_hat_u, pt_hat_v, pt_pre_u, pt_pre_v, pt_old_u, pt_old_v,
                  pt_pre_diff_u_x, pt_pre_diff_u_y, pt_pre_diff_u_xx, pt_pre_diff_u_yy, pt_pre_diff_v_x,
                  pt_pre_diff_v_y, pt_pre_diff_v_xx, pt_pre_diff_v_yy, pt_phi, pt_psi, pt_pre_phi, pt_pre_psi);

    // 모든점에 대해서 반복
    for (int i = 0; i < adat.Pts_Num(); i++) {

        // 종료시각을 저장
        ptU[i].SetTerminalTime(cdat.TerminalTime());
        ptV[i].SetTerminalTime(cdat.TerminalTime());
        ptP[i].SetTerminalTime(cdat.TerminalTime());

        // 현재시각을 저장
        ptU[i].SetTime(cdat.InitialTime());
        ptV[i].SetTime(cdat.InitialTime());
        ptP[i].SetTime(cdat.InitialTime());

        // time step을 저장
        ptU[i].SetTimeStep(cdat.TimeStep());
        ptV[i].SetTimeStep(cdat.TimeStep());
        ptP[i].SetTimeStep(cdat.TimeStep());

        ptU[i].SetDt(cdat.Dt());
        ptV[i].SetDt(cdat.Dt());
        ptP[i].SetDt(cdat.Dt());

        ProgressBar(i + 1, adat.Pts_Num());
    }

    // 점의 정보를 만들기 시작했다는 알림
    printf("\n%s\n", "Read point data:");

    // 모든점의 대해서 반복
    for (int i = 0; i < adat.Pts_Num(); i++) {

        // 각 점의 x-좌표와 y-좌표를 읽어들임
        ptP[i].SetCoordinate(adat.Pts(i, 'x'), adat.Pts(i, 'y'));

        // 각 점의 index를 읽어들임
        ptP[i].SetIndex(i);

        // 각 점의 경계조건을 읽어들임
        ptP[i].SetCondition(adat.Boundarycondition(i));

        // 각 점의 경계조건의 값을 읽어들임
        ptP[i].SetBoundaryvalue(ZeroValue);

        // 각 점의 값을 0으로 초기화
        ptP[i].SetValue(ZeroValue);

        // 각 점의 우변의 f값을 읽어들임
        ptP[i].SetF(ZeroValue);

        // 각 점의 conductivity를 읽어들임
        ptP[i].SetMaterialProperty(adat.MaterialProperty(i));

        // 각 점의 normal vector를 읽어들임
        ptP[i].SetNormal('x', adat.Normal(i, 'x'));
        ptP[i].SetNormal('y', adat.Normal(i, 'y'));

        // 모든점에 대한 현재 진행중인 점의 index를 알림
        ProgressBar(i + 1, adat.Pts_Num());
    }
    printf("\n");

    // 모든점의 대해서 반복
    for (int i = 0; i < adat.Pts_Num(); i++) {

        // 각점의 주의의 점들을 찾음
        ptP[i].FindAxialElement(&adat, ptP);

        // 모든점에 대한 현재 진행중인 점의 index를 알림
        ProgressBar(i + 1, adat.Pts_Num());
    }
    printf("\n");

    // 모든점의 대해서 반복
    for (int i = 0; i < adat.Pts_Num(); i++) {

        // 각점의 두번째 주의의 점들을 찾음
        ptP[i].Find2ndAxialElement(&adat, ptP);

        // 모든점에 대한 현재 진행중인 점의 index를 알림
        ProgressBar(i + 1, adat.Pts_Num());
    }
    printf("\n");

    // 모든점의 대해서 반복
    for (int i = 0; i < adat.Pts_Num(); i++) {

        // 각 점을 포함하는 국소축선의 시작점과 끝점의 좌표를 찾음
        ptP[i].SetMinMax();

        // 경계의 주위점들을 찾음
        ptP[i].FindBoundaryElement();

        // Interface위의 점인 경우, 양쪽의 conductivity의 동일여부 확인
        ptP[i].IsInterface();

        // 모든점에 대한 현재 진행중인 점의 index를 알림
        ProgressBar(i + 1, adat.Pts_Num());
    }
    printf("\n");

    // 모든점의 대해서 반복
    for (int i = 0; i < adat.Pts_Num(); i++) {

        // 각 점을 포함하는 국소축선의 시작점과 끝점의 좌표를 다시 찾음
        ptP[i].SetMinMax();

        // 모든점에 대한 현재 진행중인 점의 index를 알림
        ProgressBar(i + 1, adat.Pts_Num());
    }
    printf("\n");

    // 모든점의 index와 내부점의 index, 모든점의 index와 phi의 값을 계산하는 점의 index를 저장하는 변수의 allocation
    adat.AllocatePhipts(adat.In_Pts_Num() - CountInterface(&adat, ptP));

    // 점들을 모든점과 내부점, phi의 값을 계산하는 점들로 분류
    sortPts(&adat, ptP);

    // // 점과 축선의 정보를 내보내기
    // adat.ExportAxialData (&cdat);

    // 모든 점에 대해서 반복
    for (int i = 0; i < adat.Pts_Num(); i++) {

        // u값을 저장하는 점들의 좌표를 저장
        ptU[i].SetCoordinate(adat.Pts(i, 'x'), adat.Pts(i, 'y'));
        // v값을 저장하는 점들의 좌표를 저장
        ptV[i].SetCoordinate(adat.Pts(i, 'x'), adat.Pts(i, 'y'));

        // u값을 저장하는 점들의 index를 저장
        ptU[i].SetIndex(ptP[i].Index());
        // v값을 저장하는 점들의 index를 저장
        ptV[i].SetIndex(ptP[i].Index());

        // u값을 저장하는 점들의 경계조건을 읽어들임
        ptU[i].SetCondition(ptP[i].Condition());
        // v값을 저장하는 점들의 경계조건을 읽어들임
        ptV[i].SetCondition(ptP[i].Condition());

        ptU[i].Pre()->SetValue(u_ftn_Dirichlet(&ptU[i]));
        ptV[i].Pre()->SetValue(u_ftn_Dirichlet(&ptV[i]));

        ptU[i].Phi()->Pre()->SetValue(u_ftn(ptU[i].Phi()));
        ptV[i].Phi()->Pre()->SetValue(u_ftn(ptV[i].Phi()));

        ptU[i].Diff('x')->Diff('x')->Pre()->SetValue(u_ftn(ptU[i].Diff('x')->Diff('x')));
        ptU[i].Diff('y')->Diff('y')->Pre()->SetValue(u_ftn(ptU[i].Diff('y')->Diff('y')));

        ptU[i].SetTime(cdat.InitialTime() + cdat.Dt());
        ptV[i].SetTime(cdat.InitialTime() + cdat.Dt());
        ptP[i].SetTime(cdat.InitialTime() + cdat.Dt());

        // u값을 저장하는 점들의 경계조건의 값을 읽어들임
        ptU[i].SetBoundaryvalue(u_ftn(&ptU[i]));
        // v값을 저장하는 점들의 경계조건의 값을 읽어들임
        ptV[i].SetBoundaryvalue(u_ftn(&ptV[i]));

        // u값을 저장하는 점들의 값을 초기화 (경계조건의 값으로 초기화, 내부점인 경우 0으로 초기화)
        ptU[i].SetValue(u_ftn_Dirichlet(&ptU[i]));
        ptV[i].SetValue(u_ftn_Dirichlet(&ptV[i]));

        // u값을 저장하는 점들의우변의 f값을 읽어들임
        ptU[i].SetF(f_ftn(&ptU[i]));

        // v값을 저장하는 점들의 우변의 f값을 읽어들임
        ptV[i].SetF(f_ftn(&ptV[i]));

        // Heat equation의 경우 material property에 1/2를 곱한다
        double mp = ::SolverType ? ptP[i].MaterialProperty() * 5.0E-1 : ptP[i].MaterialProperty();

        // u값을 저장하는 점들의 conductivity를 읽어들임
        ptU[i].SetMaterialProperty(mp);
        // v값을 저장하는 점들의 conductivity를 읽어들임
        ptV[i].SetMaterialProperty(mp);

        // 각 점의 normal vector를 읽어들임
        ptU[i].SetNormal('x', adat.Normal(i, 'x'));
        ptU[i].SetNormal('y', adat.Normal(i, 'y'));
        ptV[i].SetNormal('x', adat.Normal(i, 'x'));
        ptV[i].SetNormal('y', adat.Normal(i, 'y'));

        // u값을 저장하는 점의 오른쪽, 왼쪽, 위쪽, 아래쪽을 대입
        SetAllEWNS(&adat, &ptP[i], &ptU[i], ptU);
        // u값을 저장하는 점의 두번재 오른쪽, 왼쪽, 위쪽, 아래쪽을 대입
        SetAllEWNS2nd(&adat, &ptP[i], &ptU[i], ptU);

        // v값을 저장하는 점의 오른쪽, 왼쪽, 위쪽, 아래쪽을 대입
        SetAllEWNS(&adat, &ptP[i], &ptV[i], ptV);
        // v값을 저장하는 점의 두번재 오른쪽, 왼쪽, 위쪽, 아래쪽을 대입
        SetAllEWNS2nd(&adat, &ptP[i], &ptV[i], ptV);

        // 모든점에 대한 현재 진행중인 점의 index를 알림
        ProgressBar(i + 1, adat.Pts_Num());
    }
    printf("\n");

    // 모든점의 대해서 반복
    for (int i = 0; i < adat.Pts_Num(); i++) {

        // 각 점을 포함하는 국소축선의 시작점과 끝점의 좌표를 다시 찾음
        ptU[i].SetMinMax();
        ptV[i].SetMinMax();

        // 모든점에 대한 현재 진행중인 점의 index를 알림
        ProgressBar(i + 1, adat.Pts_Num());

    }
    printf("\n");

    // 각 점의 정보가 다 입력되었음을 알림
    printf("\n%s\n", "All points made");

    // u-velocity의 elliptic equation을 계산
    SettingSolver("heat"), Solver(&cdat, &adat, ptU);
    // v-velocity의 elliptic equation을 계산
    // SettingSolver ("Heat"), Solver (&cdat, &adat, ptV);

    while (ptU[0].Time() < cdat.TerminalTime() - 5.0E-1 * cdat.Dt()) {

        for (int i = 0; i < adat.Pts_Num(); i++) {
            ptU[i].Pre()->SetValue(ptU[i].Value());
            // ptV[i].Pre ()->SetValue (ptV[i].Vale ());
            ptU[i].Diff('x')->Diff('x')->Pre()->SetValue(ptU[i].Diff('x')->Diff('x')->Value());
            ptU[i].Diff('y')->Diff('y')->Pre()->SetValue(ptU[i].Diff('y')->Diff('y')->Value());
            ptU[i].Phi()->Pre()->SetValue(ptU[i].Phi()->Value());
            ptV[i].Phi()->Pre()->SetValue(ptV[i].Phi()->Value());


            ptU[i].SetTime(ptU[i].Time() + cdat.Dt());
            ptV[i].SetTime(ptV[i].Time() + cdat.Dt());
            ptP[i].SetTime(ptP[i].Time() + cdat.Dt());
            // u값을 저장하는 점들의 경계조건의 값을 읽어들임
            ptU[i].SetBoundaryvalue(u_ftn(&ptU[i]));
            // ptU[i].SetBoundaryvalue (adat.Boundaryvalue (i));
            // v값을 저장하는 점들의 경계조건의 값을 읽어들임
            ptV[i].SetBoundaryvalue(u_ftn(&ptV[i]));

            ptU[i].SetF(f_ftn(&ptU[i]));
            ptV[i].SetF(f_ftn(&ptV[i]));
            // ptV[i].SetBoundaryvalue (adat.Boundaryvalue (i));
        }

        // u-velocity의 elliptic equation을 계산
        SettingSolver("Heat"), Solver(&cdat, &adat, ptU);
        // v-velocity의 elliptic equation을 계산
        // SettingSolver ("Heat"), Solver (&cdat, &adat, ptV);

        printf("Time = %.4f / %.4f\n", ptU[0].Time(), cdat.TerminalTime());
    }

    // // u-velocity의 elliptic equation을 계산
    // EllipticSolver (&cdat, &adat, ptU);
    // // v-velocity의 elliptic equation을 계산
    // EllipticSolver (&cdat, &adat, ptV);

    auto *pt_diff_u = new Point[adat.Pts_Num()];
    // Point *pt_diff_v = new Point[adat.Pts_Num ()];

    auto *pt_diff_diff_u = new Point[adat.Pts_Num()];
    // Point *pt_diff_diff_v = new Point[adat.Pts_Num ()];

    for (int i = 0; i < adat.Pts_Num(); i++) {
        pt_diff_u[i].SetMark("diff_u");
        pt_diff_u[i].SetVelocity('u', &ptU[i]);
        pt_diff_u[i].SetPressure(&ptP[i]);
        pt_diff_u[i].SetValue(sqrt(pt_diff_u_x[i].Value() * pt_diff_u_x[i].Value() +
                                   pt_diff_u_y[i].Value() * pt_diff_u_y[i].Value()));
    }

    // for (int i = 0; i < adat.Pts_Num (); i++)
    // pt_diff_v[i].SetMark ("diff_v"),
    // pt_diff_v[i].SetVelocity ('v', &ptV[i]),
    // pt_diff_v[i].SetPressure (&ptP[i]),
    // pt_diff_v[i].SetValue (sqrt (pt_diff_v_x[i].Value () * pt_diff_v_x[i].Value () + pt_diff_v_y[i].Value () * pt_diff_v_y[i].Value ()));

    for (int i = 0; i < adat.Pts_Num(); i++) {
        pt_diff_diff_u[i].SetMark("diff_diff_u");
        pt_diff_diff_u[i].SetVelocity('u', &ptU[i]);
        pt_diff_diff_u[i].SetPressure(&ptP[i]);
        pt_diff_diff_u[i].SetValue(sqrt(pt_diff_u_xx[i].Value() * pt_diff_u_xx[i].Value() +
                                        pt_diff_u_yy[i].Value() * pt_diff_u_yy[i].Value()));
    }

    // for (int i = 0; i < adat.Pts_Num (); i++)
    // pt_diff_diff_v[i].SetMark ("diff_diff_v"),
    // pt_diff_diff_v[i].SetVelocity ('v', &ptV[i]),
    // pt_diff_diff_v[i].SetPressure (&ptP[i]),
    // pt_diff_diff_v[i].SetValue (sqrt (pt_diff_v_xx[i].Value () * pt_diff_v_xx[i].Value () + pt_diff_v_yy[i].Value () * pt_diff_v_yy[i].Value ()));

    // 해의 error를 계산해서 출력
    Calc_Error(&adat, ptU);
    // Calc_Error (&adat, ptV);
    Calc_Error(&adat, pt_phi);
    // Calc_Error (&adat, pt_psi);
    Calc_Error(&adat, pt_diff_u);
    // Calc_Error (&adat, pt_diff_v);
    Calc_Error(&adat, pt_diff_diff_u);
    // Calc_Error (&adat, pt_diff_diff_v);
    // Calc_Error (&adat, pt_diff_u_x);
    // Calc_Error (&adat, pt_diff_u_y);
    // Calc_Error (&adat, pt_diff_v_x);
    // Calc_Error (&adat, pt_diff_v_y);
    // Calc_Error (&adat, pt_diff_u_xx);
    // Calc_Error (&adat, pt_diff_u_yy);
    // Calc_Error (&adat, pt_diff_v_xx);
    // Calc_Error (&adat, pt_diff_v_yy);
    // 미분의 error를 계산해서 출력
    // Calc_Derivative_Error (&adat, pt);
    // 해를 출력
    FILE *sol_output = fopen(cdat.Output_Sol().c_str(), "w");
    // fprintf (sol_output, "%s\n", "# x, y, u, ux, uy, E-filed");
    for (int i = 0; i < adat.Pts_Num(); i++) {
        fprintf(sol_output, "%9d\t", i);
        fprintf(sol_output, "%23.16e\t", ptP[i].Time());
        fprintf(sol_output, "%23.16e\t", ptP[i].Coord()[0]);
        fprintf(sol_output, "%23.16e\t", ptP[i].Coord()[1]);
        fprintf(sol_output, "%23.16e\t", ptU[i].Value());
        fprintf(sol_output, "%23.16e\t", ptV[i].Value());
        fprintf(sol_output, "%23.16e\t", ptP[i].Value());
        fprintf(sol_output, "%23.16e\t", ptU[i].Phi()->Value());
        fprintf(sol_output, "%23.16e\t", ptV[i].Phi()->Value());
        fprintf(sol_output, "%23.16e\t", ptU[i].Diff('x')->Value());
        fprintf(sol_output, "%23.16e\t", ptU[i].Diff('y')->Value());
        fprintf(sol_output, "%23.16e\t", ptV[i].Diff('x')->Value());
        fprintf(sol_output, "%23.16e\t", ptV[i].Diff('y')->Value());
        fprintf(sol_output, "%23.16e\t", ptU[i].Diff('x')->Diff('x')->Value());
        fprintf(sol_output, "%23.16e\t", ptU[i].Diff('y')->Diff('y')->Value());
        fprintf(sol_output, "%23.16e\t", ptV[i].Diff('x')->Diff('x')->Value());
        fprintf(sol_output, "%23.16e\n", ptV[i].Diff('y')->Diff('y')->Value());
    }

    return 0;

// /* Interpolation points */
//
// // Interporation을 하기 위한 점의 개수
// int nInt = 10;
// // int nInt = 1;
// // Interporation을 하기 위한 점의 선언
// Point *pt_int = new Point[nInt];
// // x-축선위에서의 해의 표현식의 계수들의 struct 선언
// xData xdat;
// // y-축선위에서의 해의 표현식의 계수들의 struct 선언
// yData ydat;
// // 원의 각도 theta
// double theta = PI / 2.0E0;
// // 원의 각도에 더할 alpha
// double alpha = ZeroValue;
//
// // 점의 개수 만큼 밤복
// for (int i = 0; i < nInt; i++) {
//   // mark의 입력
//   pt_int[i].SetMark ('a');
//   // x-좌표의 입력
//   // pt_int[i].SetCoordinate ('x', 6.0E0 - 16.0E0 * 5.0E-4);
//   // pt_int[i].SetCoordinate ('x', 5.0E0 + i * 1.0E0 / nInt);
//   // y-좌표의 입력
//   // pt_int[i].SetCoordinate ('y', sqrt (6.0E0 * 6.0E0 - (pt_int[i].Coordinate ('x') - 6.0E0) * (pt_int[i].Coordinate ('x') - 6.0E0)) + 3.5E0);
//   // pt_int[i].SetCoordinate ('y', 9.5E0 - 1.0E-5);
//   // 원의 각도에 더할 alpha
//   alpha = 1.0E-4;
//   // 원의 각도 theta
//   theta = PI / 2.0E0 + (1.0E1 - double(i)) * alpha;
//   // x-좌표의 입력
//   pt_int[i].SetCoordinate ('x', 6.0E0 * cos (theta) + 6.0E0);
//   // pt_int[i].SetCoordinate ('x', 5.9755E0);
//   // y-좌표의 입력
//   pt_int[i].SetCoordinate ('y', 6.0E0 * sin (theta) + 3.5E0);
//   // pt_int[i].SetCoordinate ('y', 9.5E0 - 1.0E-5 * (1.0E1 - double(i)));
//   // 각 점의 index를 입력
//   pt_int[i].SetIndex (i);
//   // 각 점의 경계값의 입력
//   pt_int[i].SetBoundaryvalue (ZeroValue);
//   // 각 점의 경계조건의 입력
//   pt_int[i].SetCondition ('C');
//   // 각 점의 값의 초기화
//   pt_int[i].SetValue (ZeroValue);
//   // 각 점의 우변의 f값의 입력
//   pt_int[i].SetF (f_ftn (pt[i].Coordinate ('x'), pt[i].Coordinate ('y')));
//   // 각 점의 conductivity의 입력
//   // if (pt_int[i].Coordinate ('x') < 6.0E0 - sqrt (9.0E0 - (3.0E0 - 1.0E-4) * (3.0E0 - 1.0E-4))) pt_int[i].SetMaterialProperty (1.0E0);
//   // else                                                                                         pt_int[i].SetMaterialProperty (2.0E0);
//   pt_int[i].SetMaterialProperty (1.0E0);
//   // 각 점에 대한 현재 진행중인 점의 index를 알림
//   printf ("%s%lu%s%d", "\r", i + 1, "/", nInt);
// }
// printf ("\n");
//
// // 점의 개수 만큼 밤복
// for (int i = 0; i < nInt; i++) {
//   // 각점의 주의의 점들을 찾음
//   pt_int[i].FindAxialElementIntp (&adat, pt);
//   // 각 점에 대한 현재 진행중인 점의 index를 알림
//   printf ("%s%lu%s%d", "\r", i + 1, "/", nInt);
// }
// printf ("\n");
//
// // 점의 개수 만큼 밤복
// for (int i = 0; i < nInt; i++) {
//   // 각 점을 포함하는 국소축선의 시작점과 끝점의 좌표를 찾음
//   pt_int[i].SetMinMax ();
//   // 모든점에 대한 현재 진행중인 점의 index를 알림
//   printf ("%s%lu%s%d", "\r", i + 1, "/", nInt);
// }
// printf ("\n");
//
// // 점의 개수 만큼 밤복
// for (int i = 0; i < nInt; i++) {
//   // phi를 저장
//   pt_int[i].SetPhi (PhiInterp (&pt_int[i]));
//   // 값을 계산
//   // pt_int[i].SetValue (5.0E-1 * (CalcValue ('x', &pt_int[i], &xdat, &ydat) + CalcValue ('y', &pt_int[i], &xdat, &ydat)));
//   pt_int[i].SetValue (CalcValue ('y', &pt_int[i], &xdat, &ydat));
//   // x-성분의 미분을 계산
//   pt_int[i].SetDiff ('x', CalcDiff ('x', &pt_int[i], &xdat, &ydat));
//   // y-성분의 미분을 계산
//   pt_int[i].SetDiff ('y', CalcDiff ('y', &pt_int[i], &xdat, &ydat));
//   // 모든점에 대한 현재 진행중인 점의 index를 알림
//   printf ("%s%lu%s%d", "\r", i + 1, "/", nInt);
// }
// printf ("\n");
//
// FILE *sample_output = fopen ("SamplePoints", "w");
// for (int i = 0; i < 10; i++) {
//   fprintf (sample_output, "%9lu\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%23.16e\n", i, pt_int[i].Coordinate ('x'), pt_int[i].Coordinate ('y'),
//   pt_int[i].Value (), pt_int[i].Diff ('x'), pt_int[i].Diff ('y'), sqrt(pt_int[i].Diff ('x') * pt_int[i].Diff ('x') + pt_int[i].Diff ('y') * pt_int[i].Diff ('y')));
// }
// fclose (sample_output);
//
// // 경계조건을 출력
// // FILE *condition_output = fopen (cdat.Boundary ().c_str (), "w");
// // for (int i = 0; i < adat.Pts_Num (); i++) {
// //   fprintf (condition_output, "%9lu\t%c\n", i, pt[i].Condition ());
// // }
// // fclose (condition_output);
// return 0;
}

#endif
