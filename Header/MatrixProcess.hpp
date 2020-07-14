#ifndef MATRIXPROCESS_H
#define MATRIXPROCESS_H

#include "CalcRepresenSol.hpp"
#include "mkl_pardiso.h"

/* 생성자 */
MatrixProcess::MatrixProcess (AxialData *adat) {
  // 한 row에서 0이아닌 원소의 개수
  this->Int_num = 0;
  // 행렬에서 0이 아닌 원소의 개수
  this->ent_num = 0;
  // 오른쪽, 왼쪽, 위쪽, 아래쪽의 점들의 주소를 저장하기 위한 array
  this->arrInt = new Point*[26];
  // 오른쪽, 왼쪽, 위쪽, 아래쪽의 점에서의 값을 저장하기 위한 array
  this->arrEnt = new double[26];
  // arrInt에서 중복되는 점들을 제거한 index를 저장하기 위한 array
  this->uniqInt = new int[26];
  // arrInt에서 중복되는 점들을 제거한 점에서의 값을 저장하기 위한 array
  this->uniqEnt = new double[26];
  // uniqInt에서 크기순서로 정렬한 index를 저장하기 위한 array
  this->rowsInt = new int[26];
  // uniqEnt에서 크기순서로 정렬한 점에서의 값을 저장하기 위한 array
  this->rowsEnt = new double[26];
  // 크기순서로 정렬하기 위한 array
  this->sortInt = new int[26];
}

/* 소멸자 */
MatrixProcess::~MatrixProcess () {

}

/* class의 변수의 초기화 */
MatrixProcess & MatrixProcess::initialization (AxialData *adat, int k) {
  // 행렬의 크기 (k x k)
  this->matrixSize = k;
  // 행렬에서 0이 아닌 원소의 column의 index
  this->ja = new int[this->ent_num];
  // 행렬에서 각 row에서 처음으로 0이 아닌 원소의 ja에서의 index
  this->ia = new int[this->matrixSize + 1];
  // 행렬에서 0이 아닌 원소의 값
  this->ent = new double[this->ent_num];
  // 행렬식에서 오른쪽의 항
  this->rb = new double[this->matrixSize];
  // 행렬에서 0이 아닌 원소의 개수
  this->ent_num = 0;
  return *this;
}

/* 행렬의 0이 아닌 원소의 개수를 세는 모듈 */
MatrixProcess & MatrixProcess::countEnt_num (int node, AxialData *adat, Point *pt, xData *xdat, yData *ydat) {
  // 한 개의 row에서 0이 아닌 원소의 개수
  this->Int_num = 0;
  // 각 array의 초기화
  for (size_t i = 0; i < 26; i++) {
    // 오른쪽, 왼쪽, 위쪽, 아래쪽의 점들의 주소를 저장하기 위한 array
    this->uniqInt[i] = -1;
    // 오른쪽, 왼쪽, 위쪽, 아래쪽의 점에서의 값을 저장하기 위한 array
    this->uniqEnt[i] = 0.0;
    // arrInt에서 중복되는 점들을 제거한 index를 저장하기 위한 array
    this->arrInt[i] = NULL;
    // arrInt에서 중복되는 점들을 제거한 점에서의 값을 저장하기 위한 array
    this->arrEnt[i] = 0.0;
  }
  this->SettingAzimuth (pt, adat, node);
  // // arrInt의 원소가 경계점이면 초기화
  // for (size_t i = 0; i < 13; i++) {
  //   // arrInt의 원소가 존재하는 경우
  //   if (this->arrInt[i] != NULL) {
  //     // arrInt의 원소가 경계점인 경우
  //     if (this->arrInt[i]->Condition () == 'D') {
  //       // arrInt에 NULL값을 저장
  //       this->arrInt[i] = NULL;
  //     }
  //   }
  // }
  // 행렬이 2n x 2n 이기 때문에, 14 - 26번째 원소는 1 - 13번째 원소의 반복
  for (size_t i = 13; i < 26; i++) {
    // arrInt의 원소가 존재하는 경우
    if (this->arrInt[i - 13] != NULL) {
      // arrInt가 저장한 주소를 그대로 저장
      this->arrInt[i] = this->arrInt[i - 13];
    }
  }

  // x-축선에서의 식과 y-축선에서의 식의 합을 계산
  // 현재점에서의 u의 계수
  this->arrEnt[0 ] = xdat->Cu + ydat->Cu;
  // 오른쪽점에서의 u의 계수
  this->arrEnt[1 ] = xdat->Eu;
  // 왼쪽점에서의 u의 계수
  this->arrEnt[2 ] = xdat->Wu;
  // 위쪽점에서의 u의 계수
  this->arrEnt[3 ] = ydat->Nu;
  // 아래쪽점에서의 u의 계수
  this->arrEnt[4 ] = ydat->Su;
  // 오른쪽의 위쪽점에서의 u의 계수
  this->arrEnt[5 ] = xdat->ENu;
  // 오른쪽의 아래쪽점에서의 u의 계수
  this->arrEnt[6 ] = xdat->ESu;
  // 왼쪽의 위쪽점에서의 u의 계수
  this->arrEnt[7 ] = xdat->WNu;
  // 왼쪽의 아래쪽점에서의 u의 계수
  this->arrEnt[8 ] = xdat->WSu;
  // 위쪽의 오른쪽점에서의 u의 계수
  this->arrEnt[9 ] = ydat->NEu;
  // 위쪽의 왼쪽점에서의 u의 계수
  this->arrEnt[10] = ydat->NWu;
  // 아래쪽의 오른쪽점에서의 u의 계수
  this->arrEnt[11] = ydat->SEu;
  // 아래쪽의 왼쪽점에서의 u의 계수
  this->arrEnt[12] = ydat->SWu;
  // 현재점에서의 phi의 계수
  this->arrEnt[13] = xdat->Cphi + ydat->Cphi;
  // 오른쪽점에서의 phi의 계수
  this->arrEnt[14] = xdat->Ephi;
  // 왼쪽점에서의 phi의 계수
  this->arrEnt[15] = xdat->Wphi;
  // 위쪽점에서의 phi의 계수
  this->arrEnt[16] = ydat->Nphi;
  // 아래쪽점에서의 phi의 계수
  this->arrEnt[17] = ydat->Sphi;
  // 오른쪽의 위쪽점에서의 phi의 계수
  this->arrEnt[18] = xdat->ENphi;
  // 오른쪽의 아래쪽점에서의 phi의 계수
  this->arrEnt[19] = xdat->ESphi;
  // 왼쪽의 위쪽점에서의 phi의 계수
  this->arrEnt[20] = xdat->WNphi;
  // 왼쪽의 아래쪽점에서의 phi의 계수
  this->arrEnt[21] = xdat->WSphi;
  // 위쪽의 오른쪽점에서의 phi의 계수
  this->arrEnt[22] = ydat->NEphi;
  // 위쪽의 왼쪽점에서의 phi의 계수
  this->arrEnt[23] = ydat->NWphi;
  // 아래쪽의 오른쪽점에서의 phi의 계수
  this->arrEnt[24] = ydat->SEphi;
  // 아래쪽의 왼쪽점에서의 phi의 계수
  this->arrEnt[25] = ydat->SWphi;
  // arrEnt의 원소가 0인 경우 arrInt의 원소의 값을 초기화
  for (size_t i = 0; i < 26; i++) {
    // arrEnt의 원소가 0인 경우
    if (IsEqualDouble (this->arrEnt[i], ZeroValue)) {
      // arrInt의 원소의 값을 초기화
      this->arrInt[i] = NULL;
    }
  }
  // arrInt의 주소에 해당하는 index의 중복된 값을 제거하기 위한 변수선언 (j: 중복되는 점의 arrInt에서의 위치, k: 중복되지 않는 원소의 개수)
  int j = -1, k = 0;
  // arrInt의 주소에 해당하는 index의 중복된 값을 제거
  for (size_t i = 0; i < 26; i++) {
    // arrInt의 원소가 존재하는 경우
    if (this->arrInt[i] != NULL) {
      // x-축선에서의 식과 y-축선에서의 식의 합인 경우
      if (i < 13) {
        // arrInt의 원소가 중복되는지의 여부와 중복된다면, 중복되는 점의 arrInt의 원소의 내부점의 index를 저장 (중복x: -1)
        j = is_member (adat->PtsTOPts ('P', this->arrInt[i]->Index ()), this->uniqInt, k);
        // arrInt의 원소가 중복되지 않는 경우
        if (j == -1) {
          // uniqInt의 원소에 arrInt의 원소의 내부점의 index를 저장
          this->uniqInt[k] = adat->PtsTOPts ('P', this->arrInt[i]->Index ());
          // 중복되지 않는 원소의 개수에 1을 더한다
          k += 1;
        }
      }
      // x-축선에서의 식과 y-축선에서의 식의 차인 경우
      else {
        // arrInt의 원소가 중복되는지의 여부와 중복된다면, 중복되는 점의 arrInt의 원소의 phi의 값을 계산하는 점의 index를 저장 (중복x: -1)
        j = is_member (adat->PtsTOPts ('T', this->arrInt[i]->Index ()) + adat->In_Pts_Num (), this->uniqInt, k);
        // arrInt의 원소가 중복되지 않는 경우
        if (j == -1) {
          // uniqInt의 원소에 arrInt의 원소의 내부점의 index를 저장
          this->uniqInt[k] = adat->PtsTOPts ('T', this->arrInt[i]->Index ()) + adat->In_Pts_Num ();
          // 중복되지 않는 원소의 개수에 1을 더한다
          k += 1;
        }
      }
    }
  }
  // 중복되지 않는 원소의 개수 count
  for (size_t i = 0; i < 26; i++) {
    // uniqInt의 원소가 존재한다면 (중복되지 않는 원소라면) Int_num에 1을 더한다
    if (this->uniqInt[i] > -1) this->Int_num += 1;
  }
  // 행렬에서 0이 아닌 원소의 개수에 Int_num을 더한다
  this->ent_num += this->Int_num;
  return *this;
}

/* 행렬을 구성하는 모듈 */
MatrixProcess & MatrixProcess::MakeMatrixSystem (int node, AxialData *adat, Point *pt, xData *xdat, yData *ydat) {
  // 한 개의 row에서 0이 아닌 원소의 개수
  this->Int_num = 0;
  // 각 array의 초기화
  for (size_t i = 0; i < 26; i++) {
    // 오른쪽, 왼쪽, 위쪽, 아래쪽의 점들의 주소를 저장하기 위한 array
    this->arrInt[i] = NULL;
    // 오른쪽, 왼쪽, 위쪽, 아래쪽의 점에서의 값을 저장하기 위한 array
    this->arrEnt[i] = ZeroValue;
    // arrInt에서 중복되는 점들을 제거한 index를 저장하기 위한 array
    this->uniqInt[i] = -1;
    // arrInt에서 중복되는 점들을 제거한 점에서의 값을 저장하기 위한 array
    this->uniqEnt[i] = ZeroValue;
    // uniqInt에서 크기순서로 정렬한 index를 저장하기 위한 array
    this->rowsInt[i] = -1;
    // uniqEnt에서 크기순서로 정렬한 점에서의 값을 저장하기 위한 array
    this->rowsEnt[i] = ZeroValue;
    // 크기순서로 정렬하기 위한 array
    this->sortInt[i] = -1;
  }
  this->SettingAzimuth (pt, adat, node);
  // // arrInt의 원소가 경계점이면 초기화
  // for (size_t i = 0; i < 13; i++) {
  //   // arrInt의 원소가 존재하는 경우
  //   if (this->arrInt[i] != NULL) {
  //     // arrInt의 원소가 경계점인 경우
  //     if (this->arrInt[i]->Condition () == 'D') {
  //       // arrInt에 NULL값을 저장
  //       this->arrInt[i] = NULL;
  //     }
  //   }
  // }
  // 행렬이 2n x 2n 이기 때문에, 14 - 26번째 원소는 1 - 13번째 원소의 반복
  for (size_t i = 13; i < 26; i++) {
    // arrInt의 원소가 존재하는 경우
    if (this->arrInt[i - 13] != NULL) {
      // arrInt가 저장한 주소를 그대로 저장
      this->arrInt[i] = this->arrInt[i - 13];
    }
  }
  // x-축선에서의 식과 y-축선에서의 식의 합을 계산
  // 현재점에서의 u의 계수
  this->arrEnt[0 ] = xdat->Cu + ydat->Cu;
  // 오른쪽점에서의 u의 계수
  this->arrEnt[1 ] = xdat->Eu;
  // 왼쪽점에서의 u의 계수
  this->arrEnt[2 ] = xdat->Wu;
  // 위쪽점에서의 u의 계수
  this->arrEnt[3 ] = ydat->Nu;
  // 아래쪽점에서의 u의 계수
  this->arrEnt[4 ] = ydat->Su;
  // 오른쪽의 위쪽점에서의 u의 계수
  this->arrEnt[5 ] = xdat->ENu;
  // 오른쪽의 아래쪽점에서의 u의 계수
  this->arrEnt[6 ] = xdat->ESu;
  // 왼쪽의 위쪽점에서의 u의 계수
  this->arrEnt[7 ] = xdat->WNu;
  // 왼쪽의 아래쪽점에서의 u의 계수
  this->arrEnt[8 ] = xdat->WSu;
  // 위쪽의 오른쪽점에서의 u의 계수
  this->arrEnt[9 ] = ydat->NEu;
  // 위쪽의 왼쪽점에서의 u의 계수
  this->arrEnt[10] = ydat->NWu;
  // 아래쪽의 오른쪽점에서의 u의 계수
  this->arrEnt[11] = ydat->SEu;
  // 아래쪽의 왼쪽점에서의 u의 계수
  this->arrEnt[12] = ydat->SWu;
  // 현재점에서의 phi의 계수
  this->arrEnt[13] = xdat->Cphi + ydat->Cphi;
  // 오른쪽점에서의 phi의 계수
  this->arrEnt[14] = xdat->Ephi;
  // 왼쪽점에서의 phi의 계수
  this->arrEnt[15] = xdat->Wphi;
  // 위쪽점에서의 phi의 계수
  this->arrEnt[16] = ydat->Nphi;
  // 아래쪽점에서의 phi의 계수
  this->arrEnt[17] = ydat->Sphi;
  // 오른쪽의 위쪽점에서의 phi의 계수
  this->arrEnt[18] = xdat->ENphi;
  // 오른쪽의 아래쪽점에서의 phi의 계수
  this->arrEnt[19] = xdat->ESphi;
  // 왼쪽의 위쪽점에서의 phi의 계수
  this->arrEnt[20] = xdat->WNphi;
  // 왼쪽의 아래쪽점에서의 phi의 계수
  this->arrEnt[21] = xdat->WSphi;
  // 위쪽의 오른쪽점에서의 phi의 계수
  this->arrEnt[22] = ydat->NEphi;
  // 위쪽의 왼쪽점에서의 phi의 계수
  this->arrEnt[23] = ydat->NWphi;
  // 아래쪽의 오른쪽점에서의 phi의 계수
  this->arrEnt[24] = ydat->SEphi;
  // 아래쪽의 왼쪽점에서의 phi의 계수
  this->arrEnt[25] = ydat->SWphi;
  // arrEnt의 원소가 0인 경우 arrInt의 원소의 값을 초기화
  for (size_t i = 0; i < 26; i++) {
    // arrEnt의 원소가 0인 경우
    if (IsEqualDouble (this->arrEnt[i], ZeroValue)) {
      // arrInt의 원소의 값을 초기화
      this->arrInt[i] = NULL;
    }
  }
  // for (size_t i = 0; i < 26; i++) {
  //   if (this->arrInt[i] == NULL) continue;
  //   printf ("this->arrInt[%lu] = ", i);
  //   printf ("%d, ", this->arrInt[i]->Index ());
  // }
  // printf ("\n");

  // double temp_matrix_confirm = 0;
  //
  //
  // temp_matrix_confirm = 0;
  // for (size_t i = 0; i < 13; i++) {
  //   if (this->arrInt[i]) temp_matrix_confirm += this->arrEnt[i] * u_ftn (this->arrInt[i]->Coord ().Value ('x'), this->arrInt[i]->Coord ().Value ('y'));
  //   if (this->arrInt[i]) temp_matrix_confirm += this->arrEnt[i + 13] * phi_ftn (this->arrInt[i]->Coord ().Value ('x'), this->arrInt[i]->Coord ().Value ('y'));
  // }

  // if (pt->Index () == 0) {
  //   this->printPointDataError (node, pt);
  //   for (size_t i = 0; i < 13; i++) {
  //     if (this->arrInt[i] != NULL) printf ("arrEnt[%2zu] = %23.15e\n", i,  this->arrEnt[i] * u_ftn (this->arrInt[i]->Coordinate ('x'), this->arrInt[i]->Coordinate ('y')));
  //     if (this->arrInt[i] != NULL) printf ("arrEnt[%2zu] = %23.15e\n", i + 13,  this->arrEnt[i + 13] * phi_ftn (this->arrInt[i]->Coordinate ('x'), this->arrInt[i]->Coordinate ('y')));
  //   }
  //   printf ("%f\n", (xdat->F + ydat->F));
  //   printf ("%f\n", (xdat->F));
  //   printf ("%f\n", (ydat->F));
  //   printf ("%23.16e\n", xdat->Cphi * phi_ftn (this->arrInt[0]->Coordinate ('x'), this->arrInt[0]->Coordinate ('y')));
  //   printf ("%23.16e\n", xdat->Cphi);
  //   printf ("%23.16e\n", phi_ftn (this->arrInt[0]->Coordinate ('x'), this->arrInt[0]->Coordinate ('y')));
  //   exit (1);
  // }

  // // node가 내부점의 개수보다 작은 경우
  // if (node < adat->In_Pts_Num ()) {
  //   // x-축선에서의 식과 y-축선에서의 우변의 합을 계산
  //   temp_matrix_confirm -= (xdat->F + ydat->F);
  // }
  // // node가 내부점의 개수보다 큰 경우
  // else {
  //   // x-축선에서의 식과 y-축선에서의 우변의 차을 계산
  //   temp_matrix_confirm -= (xdat->F - ydat->F);
  // }
  //
  // if (temp_matrix_confirm > 1.0E-2) printf ("xdat->F = %23.16e, ydat->F = %23.16e\n", xdat->F, ydat->F), this->printPointDataError (node, pt), printf ("\nmatrix_confirm = %23.16e\n", temp_matrix_confirm), pt->PrintDebuggingData ("11111", xdat, ydat), exit (111);

  // FILE *matrix_confirm;
  // matrix_confirm = fopen ("matrix_confirm.dat","a");
  // if (matrix_confirm!=NULL) {
  //   fprintf (matrix_confirm, "%9d\t%23.16e\n", node, temp_matrix_confirm);
  //   fclose (matrix_confirm);
  // }
  // arrInt의 주소에 해당하는 index의 중복된 값을 제거하기 위한 변수선언 (j: 중복되는 점의 arrInt에서의 위치, k: 중복되지 않는 원소의 개수)
  int j = -1, k = 0;
  // arrInt의 주소에 해당하는 index의 중복된 값을 제거하고 값을 더함
  for (size_t i = 0; i < 26; i++) {
    // arrInt의 원소가 존재하는 경우
    if (this->arrInt[i] != NULL) {
      // x-축선에서의 식과 y-축선에서의 식의 합인 경우
      if (i < 13) {
        // arrInt의 원소가 중복되는지의 여부와 중복된다면, 중복되는 점의 arrInt의 원소의 내부점의 index를 저장 (중복x: -1)
        j = is_member (adat->PtsTOPts ('P', this->arrInt[i]->Index ()), this->uniqInt, k);
        // arrInt의 원소가 중복되지 않는 경우
        if (j == -1) {
          // uniqInt의 원소에 arrInt의 원소의 내부점의 index를 저장
          this->uniqInt[k] = adat->PtsTOPts ('P', this->arrInt[i]->Index ());
          // uniqEnt의 원소에 arrEnt의 원소를 저장
          this->uniqEnt[k] = this->arrEnt[i];
          // 중복되지 않는 원소의 개수에 1을 더한다
          k += 1;
        }
        // arrInt의 원소가 중복되는 경우
        else {
          // 중복되는 점의 위치에 해당하는 uniqEnt에 arrEnt원소를 더함
          this->uniqEnt[j] += this->arrEnt[i];
        }
      }
      // x-축선에서의 식과 y-축선에서의 식의 차인 경우
      else {
        // arrInt의 원소가 중복되는지의 여부와 중복된다면, 중복되는 점의 arrInt의 원소의 phi의 값을 계산하는 점의 index를 저장 (중복x: -1)
        j = is_member (adat->PtsTOPts ('T', this->arrInt[i]->Index ()) + adat->In_Pts_Num (), this->uniqInt, k);
        // arrInt의 원소가 중복되지 않는 경우
        if (j == -1) {
          // uniqInt의 원소에 arrInt의 원소의 내부점의 index를 저장
          this->uniqInt[k] = adat->PtsTOPts ('T', this->arrInt[i]->Index ()) + adat->In_Pts_Num ();
          // uniqEnt의 원소에 arrEnt의 원소를 저장
          this->uniqEnt[k] = this->arrEnt[i];
          // 중복되지 않는 원소의 개수에 1을 더한다
          k += 1;
        }
        // arrInt의 원소가 중복되는 경우
        else {
          // 중복되는 점의 위치에 해당하는 uniqEnt에 arrEnt원소를 더함
          this->uniqEnt[j] += this->arrEnt[i];
        }
      }
    }
  }
  // 중복되지 않는 원소의 개수 count
  for (size_t i = 0; i < 26; i++) {
    // uniqInt의 원소가 존재한다면 (중복되지 않는 원소라면) Int_num에 1을 더한다
    if (this->uniqInt[i] > -1) this->Int_num += 1;
  }
  // uniqInt의 원소를 크기순서로 정렬
  for (size_t i = 0; i < 26; i++) {
    // uniqInt의 원소가 -1인 경우 다음으로 넘어감
    if (this->uniqInt[i] == -1) continue;
    // uniqInt의 원소의 크기를 비교
    for (size_t j = 0; j < 26; j++) {
      // uniqInt의 원소가 -1인 경우, 종료
      if (this->uniqInt[j] == -1) break;
      // uniqInt의 i번째 원소가 j번째 원소보다 크거나 같은 경우, sortInt에 1을 더한다.
      if (this->uniqInt[i] >= this->uniqInt[j]) this->sortInt[i] += 1;
    }
  }
  // rowsInt와 rowsEnt에 중복을 제거하고 크기순서대로 정렬을 한 row에서의 index
  for (size_t i = 0; i < 26; i++) {
    // sortInt의 원소가 -1이 아니라면
    if (this->sortInt[i] != -1) {
      // sortInt의 원소에 해당하는 rowsInt의 위치에 uniqInt의 원소를 저장
      this->rowsInt[this->sortInt[i]] = this->uniqInt[i];
      // sortInt의 원소에 해당하는 rowsEnt의 위치에 uniqEnt의 원소를 저장
      this->rowsEnt[this->sortInt[i]] = this->uniqEnt[i];
    }
  }
  // ia에 지금까지의 0이 아닌 원소의 개수를 저장
  this->ia[node] = this->ent_num;
  this->rb[node] = xdat->F + ydat->F;
  if (this->rb[node] != this->rb[node]) printf ("xdat->F = %23.16e\n", xdat->F), printf ("ydat->F = %23.16e\n", ydat->F), this->printPointDataError (node, pt), exit (123);
  // ja의 array의 index에 지금까지의 0이 아닌 원소의 개수를 저장
  this->ja_num = this->ent_num;
  // ja와 ent array의 원소를 저장
  for (size_t i = 0; i < 26; i++) {
    // rowsInt의 원소가 -1보다 큰 경우 (점의 index가 존재하는경우)
    if (rowsInt[i] > -1) {
      // ja에 rowsInt의 원소를 저장
      this->ja[this->ja_num] = this->rowsInt[i];
      // entdp rowsEnt의 원소를 저장
      this->ent[this->ja_num] = this->rowsEnt[i];
      // ent에 저장된 값이 NaN값인 경우
      if (this->ent[this->ja_num] != this->ent[this->ja_num]) {
        // 해당하는 점의 정보를 출력하고 종료
        this->printPointDataError (node, pt);
      }
      // ja의 array의 index에 1을 더함
      this->ja_num += 1;
    }
  }
  // 행렬에서 0이 아닌 원소의 개수에 Int_num을 더한다
  this->ent_num += this->Int_num;
  return *this;
}

/* 행렬식을 계산하는 모듈 */
MatrixProcess & MatrixProcess::calcMatrix (ControlData *cdat, AxialData *adat, Point *pts) {
  // 행렬이 크기
  int      n = matrixSize;
  // 우변의 개수
  int      nrhs = 1;
  void    *pt[64];
  int      iparm[64];
  int      mtype, maxfct, mnum, phase, error = 0, msglvl;
  int      idum;
  // 해
  double   x[matrixSize];

  this->ia[this->matrixSize] = this->ent_num;

  for (size_t i = 0; i < this->matrixSize + 1; i++) this->ia[i] += 1;
  for (size_t i = 0; i < this->ent_num; i++) this->ja[i] += 1;
  for (size_t i = 0; i < this->matrixSize; i++) x[i] = 0.0;
  for (size_t i = 0; i < 64; i++) iparm[i] = 0;
  iparm[0] = 1;         /* No solver default */
  iparm[1] = 3;         /* Fill-in reordering from METIS */
  iparm[3] = 0;         /* No iterative-direct algorithm */
  iparm[4] = 0;         /* No user fill-in reducing permutation */
  iparm[5] = 0;         /* Write solution into x */
  iparm[6] = 0;         /* Not in use */
  iparm[7] = 2;         /* Max numbers of iterative refinement steps */
  iparm[8] = 0;         /* Not in use */
  iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;        /* Conjugate transposed/transpose solve */
  iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
  iparm[13] = 0;        /* Output: Number of perturbed pivots */
  iparm[14] = 0;        /* Not in use */
  iparm[15] = 0;        /* Not in use */
  iparm[16] = 0;        /* Not in use */
  iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1;       /* Output: Mflops for LU factorization */
  iparm[19] = 0;        /* Output: Numbers of CG Iterations */
  maxfct = 1;           /* Maximum number of numerical factorizations. */
  mnum = 1;             /* Which factorization to use. */
  msglvl = 0;           /* Print statistical information in file */
  error = 0;            /* Initialize error flag */
  mtype = 11;

  iparm[28] = 0;


  for (size_t i = 0; i < 64; i++) pt[i] = 0;

  phase = 13;
  pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, this->ent, this->ia, this->ja, &idum, &nrhs, iparm, &msglvl, this->rb, x, &error);
  if ( error != 0 ) {
    printf ("\nERROR during solution: %d", error);
    exit (3);
  }

  for (size_t i = 0; i < matrixSize; i++) {
    // printf ("i = %zu, x[%zu] = %23.16e\n", i, i, x[i]);
    if (i < adat->In_Pts_Num ()) pts[adat->PtsTOPts ('I', i)].SetValue (x[i]);
    else                         pts[adat->PtsTOPts ('H', i % adat->In_Pts_Num ())].Phi ()->SetValue (x[i]);

    // if (i > adat->In_Pts_Num ()) if (pts[adat->PtsTOPts ('H', i % adat->In_Pts_Num ())] == 'D') printf ("index = %5d, numeric value = %23.16e, analytic value = %23.16e, error = %23.16e\n", pts[adat->PtsTOPts ('H', i % adat->In_Pts_Num ())].Index (), x[i], u_ftn (pts[adat->PtsTOPts ('H', i % adat->In_Pts_Num ())].Phi ()), fabs (x[i] - u_ftn (pts[adat->PtsTOPts ('H', i % adat->In_Pts_Num ())].Phi ())));
  }

  // FILE *Solution;
  // Solution = fopen ("Solution.dat","w");
  // if (Solution!=NULL)
  // {
  //   for (size_t i = 0; i < matrixSize; i++) {
  //     fprintf(Solution, "%10lu\t%23.16e\n", i, x[i]);
  //   }
  //   fclose (Solution);
  // }

  phase = -1;
  pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, this->ent, this->ia, this->ja, &idum, &nrhs, iparm, &msglvl, this->rb, x, &error );
  ent_num = 0;
  delete [] this->ia;
  delete [] this->ja;
  delete [] this->ent;
  delete [] this->rb;

  return *this;
}
/* 행렬의 정보를 내보내는 모듈 */
MatrixProcess & MatrixProcess::ExportMatrixData (AxialData *adat) {
  // ia의 마지막 원소는 행렬에서 0이 아닌 원소의 개수
  ia[matrixSize] = ent_num;
  // ent를 내보내기
  FILE *ent_output = fopen ("ent.dat", "w");
  for (size_t i = 0; i < this->ent_num; i++) {
    fprintf (ent_output," %9lu\t%23.16e\n", i, this->ent[i]);
  }
  fclose (ent_output);
  // ja를 내보내기
  FILE *ja_output = fopen ("ja.dat", "w");
  for (size_t i = 0; i < this->ent_num; i++) {
    fprintf (ja_output, "%9lu\t%9d\n", i, this->ja[i]);
  }
  fclose (ja_output);
  // ia를 내보내기
  FILE *ia_output = fopen ("ia.dat", "w");
  for (size_t i = 0; i < this->matrixSize + 1; i++) {
    fprintf (ia_output, "%9lu\t%9d\n", i, this->ia[i]);
  }
  fclose (ia_output);
  // rb를 내보내기
  FILE *rb_output = fopen ("rb.dat", "w");
  for (size_t i = 0; i < this->matrixSize; i++) {
    fprintf (rb_output, "%9lu\t%23.16e\n", i, this->rb[i]);
  }
  fclose (rb_output);
  return *this;
}

MatrixProcess & MatrixProcess::SettingAzimuth (Point *pt, AxialData *adat, const int node) {
  // 현재점의 index
  this->arrInt[0] = pt;
  // 오른쪽점의 index
  if (pt->EWNS ('E', 'E')) this->arrInt[1] = pt->EWNS ('E', 'E');
  // 왼쪽점의 index
  if (pt->EWNS ('W', 'W')) this->arrInt[2] = pt->EWNS ('W', 'W');
  // 위쪽점의 index
  if (pt->EWNS ('N', 'N')) this->arrInt[3] = pt->EWNS ('N', 'N');
  // 아래쪽점의 index
  if (pt->EWNS ('S', 'S')) this->arrInt[4] = pt->EWNS ('S', 'S');
  // 오른쪽의 위쪽점의 index
  if (pt->EWNS ('E', 'N')) this->arrInt[5] = pt->EWNS ('E', 'N');
  // 오른쪽의 아래쪽점의 index
  if (pt->EWNS ('E', 'S')) this->arrInt[6] = pt->EWNS ('E', 'S');
  // 왼쪽의 위쪽점의 index
  if (pt->EWNS ('W', 'N')) this->arrInt[7] = pt->EWNS ('W', 'N');
  // 왼쪽의 아래쪽점의 index
  if (pt->EWNS ('W', 'S')) this->arrInt[8] = pt->EWNS ('W', 'S');
  // 위쪽의 오른쪽점의 index
  if (pt->EWNS ('N', 'E')) this->arrInt[9] = pt->EWNS ('N', 'E');
  // 위쪽의 왼쪽점의 index
  if (pt->EWNS ('N', 'W')) this->arrInt[10] = pt->EWNS ('N', 'W');
  // 아래쪽의 오른쪽점의 index
  if (pt->EWNS ('S', 'E')) this->arrInt[11] = pt->EWNS ('S', 'E');
  // 아래쪽의 왼쪽점의 index
  if (pt->EWNS ('S', 'W')) this->arrInt[12] = pt->EWNS ('S', 'W');

  if (pt->Condition () == 'N' && node >= adat->In_Pts_Num ()) {
    double XM = pt->Pressure ()->MinMaxCoordinate ('x', 'm'), XP = pt->Pressure ()->MinMaxCoordinate ('x', 'p');
    double YM = pt->Pressure ()->MinMaxCoordinate ('y', 'm'), YP = pt->Pressure ()->MinMaxCoordinate ('y', 'p');

    unordered_map<char, double*> m;
    unordered_map<char, char> opposite_azimuth;
    unordered_map<char, unordered_map<char, int>> array_int;
    unordered_map<char, int> array_int_e, array_int_w, array_int_n, array_int_s;
    unordered_map<char, char(*)[3]> sa;
    char sae[3] = {'E', 'N', 'S'}, saw[3] = {'W', 'N', 'S'}, san[3] = {'N', 'E', 'W'}, sas[3] = {'S', 'E', 'W'};
    unordered_map<char, string> coord;

    opposite_azimuth['E'] = 'W', opposite_azimuth['W'] = 'E', opposite_azimuth['N'] = 'S', opposite_azimuth['S'] = 'N';
    array_int['E'] = array_int_e, array_int['W'] = array_int_w, array_int['N'] = array_int_n, array_int['S'] = array_int_s;
    array_int['E']['E'] = 1, array_int['E']['N'] =  5, array_int['E']['S'] =  6;
    array_int['W']['W'] = 2, array_int['W']['N'] =  7, array_int['W']['S'] =  8;
    array_int['N']['N'] = 3, array_int['N']['E'] =  9, array_int['N']['W'] = 10;
    array_int['S']['S'] = 4, array_int['S']['E'] = 11, array_int['S']['W'] = 12;
    m['N'] = &YP, m['S'] = &YM;
    coord['N'] = "yp", coord['S'] = "ym";
    sa['E'] = &sae, sa['W'] = &saw, sa['N'] = &san, sa['S'] = &sas;

    bool sh = true;

    for (const auto &i : {'N', 'S'})
    if (pt->IsBoundary (opposite_azimuth[i]))
    for (const auto &j : *sa[i])
    if (pt->EWNS2nd (i, j))
    *m[i] = pt->EWNS2nd (i, j)->Coord ().Value ('y'),
    this->arrInt[array_int[i][j]] = pt->EWNS2nd (i, j),
    sh = false;


    if (IsEqualDouble (XP - XM, YP - YM) || sh)
    for (const auto &i : {'E', 'W'})
    if (pt->IsBoundary (opposite_azimuth[i]))
    for (const auto &j : *sa[i])
    if (pt->EWNS2nd (i, j))
    this->arrInt[array_int[i][j]] = pt->EWNS2nd (i, j);
  }

  if (pt->Condition () == 'D' && node >= adat->In_Pts_Num ()) {
    char errorMassage[256];
    Point *calc_pt = NULL;
    unordered_map<char, int> array_int;
    unordered_map<char, char> opposite_azimuth;
    unordered_map<char, char> coord;

    coord['E'] = 'x', coord['W'] = 'x', coord['N'] = 'y', coord['S'] = 'y';
    array_int['E'] = 2, array_int['W'] = 1, array_int['N'] = 4, array_int['S'] = 3;
    opposite_azimuth['E'] = 'W', opposite_azimuth['W'] = 'E', opposite_azimuth['N'] = 'S', opposite_azimuth['S'] = 'N';

    for (const auto &i : {'E', 'W', 'N', 'S'})
    if (pt->Axis (coord[i]) > -1)
    if (pt->IsBoundary (opposite_azimuth[i]))
    if (pt->EWNS (i, i))
    if (pt->EWNS (i, i)->EWNS (i, i))
    if (pt->EWNS (i, i)->EWNS (i, i)->Condition () == 'C') {calc_pt = pt->EWNS (i, i)->EWNS (i, i); break;}

    if (!calc_pt) sprintf (errorMassage, "SettingAzimuth with Dirichlet boundary, calc_pt = NULL point"), pt->PrintDebuggingData ("11110", "SettingAzimuth"), PrintError (errorMassage);

    SettingAzimuth (calc_pt, adat, node);

    for (const auto &i : {'E', 'W', 'N', 'S'})
    if (pt->Axis (coord[i]) > -1)
    if (pt->IsBoundary (opposite_azimuth[i]))
    if (pt->EWNS (i, i))
    if (pt->EWNS (i, i)->EWNS (i, i))
    if (pt->EWNS (i, i)->EWNS (i, i)->Condition () == 'C') {this->arrInt[array_int[i]] = pt; break;}
  }

  return *this;
}

bool is_error_tol (double tol, AxialData *adat, Point *pt, int *k) {
  //
  //   double error_u = pt[0].ReturnError ('u'), error_v = pt[0].ReturnError ('v'), error_p = pt[0].ReturnError ('p');
  //   double value_u = fabs (pt[0].ReturnValue ('u')), value_v = fabs (pt[0].ReturnValue ('v')), value_p = fabs (pt[0].ReturnValue ('p'));
  //   double error = 0.0;
  //
  //   pt[0].SetPreviousvalue ('u', pt[0].ReturnValue ('u'));
  //   pt[0].SetPreviousvalue ('v', pt[0].ReturnValue ('v'));
  //   pt[0].SetPreviousvalue ('p', pt[0].ReturnValue ('p'));
  //
  //   for (size_t i = 1; i < adat->ReturnPts_num (); i++) {
  //
  //     if (error_u < pt[i].ReturnError ('u')) error_u = pt[i].ReturnError ('u');
  //     if (error_v < pt[i].ReturnError ('v')) error_v = pt[i].ReturnError ('v');
  //     if (error_p < pt[i].ReturnError ('p')) error_p = pt[i].ReturnError ('p');
  //
  //     if (value_u < fabs (pt[i].ReturnValue ('u'))) value_u = fabs (pt[i].ReturnValue ('u'));
  //     if (value_v < fabs (pt[i].ReturnValue ('v'))) value_v = fabs (pt[i].ReturnValue ('v'));
  //     if (value_p < fabs (pt[i].ReturnValue ('p'))) value_p = fabs (pt[i].ReturnValue ('p'));
  //
  //     pt[i].SetPreviousvalue ('u', pt[i].ReturnValue ('u'));
  //     pt[i].SetPreviousvalue ('v', pt[i].ReturnValue ('v'));
  //     pt[i].SetPreviousvalue ('p', pt[i].ReturnValue ('p'));
  //
  //   }
  //
  //   error = error_u / value_u;
  //
  //   if (error < error_v / value_v) error = error_v / value_v;
  //   if (error < error_p / value_p) error = error_p / value_p;
  //
  //   *k += 1;
  //
  //   printf ("============ k = %06d ============\n", *k);
  //   printf ("Error of u = %23.16e\n", error_u / value_u);
  //   printf ("Error of v = %23.16e\n", error_v / value_v);
  //   printf ("Error of p = %23.16e\n", error_p / value_p);
  //   printf ("Error      = %23.16e\n", error);
  //   printf ("====================================\n\n");
  //
  //   if (error < tol) return false;
  return true;

}
// ent에 저장된 값이 NaN값인 경우 점의 정보를 출력
MatrixProcess & MatrixProcess::printPointDataError (int node, Point * pt) {
  // ent에 NaN값이 저장되었음을 알림
  printf ("\n\n%s\n", "NaN value stored in ent array");
  // 점의 index를 출력
  printf ("\n\n%s%d\n", "pt->Index () = ", pt->Index ());
  // 점의 경계조건을 출력
  printf ("%s%c\n\n", "pt->Condition () = ", pt->Condition ());
  // 점의 좌표를 출력
  printf ("%s ( %23.16e, %23.16e)\n", "pt->Coordinate = ", pt->Coord ().Value ('x'), pt->Coord ().Value ('y'));
  // 점을 포함하는 국소 x-축선의 왼쪽끝의 좌표와 오른쪽끝의 좌표를 출력
  printf ("xm = %23.16e, xp = %23.16e\n", pt->MinMaxCoordinate ('x', 'm'), pt->MinMaxCoordinate ('x', 'p'));
  // 점을 포함하는 국소 y-축선의 아래쪽끝의 좌표와 위쪽끝의 좌표를 출력
  printf ("ym = %23.16e, yp = %23.16e\n\n", pt->MinMaxCoordinate ('y', 'm'), pt->MinMaxCoordinate ('y', 'p'));
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
  // 점의 오른쪽점이 존재하는 경우, 점의 오른쪽점의 index를 출력
  if (pt->EWNS ('E', 'E') != NULL) printf ("%-6s%d\n", "E = ", pt->EWNS ('E', 'E')->Index ());
  // 점의 왼쪽점이 존재하는 경우, 점의 왼쪽점의 index를 출력
  if (pt->EWNS ('W', 'W') != NULL) printf ("%-6s%d\n", "W = ", pt->EWNS ('W', 'W')->Index ());
  // 점의 위쪽점이 존재하는 경우, 점의 위쪽점의 index를 출력
  if (pt->EWNS ('N', 'N') != NULL) printf ("%-6s%d\n", "N = ", pt->EWNS ('N', 'N')->Index ());
  // 점의 아래쪽점이 존재하는 경우, 점의 아래쪽점의 index를 출력
  if (pt->EWNS ('S', 'S') != NULL) printf ("%-6s%d\n", "S = ", pt->EWNS ('S', 'S')->Index ());
  // 점의 오른쪽의 위쪽점이 존재하는 경우, 점의 오른쪽의 위쪽점의 index를 출력
  if (pt->EWNS ('E', 'N') != NULL) printf ("%-6s%d\n", "EN = ", pt->EWNS ('E', 'N')->Index ());
  // 점의 오른쪽의 아래쪽점이 존재하는 경우, 점의 오른쪽의 아래쪽점의 index를 출력
  if (pt->EWNS ('E', 'S') != NULL) printf ("%-6s%d\n", "ES = ", pt->EWNS ('E', 'S')->Index ());
  // 점의 왼쪽의 위쪽점이 존재하는 경우, 점의 왼쪽의 위쪽점의 index를 출력
  if (pt->EWNS ('W', 'N') != NULL) printf ("%-6s%d\n", "WN = ", pt->EWNS ('W', 'N')->Index ());
  // 점의 왼쪽의 아래쪽점이 존재하는 경우, 점의 왼쪽의 아래쪽점의 index를 출력
  if (pt->EWNS ('W', 'S') != NULL) printf ("%-6s%d\n", "WS = ", pt->EWNS ('W', 'S')->Index ());
  // 점의 위쪽의 오른쪽점이 존재하는 경우, 점의 위쪽의 오른쪽점의 index를 출력
  if (pt->EWNS ('N', 'E') != NULL) printf ("%-6s%d\n", "NE = ", pt->EWNS ('N', 'E')->Index ());
  // 점의 위쪽의 왼쪽점이 존재하는 경우, 점의 위쪽의 왼쪽점의 index를 출력
  if (pt->EWNS ('N', 'W') != NULL) printf ("%-6s%d\n", "NW = ", pt->EWNS ('N', 'W')->Index ());
  // 점의 아래쪽의 오른쪽점이 존재하는 경우, 점의 아래쪽의 오른쪽점의 index를 출력
  if (pt->EWNS ('S', 'E') != NULL) printf ("%-6s%d\n", "SE = ", pt->EWNS ('S', 'E')->Index ());
  // 점의 아래쪽의 왼쪽점이 존재하는 경우, 점의 아래쪽의 왼쪽점의 index를 출력
  if (pt->EWNS ('S', 'W') != NULL) printf ("%-6s%d\n", "SW = ", pt->EWNS ('S', 'W')->Index ());
  // arrEnt의 원소를 출력
  for (size_t j = 0; j < 26; j++) {
    // arrEnt의 원소를 출력
    printf ("arrEnt[%2zu] = %f\n", j, this->arrEnt[j]);
  }
  printf ("rb = %23.16e\n", this->rb[node]);
  // node가 내부점의 개수보다 작은 경우
  return *this;
}

MatrixProcess & MatrixProcess::SettingArray (int node, AxialData *adat, Point *pt, xData *xdat, yData *ydat) {
  for (size_t i = 0; i < 26; i++) {
    // 오른쪽, 왼쪽, 위쪽, 아래쪽의 점들의 주소를 저장하기 위한 array
    this->uniqInt[i] = -1;
    // 오른쪽, 왼쪽, 위쪽, 아래쪽의 점에서의 값을 저장하기 위한 array
    this->uniqEnt[i] = 0.0;
    // arrInt에서 중복되는 점들을 제거한 index를 저장하기 위한 array
    this->arrInt[i] = NULL;
    // arrInt에서 중복되는 점들을 제거한 점에서의 값을 저장하기 위한 array
    this->arrEnt[i] = 0.0;
  }
  this->SettingAzimuth (pt, adat, node);

  // 행렬이 2n x 2n 이기 때문에, 14 - 26번째 원소는 1 - 13번째 원소의 반복
  for (size_t i = 13; i < 26; i++) {
    // arrInt의 원소가 존재하는 경우
    if (this->arrInt[i - 13] != NULL) {
      // arrInt가 저장한 주소를 그대로 저장
      this->arrInt[i] = this->arrInt[i - 13];
    }
  }

  // x-축선에서의 식과 y-축선에서의 식의 합을 계산
  // 현재점에서의 u의 계수
  this->arrEnt[0 ] = xdat->Cu + ydat->Cu;
  // 오른쪽점에서의 u의 계수
  this->arrEnt[1 ] = xdat->Eu;
  // 왼쪽점에서의 u의 계수
  this->arrEnt[2 ] = xdat->Wu;
  // 위쪽점에서의 u의 계수
  this->arrEnt[3 ] = ydat->Nu;
  // 아래쪽점에서의 u의 계수
  this->arrEnt[4 ] = ydat->Su;
  // 오른쪽의 위쪽점에서의 u의 계수
  this->arrEnt[5 ] = xdat->ENu;
  // 오른쪽의 아래쪽점에서의 u의 계수
  this->arrEnt[6 ] = xdat->ESu;
  // 왼쪽의 위쪽점에서의 u의 계수
  this->arrEnt[7 ] = xdat->WNu;
  // 왼쪽의 아래쪽점에서의 u의 계수
  this->arrEnt[8 ] = xdat->WSu;
  // 위쪽의 오른쪽점에서의 u의 계수
  this->arrEnt[9 ] = ydat->NEu;
  // 위쪽의 왼쪽점에서의 u의 계수
  this->arrEnt[10] = ydat->NWu;
  // 아래쪽의 오른쪽점에서의 u의 계수
  this->arrEnt[11] = ydat->SEu;
  // 아래쪽의 왼쪽점에서의 u의 계수
  this->arrEnt[12] = ydat->SWu;
  // 현재점에서의 phi의 계수
  this->arrEnt[13] = xdat->Cphi + ydat->Cphi;
  // 오른쪽점에서의 phi의 계수
  this->arrEnt[14] = xdat->Ephi;
  // 왼쪽점에서의 phi의 계수
  this->arrEnt[15] = xdat->Wphi;
  // 위쪽점에서의 phi의 계수
  this->arrEnt[16] = ydat->Nphi;
  // 아래쪽점에서의 phi의 계수
  this->arrEnt[17] = ydat->Sphi;
  // 오른쪽의 위쪽점에서의 phi의 계수
  this->arrEnt[18] = xdat->ENphi;
  // 오른쪽의 아래쪽점에서의 phi의 계수
  this->arrEnt[19] = xdat->ESphi;
  // 왼쪽의 위쪽점에서의 phi의 계수
  this->arrEnt[20] = xdat->WNphi;
  // 왼쪽의 아래쪽점에서의 phi의 계수
  this->arrEnt[21] = xdat->WSphi;
  // 위쪽의 오른쪽점에서의 phi의 계수
  this->arrEnt[22] = ydat->NEphi;
  // 위쪽의 왼쪽점에서의 phi의 계수
  this->arrEnt[23] = ydat->NWphi;
  // 아래쪽의 오른쪽점에서의 phi의 계수
  this->arrEnt[24] = ydat->SEphi;
  // 아래쪽의 왼쪽점에서의 phi의 계수
  this->arrEnt[25] = ydat->SWphi;
  // arrEnt의 원소가 0인 경우 arrInt의 원소의 값을 초기화
  for (size_t i = 0; i < 26; i++) {
    // arrEnt의 원소가 0인 경우
    if (IsEqualDouble (this->arrEnt[i], ZeroValue)) {
      // arrInt의 원소의 값을 초기화
      this->arrInt[i] = NULL;
    }
  }

  return *this;
}

MatrixProcess & MatrixProcess::PrintDebuggingData (AxialData *adat, Point *pt, xData *xdat, yData *ydat) {
  // int node = pt->Index ();
  int node = pt->Index () + adat->In_Pts_Num ();
  CalcRepresenCoef (pt, xdat, ydat, false);
  TransposeBoundaryData (pt, xdat, ydat, false);
  this->SettingArray (node, adat, pt, xdat, ydat);

  if (*pt == 4005) {
    // pt->PrintDebuggingData ("11100", "pt");
    for (size_t i = 0 ; i < 13; i++) if (this->arrInt[i]) printf ("\n"),
    printf ("arrEnt[%02zu] * u_ftn   = %23.16e, ", i, this->arrEnt[i] * u_ftn (this->arrInt[i])),
    printf ("arrEnt[%02zu] = %23.16e, "          , i, this->arrEnt[i]),
    printf ("u_ftn   = %23.16e"                , u_ftn (this->arrInt[i]));
    // this->arrInt[i]->PrintDebuggingData ("11100", "-----");

    for (size_t i = 13; i < 26; i++) if (this->arrInt[i]) printf ("\n"),
    printf ("arrEnt[%02zu] * phi_ftn = %23.16e, ", i, this->arrEnt[i] * u_ftn_Dirichlet (this->arrInt[i]->Phi ())),
    printf ("arrEnt[%02zu] = %23.16e, "          , i, this->arrEnt[i]),
    printf ("phi_ftn = %23.16e"                , u_ftn_Dirichlet (this->arrInt[i]->Phi ()));
    // this->arrInt[i]->PrintDebuggingData ("11100", "-----");
    // exit (112);
  }

  // printf ("xdat->Cu   = %23.16e\n", xdat->Cu);
  // printf ("xdat->Eu   = %23.16e\n", xdat->Eu);
  // printf ("xdat->Wu   = %23.16e\n", xdat->Wu);
  // printf ("xdat->Cphi = %23.16e\n", xdat->Cphi);
  // printf ("xdat->Ephi = %23.16e\n", xdat->Ephi);
  // printf ("xdat->wphi = %23.16e\n", xdat->Wphi);
  // printf ("xdat->F    = %23.16e\n", xdat->F);
  // printf ("\n");
  // printf ("ydat->Cu   = %23.16e\n", ydat->Cu);
  // printf ("ydat->Nu   = %23.16e\n", ydat->Nu);
  // printf ("ydat->Su   = %23.16e\n", ydat->Su);
  // printf ("ydat->Cphi = %23.16e\n", ydat->Cphi);
  // printf ("ydat->Nphi = %23.16e\n", ydat->Nphi);
  // printf ("ydat->Sphi = %23.16e\n", ydat->Sphi);
  // printf ("ydat->F    = %23.16e\n", ydat->F);
  // printf ("\n");
  // printf ("value      = %23.16e\n", u_ftn (pt));
  // printf ("\n");
  //
  //
  // pt->PrintDebuggingData ("11100", "pt");
  // for (size_t i = 0 ; i < 13; i++) if (this->arrInt[i]) printf ("\n"),
  // printf ("arrEnt[%02zu] * u_ftn   = %23.16e, ", i, this->arrEnt[i] * u_ftn (this->arrInt[i])),
  // printf ("arrEnt[%02zu] = %23.16e, "          , i, this->arrEnt[i]),
  // printf ("u_ftn   = %23.16e"                , u_ftn (this->arrInt[i])),
  // this->arrInt[i]->PrintDebuggingData ("11100", "-----");
  //
  // for (size_t i = 13; i < 26; i++) if (this->arrInt[i]) printf ("\n"),
  // printf ("arrEnt[%02zu] * phi_ftn = %23.16e, ", i, this->arrEnt[i] * u_ftn_Dirichlet (this->arrInt[i]->Phi ())),
  // printf ("arrEnt[%02zu] = %23.16e, "          , i, this->arrEnt[i]),
  // printf ("phi_ftn = %23.16e"                , u_ftn_Dirichlet (this->arrInt[i]->Phi ())),
  // this->arrInt[i]->PrintDebuggingData ("11100", "-----");

  double result = ZeroValue;

  for (size_t i = 0 ; i < 13; i++) if (this->arrInt[i]) result += this->arrEnt[i] * u_ftn_Dirichlet (this->arrInt[i]);

  for (size_t i = 13; i < 26; i++) if (this->arrInt[i]) result += this->arrEnt[i] * u_ftn_Dirichlet (this->arrInt[i]->Phi ());

  result -= xdat->F + ydat->F;

  // printf ("xdat->F = %23.16e\n", xdat->F);
  // printf ("ydat->F = %23.16e\n", ydat->F);

  // FILE *verify1;
  // verify1 = fopen ("verify1.txt", "a");
  // if (verify1!=NULL) {
  //   if (*pt == 'C')
  //   fprintf (verify1, "%c, "     , pt->Condition ()),
  //   fprintf (verify1, "%23.16e\n", result);
  //   fclose (verify1);
  // }
  //
  // FILE *verify2;
  // verify2 = fopen ("verify2.txt", "a");
  // if (verify2!=NULL) {
  //   if (*pt == 'D') if (fabs (result) < 0.01)
  //   fprintf (verify2, "%d, ", pt->Index ()), fprintf (verify2, "%c, ", pt->Condition ());
  //   if (*pt == 'D') if (fabs (result) < 0.01) for (auto &i : {'E', 'W', 'N', 'S'}) if (pt->EWNS (i, i)) if (pt->EWNS (i, i) != pt) if (*pt->EWNS (i, i) == 'C')
  //   fprintf (verify2, "%c"  , i);
  //   if (*pt == 'D') if (fabs (result) < 0.01) fprintf (verify2, ", ");
  //   if (*pt == 'D') if (fabs (result) < 0.01) fprintf (verify2, "%23.16e\n", result);
  //   fclose (verify2);
  // }
  //
  // FILE *verify3;
  // verify3 = fopen ("verify3.txt", "a");
  // if (verify3!=NULL) {
  //   if (*pt == 'D') if (fabs (result) >= 1.0E-2)
  //   fprintf (verify3, "%d, ", pt->Index ()), fprintf (verify3, "%c, ", pt->Condition ());
  //   if (*pt == 'D') if (fabs (result) >= 1.0E-2) for (auto &i : {'E', 'W', 'N', 'S'}) if (pt->EWNS (i, i)) if (pt->EWNS (i, i) != pt) if (*pt->EWNS (i, i) == 'C')
  //   fprintf (verify3, "%c", i);
  //   if (*pt == 'D') if (fabs (result) >= 1.0E-2) fprintf (verify3, ", ");
  //   if (*pt == 'D') if (fabs (result) >= 1.0E-2) fprintf (verify3, "%23.16e\n", result);
  //   fclose (verify3);
  // }

  // if (fabs (result) >= 1.0E-2) printf ("result = %23.16e\n", result);

  printf ("idx = %d/%d(%d), condition = %c, node = %d, result = %23.16e\n", pt->Index (), adat->Pts_Num (), adat->In_Pts_Num (), pt->Condition (), node, result);
  // if (fabs (result) > 1.0E-5) printf ("idx = %d/%d(%d), condition = %c, node = %d, result = %23.16e\n", pt->Index (), adat->Pts_Num (), adat->In_Pts_Num (), pt->Condition (), node, result);
  return *this;
}

MatrixProcess & MatrixProcess::PrintDebuggingData (AxialData *adat, Point *pt, xData *xdat, yData *ydat, double *returnValue) {
  int node = pt->Index ();
  // int node = pt->Index () + adat->In_Pts_Num ();
  CalcRepresenCoef (pt, xdat, ydat, true);
  TransposeBoundaryData (pt, xdat, ydat, true);
  this->SettingArray (node, adat, pt, xdat, ydat);

  // if (*pt == 8161) {
  //   pt->PrintDebuggingData ("11100", "pt");
  //   for (size_t i = 0 ; i < 13; i++) if (this->arrInt[i]) printf ("\n"),
  //   printf ("arrEnt[%02zu] * u_ftn   = %23.16e, ", i, this->arrEnt[i] * u_ftn_Dirichlet (this->arrInt[i])),
  //   printf ("arrEnt[%02zu] = %23.16e, "          , i, this->arrEnt[i]),
  //   printf ("u_ftn   = %23.16e"                , u_ftn_Dirichlet (this->arrInt[i]));
  //   // this->arrInt[i]->PrintDebuggingData ("11100", "-----");
  //
  //   for (size_t i = 13; i < 26; i++) if (this->arrInt[i]) printf ("\n"),
  //   printf ("arrEnt[%02zu] * phi_ftn = %23.16e, ", i, this->arrEnt[i] * u_ftn_Dirichlet (this->arrInt[i]->Phi ())),
  //   printf ("arrEnt[%02zu] = %23.16e, "          , i, this->arrEnt[i]),
  //   printf ("phi_ftn = %23.16e"                , u_ftn_Dirichlet (this->arrInt[i]->Phi ()));
  //   // this->arrInt[i]->PrintDebuggingData ("11100", "-----");
  //   printf ("\nrb = %23.16e\n", xdat->F + ydat->F);
  //   // exit (112);
  // }

  // printf ("xdat->Cu   = %23.16e\n", xdat->Cu);
  // printf ("xdat->Eu   = %23.16e\n", xdat->Eu);
  // printf ("xdat->Wu   = %23.16e\n", xdat->Wu);
  // printf ("xdat->Cphi = %23.16e\n", xdat->Cphi);
  // printf ("xdat->Ephi = %23.16e\n", xdat->Ephi);
  // printf ("xdat->wphi = %23.16e\n", xdat->Wphi);
  // printf ("xdat->F    = %23.16e\n", xdat->F);
  // printf ("\n");
  // printf ("ydat->Cu   = %23.16e\n", ydat->Cu);
  // printf ("ydat->Nu   = %23.16e\n", ydat->Nu);
  // printf ("ydat->Su   = %23.16e\n", ydat->Su);
  // printf ("ydat->Cphi = %23.16e\n", ydat->Cphi);
  // printf ("ydat->Nphi = %23.16e\n", ydat->Nphi);
  // printf ("ydat->Sphi = %23.16e\n", ydat->Sphi);
  // printf ("ydat->F    = %23.16e\n", ydat->F);
  // printf ("\n");
  // printf ("value      = %23.16e\n", u_ftn (pt));
  // printf ("\n");
  //
  // pt->PrintDebuggingData ("11111", xdat, ydat, true);
  //
  //
  // pt->PrintDebuggingData ("11100", "pt");
  // for (size_t i = 0 ; i < 13; i++) if (this->arrInt[i]) printf ("\n"),
  // printf ("arrEnt[%02zu] * u_ftn   = %23.16e, ", i, this->arrEnt[i] * u_ftn (this->arrInt[i])),
  // printf ("arrEnt[%02zu] = %23.16e, "          , i, this->arrEnt[i]),
  // printf ("u_ftn   = %23.16e"                , u_ftn (this->arrInt[i])),
  // this->arrInt[i]->PrintDebuggingData ("11100", "-----");
  //
  // for (size_t i = 13; i < 26; i++) if (this->arrInt[i]) printf ("\n"),
  // printf ("arrEnt[%02zu] * phi_ftn = %23.16e, ", i, this->arrEnt[i] * u_ftn_Dirichlet (this->arrInt[i]->Phi ())),
  // printf ("arrEnt[%02zu] = %23.16e, "          , i, this->arrEnt[i]),
  // printf ("phi_ftn = %23.16e"                , u_ftn_Dirichlet (this->arrInt[i]->Phi ())),
  // this->arrInt[i]->PrintDebuggingData ("11100", "-----");

  double result = ZeroValue;

  for (size_t i = 0 ; i < 13; i++) if (this->arrInt[i]) result += this->arrEnt[i] * u_ftn_Dirichlet (this->arrInt[i]);

  for (size_t i = 13; i < 26; i++) if (this->arrInt[i]) result += this->arrEnt[i] * u_ftn_Dirichlet (this->arrInt[i]->Phi ());

  result -= xdat->F + ydat->F;

  // printf ("result = %23.16e\n", result);
  *returnValue = result;

  // printf ("xdat->F = %23.16e\n", xdat->F);
  // printf ("ydat->F = %23.16e\n", ydat->F);

  // FILE *verify1;
  // verify1 = fopen ("verify1.txt", "a");
  // if (verify1!=NULL) {
  //   if (*pt == 'C')
  //   fprintf (verify1, "%c, "     , pt->Condition ()),
  //   fprintf (verify1, "%23.16e\n", result);
  //   fclose (verify1);
  // }
  //
  // FILE *verify2;
  // verify2 = fopen ("verify2.txt", "a");
  // if (verify2!=NULL) {
  //   if (*pt == 'D') if (fabs (result) < 0.01)
  //   fprintf (verify2, "%d, ", pt->Index ()), fprintf (verify2, "%c, ", pt->Condition ());
  //   if (*pt == 'D') if (fabs (result) < 0.01) for (auto &i : {'E', 'W', 'N', 'S'}) if (pt->EWNS (i, i)) if (pt->EWNS (i, i) != pt) if (*pt->EWNS (i, i) == 'C')
  //   fprintf (verify2, "%c"  , i);
  //   if (*pt == 'D') if (fabs (result) < 0.01) fprintf (verify2, ", ");
  //   if (*pt == 'D') if (fabs (result) < 0.01) fprintf (verify2, "%23.16e\n", result);
  //   fclose (verify2);
  // }
  //
  // FILE *verify3;
  // verify3 = fopen ("verify3.txt", "a");
  // if (verify3!=NULL) {
  //   if (*pt == 'D') if (fabs (result) >= 1.0E-2)
  //   fprintf (verify3, "%d, ", pt->Index ()), fprintf (verify3, "%c, ", pt->Condition ());
  //   if (*pt == 'D') if (fabs (result) >= 1.0E-2) for (auto &i : {'E', 'W', 'N', 'S'}) if (pt->EWNS (i, i)) if (pt->EWNS (i, i) != pt) if (*pt->EWNS (i, i) == 'C')
  //   fprintf (verify3, "%c", i);
  //   if (*pt == 'D') if (fabs (result) >= 1.0E-2) fprintf (verify3, ", ");
  //   if (*pt == 'D') if (fabs (result) >= 1.0E-2) fprintf (verify3, "%23.16e\n", result);
  //   fclose (verify3);
  // }

  // if (fabs (result) >= 1.0E-2) printf ("result = %23.16e\n", result);

  // printf ("idx = %d/%d(%d), condition = %c, node = %d, result = %23.16e\n", pt->Index (), adat->Pts_Num (), adat->In_Pts_Num (), pt->Condition (), node, result);
  // if (fabs (result) > 1.0E-5) printf ("idx = %d/%d(%d), condition = %c, node = %d, result = %23.16e\n", pt->Index (), adat->Pts_Num (), adat->In_Pts_Num (), pt->Condition (), node, result);
  return *this;
}

#endif
