#ifndef GREENSFUNCTIONS_H
#define GREENSFUNCTIONS_H

#include "util.hpp"

/* 그린함수 */
double greens_function (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
  /*
  nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  2차원의 무한대에서의 점근적함수는 u (x, y) = log (x^2 + y^2)
  3차원축대칭의 무한대에서의 점근적함수는 u (r, z) = 1 / (r^2 + z^2)

  node: x축 성분인지 y축 성분인지를 나타내는 변수
  1: x축 성분
  2: y축 성분

  bdc: boundary condition
  0: Dirichlet - Dirichlet
  1: Neumann - Dirichlet
  2: Dirichlet - Neumann
  3: Infinity - Dirichlet
  4: Dirichlet - Infinity
  5: Singularity - Dirichlet
  6: Dirichlet - Singularity

  mp1: 왼쪽(아래쪽) 부분의 conductivity
  mp2: 오른쪽(위쪽) 부분의 conductivity
  */

  // 차원
  int nD = 2;

  // node 값의 에러 출력
  if (node > 2) {printf ("%s%d\n", "greens_function node error, node = ", node); exit (1);}

  // 2차원
  if (nD == 2) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      if (IsEqualDouble (t, tm) || t < tb) return ((tb-tp)*(tm-t))/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
      else        return -((t-tp)*(tb-tm))/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      if (IsEqualDouble (t, tm) || t < tb) return -(tb-tp)/mp2;
      else        return -(t-tp)/mp2;
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      if (IsEqualDouble (t, tm) || t < tb) return (t-tm)/mp1;
      else        return (tb-tm)/mp1;
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      if (IsEqualDouble (t, tm) || t < tb) return -(1.0/(t*t)*(tp*tp)*(tb-tp))/mp2;
      else        return (1.0/(t*t)*((t*t)*tb-tb*(tp*tp)-t*t*t+tp*tp*tp))/mp2;
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      if (IsEqualDouble (t, tm) || t < tb) return -(1.0/(t*t)*((t*t)*tb-tb*(tm*tm)-t*t*t+tm*tm*tm))/mp1;
      else        return (1.0/(t*t)*(tm*tm)*(tb-tm))/mp1;
    }
    // Singularity - Dirichlet
    else if (bdc == 5) {
      if (t < tb) return -(pow(t-tm,3.0)*(tb-tp)*1.0/pow(tm-tp,3.0))/mp2;
      else        return -((t-tp)*1.0/pow(tm-tp,3.0)*((t*t)*tb-t*(tp*tp)-(t*t)*tp+tb*(tm*tm)*3.0+tb*(tp*tp)-tm*tm*tm-t*tb*tm*3.0+t*tb*tp+t*tm*tp*3.0-tb*tm*tp*3.0))/mp2;
    }
    // Dirichlet - Singularity
    else if (bdc == 6) {
      if (t < tb) return -((t-tm)*1.0/pow(tm-tp,3.0)*((t*t)*tb-t*(tm*tm)-(t*t)*tm+tb*(tm*tm)+tb*(tp*tp)*3.0-tp*tp*tp+t*tb*tm-t*tb*tp*3.0+t*tm*tp*3.0-tb*tm*tp*3.0))/mp1;
      else        return -(pow(t-tp,3.0)*(tb-tm)*1.0/pow(tm-tp,3.0))/mp1;
    }
    // bdc값의 에러 출력
    else {
      printf ("%s%d\n", "greens_function bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 3차원 축대칭
  else if (nD == 3) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return ((t-tm)*(tb-tp))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        else        return ((t-tp)*(tb-tm))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return -((log(fabs(t))-log(fabs(tm)))*(log(fabs(tb))-log(fabs(tp))))/(mp1*mp2*((log(fabs(tb))-log(fabs(tm)))/mp1-(log(fabs(tb))-log(fabs(tp)))/mp2));
        else        return -((log(fabs(t))-log(fabs(tp)))*(log(fabs(tb))-log(fabs(tm))))/(mp1*mp2*((log(fabs(tb))-log(fabs(tm)))/mp1-(log(fabs(tb))-log(fabs(tp)))/mp2));
      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      if (node == 1) {
        // x축 성분
        if (t < tb) return -(tb-tp)/(mp2*yb);
        else        return -(t-tp)/(mp2*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return -(log(fabs(tb))-log(fabs(tp)))/mp2;
        else        return -(log(fabs(t))-log(fabs(tp)))/mp2;
      }
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return (t-tm)/(mp1*yb);
        else        return (tb-tm)/(mp1*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return (log(fabs(t))-log(fabs(tm)))/mp1;
        else        return (log(fabs(tb))-log(fabs(tm)))/mp1;
      }
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return -(1.0/(t*t)*(tp*tp)*(tb-tp))/(mp2*yb);
        else        return (1.0/(t*t)*((t*t)*tb-tb*(tp*tp)-t*t*t+tp*tp*tp))/(mp2*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return 0.0;
        else        return -(log(fabs(t))-log(fabs(tp)))/mp2;
      }
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return -(1.0/(t*t)*((t*t)*tb-tb*(tm*tm)-t*t*t+tm*tm*tm))/(mp1*yb);
        else        return (1.0/(t*t)*(tm*tm)*(tb-tm))/(mp1*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return (log(fabs(t))-log(fabs(tb)))/mp1+(1.0/(t*t*t)*((tm*tm*tm)*log(fabs(tb))-(tm*tm*tm)*log(fabs(tm))))/mp1;
        else        return (1.0/(t*t*t)*(tm*tm*tm)*(log(fabs(tb))-log(fabs(tm))))/mp1;
      }
    }
    // bdc값의 에러 출력
    else {
      printf ("%s%d\n", "greens_function bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 차원의 에러 출력
  else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

/* 그린함수의 첫번째 변수로의 미분 */
double greens_function_t (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
  /*
  nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  2차원의 무한대에서의 점근적함수는 u (x, y) = log (x^2 + y^2)
  3차원축대칭의 무한대에서의 점근적함수는 u (r, z) = 1 / (r^2 + z^2)

  node: x축 성분인지 y축 성분인지를 나타내는 변수
  1: x축 성분
  2: y축 성분

  bdc: boundary condition
  0: Dirichlet - Dirichlet
  1: Neumann - Dirichlet
  2: Dirichlet - Neumann
  3: Infinity - Dirichlet
  4: Dirichlet - Infinity
  5: Singularity - Dirichlet
  6: Dirichlet - Singularity

  mp1: 왼쪽(아래쪽) 부분의 conductivity
  mp2: 오른쪽(위쪽) 부분의 conductivity
  */

  // 차원
  int nD = 2;

  // node 값의 에러 출력
  if (node > 2) {printf ("%s%d\n", "greens_function_t node error, node = ", node); exit (1);}

  // 2차원
  if (nD == 2) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      if (IsEqualDouble (t, tm) || t < tb) return -(tb-tp)/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
      else        return -(tb-tm)/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      if (IsEqualDouble (t, tm) || t < tb) return 0.0;
      else        return -1.0/mp2;
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      if (IsEqualDouble (t, tm) || t < tb) return 1.0/mp1;
      else        return 0.0;
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      if (IsEqualDouble (t, tm) || t < tb) return (1.0/(t*t*t)*(tp*tp)*(tb-tp)*2.0)/mp2;
      else        return -(1.0/(t*t*t)*(tb*(tp*tp)*-2.0+t*t*t+(tp*tp*tp)*2.0))/mp2;
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      if (IsEqualDouble (t, tm) || t < tb) return (1.0/(t*t*t)*(tb*(tm*tm)*-2.0+t*t*t+(tm*tm*tm)*2.0))/mp1;
      else        return (1.0/(t*t*t)*(tm*tm)*(tb-tm)*-2.0)/mp1;
    }
    // Singularity - Dirichlet
    else if (bdc == 5) {
      if (t < tb) return (pow(t-tm,2.0)*(tb-tp)*1.0/pow(tm-tp,3.0)*-3.0)/mp2;
      else        return -(1.0/pow(tm-tp,3.0)*((t*t)*tb*3.0-(t*t)*tp*3.0+tb*(tm*tm)*3.0-tm*(tp*tp)*3.0-tm*tm*tm+tp*tp*tp-t*tb*tm*6.0+t*tm*tp*6.0))/mp2;
    }
    // Dirichlet - Singularity
    else if (bdc == 6) {
      if (t < tb) return -(1.0/pow(tm-tp,3.0)*((t*t)*tb*3.0-(t*t)*tm*3.0+tb*(tp*tp)*3.0-(tm*tm)*tp*3.0+tm*tm*tm-tp*tp*tp-t*tb*tp*6.0+t*tm*tp*6.0))/mp1;
      else        return (pow(t-tp,2.0)*(tb-tm)*1.0/pow(tm-tp,3.0)*-3.0)/mp1;
    }
    // bdc값의 에러 출력
    else {
      printf ("%s%d\n", "greens_function_t bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 3차원 축대칭
  else if (nD == 3) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return (tb-tp)/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        else        return (tb-tm)/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return -(log(fabs(tb))-log(fabs(tp)))/(t*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
        else        return -(log(fabs(tb))-log(fabs(tm)))/(t*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return 0.0;
        else        return -1.0/(mp2*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return 0.0;
        else        return -1.0/(mp2*t);
      }
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return 1.0/(mp1*yb);
        else        return 0.0;
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return 1.0/(mp1*t);
        else        return 0.0;
      }
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return (1.0/(t*t*t)*(tp*tp)*(tb-tp)*2.0)/(mp2*yb);
        else        return -(1.0/(t*t*t)*(tb*(tp*tp)*-2.0+t*t*t+(tp*tp*tp)*2.0))/(mp2*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return 0.0;
        else        return -1.0/(mp2*t);
      }
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return (1.0/(t*t*t)*(tb*(tm*tm)*-2.0+t*t*t+(tm*tm*tm)*2.0))/(mp1*yb);
        else        return (1.0/(t*t*t)*(tm*tm)*(tb-tm)*-2.0)/(mp1*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return (1.0/(t*t*t*t)*((tm*tm*tm)*log(fabs(tb))*-3.0+(tm*tm*tm)*log(fabs(tm))*3.0+t*t*t))/mp1;
        else        return (1.0/(t*t*t*t)*(tm*tm*tm)*(log(fabs(tb))-log(fabs(tm)))*-3.0)/mp1;
      }
    }
    // bdc값의 에러 출력
    else {
      printf ("%s%d\n", "greens_function_t bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 차원의 에러 출력
  else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

/* 그린함수의 두번째 변수로의 미분 */
double greens_function_tau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
  /*
  nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  2차원의 무한대에서의 점근적함수는 u (x, y) = log (x^2 + y^2)
  3차원축대칭의 무한대에서의 점근적함수는 u (r, z) = 1 / (r^2 + z^2)

  node: x축 성분인지 y축 성분인지를 나타내는 변수
  1: x축 성분
  2: y축 성분

  bdc: boundary condition
  0: Dirichlet - Dirichlet
  1: Neumann - Dirichlet
  2: Dirichlet - Neumann
  3: Infinity - Dirichlet
  4: Dirichlet - Infinity
  5: Singularity - Dirichlet
  6: Dirichlet - Singularity

  mp1: 왼쪽(아래쪽) 부분의 conductivity
  mp2: 오른쪽(위쪽) 부분의 conductivity
  */

  // 차원
  int nD = 2;

  // node 값의 에러 출력
  if (node > 2) {printf ("%s%d\n", "greens_function_tau node error, node = ", node); exit (1);}

  // 2차원
  if (nD == 2) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      if (IsEqualDouble (t, tm) || t < tb) return -(t-tm)/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
      else        return -(t-tp)/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      if (IsEqualDouble (t, tm) || t < tb) return -1.0/mp2;
      else        return 0.0;
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      if (IsEqualDouble (t, tm) || t < tb) return 0.0;
      else        return 1.0/mp1;
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      if (IsEqualDouble (t, tm) || t < tb) return -(1.0/(t*t)*(tp*tp))/mp2;
      else        return (1.0/(t*t)*(t*t-tp*tp))/mp2;
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      if (IsEqualDouble (t, tm) || t < tb) return -(1.0/(t*t)*(t*t-tm*tm))/mp1;
      else        return (1.0/(t*t)*(tm*tm))/mp1;
    }
    // bdc값의 에러 출력
    else {
      printf ("%s%d\n", "greens_function_tau bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 3차원 축대칭
  else if (nD == 3) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return (t-tm)/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        else        return (t-tp)/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return -(log(fabs(t))-log(fabs(tm)))/(tb*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
        else        return -(log(fabs(t))-log(fabs(tp)))/(tb*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return -1.0/(mp2*yb);
        else        return 0.0;
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return -1.0/(mp2*tb);
        else        return 0.0;
      }
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return 0.0;
        else        return 1.0/(mp1*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return 0.0;
        else        return 1.0/(mp1*tb);
      }
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return -(1.0/(t*t)*(tp*tp))/(mp2*yb);
        else        return (1.0/(t*t)*(t*t-tp*tp))/(mp2*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return 0.0;
        else        return 0.0;
      }
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return -(1.0/(t*t)*(t*t-tm*tm))/(mp1*yb);
        else        return (1.0/(t*t)*(tm*tm))/(mp1*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return -(1.0/(t*t*t)*(t*t*t-tm*tm*tm))/(mp1*tb);
        else        return (1.0/(t*t*t)*(tm*tm*tm))/(mp1*tb);
      }
    }
    // bdc값의 에러 출력
    else {
      printf ("%s%d\n", "greens_function_tau bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 차원의 에러 출력
  else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

/* 그린함수의 첫번째, 두뻔째 변수로의 미분 */
double greens_function_ttau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
  /*
  nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  2차원의 무한대에서의 점근적함수는 u (x, y) = log (x^2 + y^2)
  3차원축대칭의 무한대에서의 점근적함수는 u (r, z) = 1 / (r^2 + z^2)

  node: x축 성분인지 y축 성분인지를 나타내는 변수
  1: x축 성분
  2: y축 성분

  bdc: boundary condition
  0: Dirichlet - Dirichlet
  1: Neumann - Dirichlet
  2: Dirichlet - Neumann
  3: Infinity - Dirichlet
  4: Dirichlet - Infinity
  5: Singularity - Dirichlet
  6: Dirichlet - Singularity

  mp1: 왼쪽(아래쪽) 부분의 conductivity
  mp2: 오른쪽(위쪽) 부분의 conductivity
  */

  int nD = 2;

  // node 값의 에러 출력
  if (node > 2) {printf ("%s%d\n", "greens_function_ttau node error, node = ", node); exit (1);}

  // 2차원
  if (nD == 2) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      if (IsEqualDouble (t, tm) || t < tb) return 1.0/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      else        return 1.0/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      if (IsEqualDouble (t, tm) || t < tb) return 0.0;
      else        return 0.0;
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      if (IsEqualDouble (t, tm) || t < tb) return 0.0;
      else        return 0.0;
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      if (IsEqualDouble (t, tm) || t < tb) return (1.0/(t*t*t)*(tp*tp)*2.0)/mp2;
      else        return (1.0/(t*t*t)*(tp*tp)*2.0)/mp2;
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      if (IsEqualDouble (t, tm) || t < tb) return (1.0/(t*t*t)*(tm*tm)*-2.0)/mp1;
      else        return (1.0/(t*t*t)*(tm*tm)*-2.0)/mp1;
    }
    // bdc값의 에러 출력
    else {
      printf ("%s%d\n", "greens_function_ttau bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 3차원 축대칭
  else if (nD == 3) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return 1.0/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        else        return 1.0/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return -1.0/(t*tb*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
        else        return -1.0/(t*tb*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      if (t < tb) return 0.0;
      else        return 0.0;
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      if (t < tb) return 0.0;
      else        return 0.0;
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return (1.0/(t*t*t)*(tp*tp)*2.0)/(mp2*yb);
        else        return (1.0/(t*t*t)*(tp*tp)*2.0)/(mp2*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return 0.0;
        else        return 0.0;
      }
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      // x축 성분
      if (node == 1) {
        if (t < tb) return (1.0/(t*t*t)*(tm*tm)*-2.0)/(mp1*yb);
        else        return (1.0/(t*t*t)*(tm*tm)*-2.0)/(mp1*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (t < tb) return (1.0/(t*t*t*t)*(tm*tm*tm)*-3.0)/(mp1*tb);
        else        return (1.0/(t*t*t*t)*(tm*tm*tm)*-3.0)/(mp1*tb);
      }
    }
    // bdc값의 에러 출력
    else {
      printf ("%s%d\n", "greens_function_ttau bdc error, bdc = ", bdc);
      exit (1);
    }
    // 차원의 에러 출력
  } else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

/* 그린함수와 basis function의 적분 */
double greens_integral (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
  /*
  nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  2차원의 무한대에서의 점근적함수는 u (x, y) = log (x^2 + y^2)
  3차원축대칭의 무한대에서의 점근적함수는 u (r, z) = 1 / (r^2 + z^2)

  index: 그린함수와 basis function의 적분을 가리키는 변수
  1: 왼쪽(아래쪽) 점에서 1인 piecewise linear basis function
  2: 가운데 점에서 1인 tau보다 큰 부분에서는 0인 basis function
  3: 가운데 점에서 1인 tau보다 작은 부분에서는 0인 basis function
  4: 오른쪽(위쪽) 점에서 1인 piecewise linear basis function

  node: x축 성분인지 y축 성분인지를 나타내는 변수
  1: x축 성분
  2: y축 성분

  bdc: boundary condition
  0: Dirichlet - Dirichlet
  1: Neumann - Dirichlet
  2: Dirichlet - Neumann
  3: Infinity - Dirichlet
  4: Dirichlet - Infinity
  5: Singularity - Dirichlet
  6: Dirichlet - Singularity

  mp1: 왼쪽(아래쪽) 부분의 conductivity
  mp2: 오른쪽(위쪽) 부분의 conductivity
  */

  // 차원
  int nD = 2;

  // node 값의 에러 출력
  if (node > 2) {printf ("%s%d\n", "greens_integral node error, node = ", node); exit (1);}

  // bdc 값의 에러 출력
  if (bdc > 6) {printf ("%s%d\n", "greens_integral bdc error, bdc = ", bdc); exit (1);}

  // 2차원
  if (nD ==2) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // 양 쪽의 conductivity가 같은 경우
      if (IsEqualDouble (mp1, mp2)) {
        if (index == 1) return (pow(tb-tm,2.0)*(tb-tp)*(1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 2) return (pow(tb-tm,2.0)*(tb-tp)*(1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 3) return ((tb-tm)*pow(tb-tp,2.0)*(-1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 4) return ((tb-tm)*pow(tb-tp,2.0)*(-1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      }
      // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
      else {
        if (index == 1) return (pow(tb-tm,2.0)*(tb-tp)*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return ((tb-tm)*pow(tb-tp,2.0)*(-1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);

      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      if (index == 1) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/mp2;
      if (index == 2) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/mp2;
      if (index == 3) return (pow(tb-tp,2.0)*(1.0/3.0))/mp2;
      if (index == 4) return (pow(tb-tp,2.0)*(1.0/6.0))/mp2;
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      if (index == 1) return (pow(tb-tm,2.0)*(1.0/6.0))/mp1;
      if (index == 2) return (pow(tb-tm,2.0)*(1.0/3.0))/mp1;
      if (index == 3) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/mp1;
      if (index == 4) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/mp1;
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return (tb*(tp*tp)*(-7.0/6.0)+(tb*tb)*tp*(1.0/3.0)-(tb*tb*tb)*(1.0/6.0)+tp*tp*tp+tb*(tp*tp)*log(fabs(tb))-tb*(tp*tp)*log(fabs(tp)))/(mp2*tb);
      if (index == 4) return (-(tp*tp)*log(fabs(tb))+(tp*tp)*log(fabs(tp))+tb*tp*(5.0/3.0)-(tb*tb)*(1.0/3.0)-(tp*tp)*(4.0/3.0))/mp2;
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      if (index == 1) return (-(tm*tm)*log(fabs(tb))+(tm*tm)*log(fabs(tm))+tb*tm*(5.0/3.0)-(tb*tb)*(1.0/3.0)-(tm*tm)*(4.0/3.0))/mp1;
      if (index == 2) return (tb*(tm*tm)*(-7.0/6.0)+(tb*tb)*tm*(1.0/3.0)-(tb*tb*tb)*(1.0/6.0)+tm*tm*tm+tb*(tm*tm)*log(fabs(tb))-tb*(tm*tm)*log(fabs(tm)))/(mp1*tb);
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    }
    // bdc 값의 에러 출력
    else {
      printf ("%s%d\n", "greens_integral bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 3차원 축대칭
  else if (nD == 3) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // x축 성분
      if (node == 1) {
        // 양 쪽의 conductivity가 같은 경우
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return (pow(tb-tm,2.0)*(tb-tp)*(1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return (pow(tb-tm,2.0)*(tb-tp)*(1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 3) return ((tb-tm)*pow(tb-tp,2.0)*(-1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 4) return ((tb-tm)*pow(tb-tp,2.0)*(-1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        }
        // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
        else {
          if (index == 1) return (pow(tb-tm,2.0)*(tb-tp)*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return ((tb-tm)*pow(tb-tp,2.0)*(-1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        }
      }
      // y축 성분
      else if (node == 2) {
        // 양 쪽의 conductivity가 같은 경우
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return ((log(fabs(tb))-log(fabs(tp)))*((tb*tb)*log(fabs(tb))*-2.0+(tb*tb)*log(fabs(tm))*2.0-tb*tm*4.0+(tb*tb)*3.0+tm*tm)*(1.0/4.0))/((tb-tm)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 2) return ((log(fabs(tb))-log(fabs(tp)))*((tb*tb)*log(fabs(tb))*-2.0+(tb*tb)*log(fabs(tm))*2.0-tb*tm*4.0+tb*tb+(tm*tm)*3.0+tb*tm*log(fabs(tb))*4.0-tb*tm*log(fabs(tm))*4.0)*(1.0/4.0))/((tb-tm)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 3) return ((log(fabs(tb))-log(fabs(tm)))*((tb*tb)*log(fabs(tb))*-2.0+(tb*tb)*log(fabs(tp))*2.0-tb*tp*4.0+tb*tb+(tp*tp)*3.0+tb*tp*log(fabs(tb))*4.0-tb*tp*log(fabs(tp))*4.0)*(-1.0/4.0))/((tb-tp)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 4) return ((log(fabs(tb))-log(fabs(tm)))*((tb*tb)*log(fabs(tb))*-2.0+(tb*tb)*log(fabs(tp))*2.0-tb*tp*4.0+(tb*tb)*3.0+tp*tp)*(-1.0/4.0))/((tb-tp)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
        }
        // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
        else {
          if (index == 1) return ((log(fabs(tb))-log(fabs(tp)))*(tb-tm-tb*log(fabs(tb))+tb*log(fabs(tm))))/(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp))));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return -((log(fabs(tb))-log(fabs(tm)))*(tb-tp+log(pow(tp,-tp))+log(pow(tp,tb))+log(pow(tp,tp))-tb*log(fabs(tb))))/(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp))));
        }
      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/(mp2*yb);
        if (index == 2) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/(mp2*yb);
        if (index == 3) return (pow(tb-tp,2.0)*(1.0/3.0))/(mp2*yb);
        if (index == 4) return (pow(tb-tp,2.0)*(1.0/6.0))/(mp2*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return ((log(fabs(tb))-log(fabs(tp)))*(tb-tm)*(-1.0/2.0))/mp2;
        if (index == 2) return ((log(fabs(tb))-log(fabs(tp)))*(tb-tm)*(-1.0/2.0))/mp2;
        if (index == 3) return ((tb*tb)*log(fabs(tb))*(1.0/2.0)-(tb*tb)*log(fabs(tp))*(1.0/2.0)+tb*tp-(tb*tb)*(1.0/4.0)-(tp*tp)*(3.0/4.0)-tb*tp*log(fabs(tb))+tb*tp*log(fabs(tp)))/(mp2*(tb-tp));
        if (index == 4) return ((tb*tb)*log(fabs(tb))*(1.0/2.0)-(tb*tb)*log(fabs(tp))*(1.0/2.0)+tb*tp-(tb*tb)*(3.0/4.0)-(tp*tp)*(1.0/4.0))/(mp2*(tb-tp));
      }
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return (pow(tb-tm,2.0)*(1.0/6.0))/(mp1*yb);
        if (index == 2) return (pow(tb-tm,2.0)*(1.0/3.0))/(mp1*yb);
        if (index == 3) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/(mp1*yb);
        if (index == 4) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/(mp1*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return ((tb*tb)*log(fabs(tb))*(1.0/2.0)-(tb*tb)*log(fabs(tm))*(1.0/2.0)+tb*tm-(tb*tb)*(3.0/4.0)-(tm*tm)*(1.0/4.0))/(mp1*(tb-tm));
        if (index == 2) return ((tb*tb)*log(fabs(tb))*(1.0/2.0)-(tb*tb)*log(fabs(tm))*(1.0/2.0)+tb*tm-(tb*tb)*(1.0/4.0)-(tm*tm)*(3.0/4.0)-tb*tm*log(fabs(tb))+tb*tm*log(fabs(tm)))/(mp1*(tb-tm));
        if (index == 3) return ((log(fabs(tb))-log(fabs(tm)))*(tb-tp)*(-1.0/2.0))/mp1;
        if (index == 4) return ((log(fabs(tb))-log(fabs(tm)))*(tb-tp)*(-1.0/2.0))/mp1;
      }
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return (tb*(tp*tp)*(-7.0/6.0)+(tb*tb)*tp*(1.0/3.0)-(tb*tb*tb)*(1.0/6.0)+tp*tp*tp+tb*(tp*tp)*log(fabs(tb))-tb*(tp*tp)*log(fabs(tp)))/(mp2*tb*yb);
        if (index == 4) return (-(tp*tp)*log(fabs(tb))+(tp*tp)*log(fabs(tp))+tb*tp*(5.0/3.0)-(tb*tb)*(1.0/3.0)-(tp*tp)*(4.0/3.0))/(mp2*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return ((tb*tb)*log(fabs(tb))*(1.0/2.0)-(tb*tb)*log(fabs(tp))*(1.0/2.0)+tb*tp-(tb*tb)*(1.0/4.0)-(tp*tp)*(3.0/4.0)-tb*tp*log(fabs(tb))+tb*tp*log(fabs(tp)))/(mp2*(tb-tp));
        if (index == 4) return ((tb*tb)*log(fabs(tb))*(1.0/2.0)-(tb*tb)*log(fabs(tp))*(1.0/2.0)+tb*tp-(tb*tb)*(3.0/4.0)-(tp*tp)*(1.0/4.0))/(mp2*(tb-tp));
      }
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return (-(tm*tm)*log(fabs(tb))+(tm*tm)*log(fabs(tm))+tb*tm*(5.0/3.0)-(tb*tb)*(1.0/3.0)-(tm*tm)*(4.0/3.0))/(mp1*yb);
        if (index == 2) return (tb*(tm*tm)*(-7.0/6.0)+(tb*tb)*tm*(1.0/3.0)-(tb*tb*tb)*(1.0/6.0)+tm*tm*tm+tb*(tm*tm)*log(fabs(tb))-tb*(tm*tm)*log(fabs(tm)))/(mp1*tb*yb);
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return (tb*(tm*tm)*(-1.0/4.0)+(tb*tb)*tm-(tb*tb*tb)*(3.0/4.0)+tm*(log(fabs(tb))-log(fabs(tm)))*(tb*tm*-3.0+(tb*tb)*3.0+tm*tm)*(1.0/2.0))/(mp1*tb*(tb-tm));
        if (index == 2) return (tb*(-1.0/4.0))/mp1+(tm*(3.0/4.0))/mp1-(1.0/(tb*tb)*((tm*tm*tm*tm)*log(fabs(tb))*(-1.0/2.0)+(tm*tm*tm*tm)*log(fabs(tm))*(1.0/2.0)+tb*((tm*tm*tm)*log(fabs(tb))-(tm*tm*tm)*log(fabs(tm)))))/(mp1*(tb-tm));
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    }
    // bdc 값의 에러 출력
    else {
      printf ("%s%d\n", "greens_integral bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 차원의 에러 출력
  else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

/* 그린함수의 첫번째 변수로의 미분과 basis function의 적분 */
double greens_integral_t (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
  /*
  nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  2차원의 무한대에서의 점근적함수는 u (x, y) = log (x^2 + y^2)
  3차원축대칭의 무한대에서의 점근적함수는 u (r, z) = 1 / (r^2 + z^2)

  index: 그린함수와 basis function의 적분을 가리키는 변수
  1: 왼쪽(아래쪽) 점에서 1인 piecewise linear basis function
  2: 가운데 점에서 1인 tau보다 큰 부분에서는 0인 basis function
  3: 가운데 점에서 1인 tau보다 작은 부분에서는 0인 basis function
  4: 오른쪽(위쪽) 점에서 1인 piecewise linear basis function

  node: x축 성분인지 y축 성분인지를 나타내는 변수
  1: x축 성분
  2: y축 성분

  bdc: boundary condition
  0: Dirichlet - Dirichlet
  1: Neumann - Dirichlet
  2: Dirichlet - Neumann
  3: Infinity - Dirichlet
  4: Dirichlet - Infinity
  5: Singularity - Dirichlet
  6: Dirichlet - Singularity

  mp1: 왼쪽(아래쪽) 부분의 conductivity
  mp2: 오른쪽(위쪽) 부분의 conductivity
  */

  // 차원
  int nD = 2;

  // node 값의 에러 출력
  if (node > 2) {printf ("%s%d\n", "greens_integral node error, node = ", node); exit (1);}

  // bdc 값의 에러 출력
  if (bdc > 6) {printf ("%s%d\n", "greens_integral bdc error, bdc = ", bdc); exit (1);}

  // 2차원
  if (nD ==2) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // 양 쪽의 conductivity가 같은 경우
      if (IsEqualDouble (mp1, mp2)) {
        if (index == 1) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 2) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 3) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 4) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      }
      // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
      else {
        if (index == 1) return -((tb-tm)*(tb-tp))/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return ((tb-tm)*(tb-tp))/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return (tb*(1.0/2.0)-tp*(1.0/2.0))/mp2;
      if (index == 4) return (tb*(1.0/2.0)-tp*(1.0/2.0))/mp2;
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      if (index == 1) return (tb*(1.0/2.0)-tm*(1.0/2.0))/mp1;
      if (index == 2) return (tb*(1.0/2.0)-tm*(1.0/2.0))/mp1;
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return (1.0/(tb*tb)*(tb*(tp*tp)*4.0-(tb*tb)*tp*3.0+tb*tb*tb-(tp*tp*tp)*2.0)*(1.0/2.0))/mp2;
      if (index == 4) return (tb*tp*(3.0/2.0)-(tb*tb)*(1.0/2.0)-tp*tp)/(mp2*tb);
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      if (index == 1) return (tb*tm*(3.0/2.0)-(tb*tb)*(1.0/2.0)-tm*tm)/(mp1*tb);
      if (index == 2) return (1.0/(tb*tb)*(tb*(tm*tm)*4.0-(tb*tb)*tm*3.0+tb*tb*tb-(tm*tm*tm)*2.0)*(1.0/2.0))/mp1;
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    }
    // bdc 값의 에러 출력
    else {
      printf ("%s%d\n", "greens_integral bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 3차원 축대칭
  else if (nD == 3) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // x축 성분
      if (node == 1) {
        // 양 쪽의 conductivity가 같은 경우
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 3) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 4) return ((tb-tm)*(tb-tp)*(-1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        }
        // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
        else {
          if (index == 1) return ((tb-tm)*(tb-tp))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return -((tb-tm)*(tb-tp))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        }
      }
      // y축 성분
      else if (node == 2) {
        // 양 쪽의 conductivity가 같은 경우
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return ((log(fabs(tb))-log(fabs(tp)))*(tb-tm-tb*log(fabs(tb))+tb*log(fabs(tm))))/((tb-tm)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 2) return -((log(fabs(tb))-log(fabs(tp)))*(tb-tm-tm*log(fabs(tb))+tm*log(fabs(tm))))/((tb-tm)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 3) return ((log(fabs(tb))-log(fabs(tm)))*(tb-tp-tp*log(fabs(tb))+tp*log(fabs(tp))))/((tb-tp)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 4) return -((log(fabs(tb))-log(fabs(tm)))*(tb-tp-tb*log(fabs(tb))+tb*log(fabs(tp))))/((tb-tp)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
        }
        // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
        else {
          if (index == 1) return -((log(fabs(tb))-log(fabs(tm)))*(log(fabs(tb))-log(fabs(tp))))/(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp))));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return ((log(fabs(tb))-log(fabs(tm)))*(log(fabs(tb))-log(fabs(tp))))/(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp))));
        }
      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return (tb*(1.0/2.0)-tp*(1.0/2.0))/(mp2*yb);
        if (index == 4) return (tb*(1.0/2.0)-tp*(1.0/2.0))/(mp2*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return 1.0/mp2-(tp*(log(fabs(tb))-log(fabs(tp))))/(mp2*(tb-tp));
        if (index == 4) return -1.0/mp2+(tb*(log(fabs(tb))-log(fabs(tp))))/(mp2*(tb-tp));
      }
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp1*yb);
        if (index == 2) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp1*yb);
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return -1.0/mp1+(tb*(log(fabs(tb))-log(fabs(tm))))/(mp1*(tb-tm));
        if (index == 2) return 1.0/mp1-(tm*(log(fabs(tb))-log(fabs(tm))))/(mp1*(tb-tm));
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return (1.0/(tb*tb)*(tb*(tp*tp)*4.0-(tb*tb)*tp*3.0+tb*tb*tb-(tp*tp*tp)*2.0)*(1.0/2.0))/(mp2*yb);
        if (index == 4) return (tb*tp*(3.0/2.0)-(tb*tb)*(1.0/2.0)-tp*tp)/(mp2*tb*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return 1.0/mp2-(tp*(log(fabs(tb))-log(fabs(tp))))/(mp2*(tb-tp));
        if (index == 4) return -1.0/mp2+(tb*(log(fabs(tb))-log(fabs(tp))))/(mp2*(tb-tp));
      }
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return (tb*tm*(3.0/2.0)-(tb*tb)*(1.0/2.0)-tm*tm)/(mp1*tb*yb);
        if (index == 2) return (1.0/(tb*tb)*(tb*(tm*tm)*4.0-(tb*tb)*tm*3.0+tb*tb*tb-(tm*tm*tm)*2.0)*(1.0/2.0))/(mp1*yb);
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return -1.0/mp1+(1.0/(tb*tb)*((tm*tm*tm)*log(fabs(tb))*(-1.0/2.0)+(tm*tm*tm)*log(fabs(tm))*(1.0/2.0)+(tb*tb)*(tm*log(fabs(tb))*(3.0/2.0)-tm*log(fabs(tm))*(3.0/2.0))))/(mp1*(tb-tm));
        if (index == 2) return 1.0/mp1-(1.0/(tb*tb*tb)*((tm*tm*tm*tm)*log(fabs(tb))-(tm*tm*tm*tm)*log(fabs(tm))-tb*((tm*tm*tm)*log(fabs(tb))*(3.0/2.0)-(tm*tm*tm)*log(fabs(tm))*(3.0/2.0))+(tb*tb*tb)*(tm*log(fabs(tb))*(3.0/2.0)-tm*log(fabs(tm))*(3.0/2.0))))/(mp1*(tb-tm));
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    }
    // bdc 값의 에러 출력
    else {
      printf ("%s%d\n", "greens_integral bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 차원의 에러 출력
  else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

/* 그린함수의 두번째변수로의 미분과 basis function의 적분 */
double greens_integral_tau (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
  /*
  nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  2차원의 무한대에서의 점근적함수는 u (x, y) = log (x^2 + y^2)
  3차원축대칭의 무한대에서의 점근적함수는 u (r, z) = 1 / (r^2 + z^2)

  index: 그린함수와 basis function의 적분을 가리키는 변수
  1: 왼쪽(아래쪽) 점에서 1인 piecewise linear basis function
  2: 가운데 점에서 1인 tau보다 큰 부분에서는 0인 basis function
  3: 가운데 점에서 1인 tau보다 작은 부분에서는 0인 basis function
  4: 오른쪽(위쪽) 점에서 1인 piecewise linear basis function

  node: x축 성분인지 y축 성분인지를 나타내는 변수
  1: x축 성분
  2: y축 성분

  bdc: boundary condition
  0: Dirichlet - Dirichlet
  1: Neumann - Dirichlet
  2: Dirichlet - Neumann
  3: Infinity - Dirichlet
  4: Dirichlet - Infinity
  5: Singularity - Dirichlet
  6: Dirichlet - Singularity

  mp1: 왼쪽(아래쪽) 부분의 conductivity
  mp2: 오른쪽(위쪽) 부분의 conductivity
  */

  // 차원
  int nD = 2;

  // node 값의 에러 출력
  if (node > 2) {printf ("%s%d\n", "greens_integral_tau node error, node = ", node); exit (1);}

  // 2차원
  if (nD == 2) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // 양 쪽의 conductivity가 같은 경우
      if (IsEqualDouble (mp1, mp2)) {
        if (index == 1) return (pow(tb-tm,2.0)*(1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 2) return (pow(tb-tm,2.0)*(1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 3) return (pow(tb-tp,2.0)*(-1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 4) return (pow(tb-tp,2.0)*(-1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      }
      // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
      else {
        if (index == 1) return (pow(tb-tm,2.0)*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return (pow(tb-tp,2.0)*(-1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      if (index == 1) return (tb*(-1.0/2.0)+tm*(1.0/2.0))/mp2;
      if (index == 2) return (tb*(-1.0/2.0)+tm*(1.0/2.0))/mp2;
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return (tb*(-1.0/2.0)+tp*(1.0/2.0))/mp1;
      if (index == 4) return (tb*(-1.0/2.0)+tp*(1.0/2.0))/mp1;
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return (tb*(tp*tp)*(-3.0/2.0)+(tb*tb)*tp-(tb*tb*tb)*(1.0/2.0)+tp*tp*tp+tb*(tp*tp)*log(fabs(tb))-tb*(tp*tp)*log(fabs(tp)))/(mp2*tb*(tb-tp));
      if (index == 4) return (-(tp*tp)*log(fabs(tb))+(tp*tp)*log(fabs(tp))+tb*tp*2.0-(tb*tb)*(1.0/2.0)-(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      if (index == 1) return (-(tm*tm)*log(fabs(tb))+(tm*tm)*log(fabs(tm))+tb*tm*2.0-(tb*tb)*(1.0/2.0)-(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
      if (index == 2) return (tb*(tm*tm)*(-3.0/2.0)+(tb*tb)*tm-(tb*tb*tb)*(1.0/2.0)+tm*tm*tm+tb*(tm*tm)*log(fabs(tb))-tb*(tm*tm)*log(fabs(tm)))/(mp1*tb*(tb-tm));
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    }
    // bdc 값의 에러 출력
    else {
      printf ("%s%d\n", "greens_integral_tau bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 3차원 축대칭
  else if (nD == 3) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // x축 성분
      if (node == 1) {
        // 양 쪽의 conductivity가 같은 경우
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return (pow(tb-tm,2.0)*(1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return (pow(tb-tm,2.0)*(1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 3) return (pow(tb-tp,2.0)*(-1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 4) return (pow(tb-tp,2.0)*(-1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        }
        // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
        else {
          if (index == 1) return (pow(tb-tm,2.0)*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return (pow(tb-tp,2.0)*(-1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        }
      }
      // y축 성분
      else if (node == 2) {
        // 양 쪽의 conductivity가 같은 경우
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return ((tb*tb)*log(fabs(tb))*(-1.0/2.0)+(tb*tb)*log(fabs(tm))*(1.0/2.0)-tb*tm+(tb*tb)*(3.0/4.0)+(tm*tm)*(1.0/4.0))/(tb*(tb-tm)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 2) return ((tb*tb)*log(fabs(tb))*(-1.0/2.0)+(tb*tb)*log(fabs(tm))*(1.0/2.0)-tb*tm+(tb*tb)*(1.0/4.0)+(tm*tm)*(3.0/4.0)+tb*tm*log(fabs(tb))-tb*tm*log(fabs(tm)))/(tb*(tb-tm)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 3) return ((tb*tb)*log(fabs(tb))*(1.0/2.0)-(tb*tb)*log(fabs(tp))*(1.0/2.0)+tb*tp-(tb*tb)*(1.0/4.0)-(tp*tp)*(3.0/4.0)-tb*tp*log(fabs(tb))+tb*tp*log(fabs(tp)))/(tb*(tb-tp)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 4) return ((tb*tb)*log(fabs(tb))*(1.0/2.0)-(tb*tb)*log(fabs(tp))*(1.0/2.0)+tb*tp-(tb*tb)*(3.0/4.0)-(tp*tp)*(1.0/4.0))/(tb*(tb-tp)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
        }
        // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
        else {
          if (index == 1) return (tb-tm-tb*log(fabs(tb))+tb*log(fabs(tm)))/(tb*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return -(tb-tp-tb*log(fabs(tb))+tb*log(fabs(tp)))/(tb*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
        }
      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      if (node == 1) {
        if (index == 1) return (tb*(-1.0/2.0)+tm*(1.0/2.0))/(mp2*yb);
        if (index == 2) return (tb*(-1.0/2.0)+tm*(1.0/2.0))/(mp2*yb);
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
      else if (node == 2) {
        if (index == 1) return (tb*(-1.0/2.0)+tm*(1.0/2.0))/(mp2*tb);
        if (index == 2) return (tb*(-1.0/2.0)+tm*(1.0/2.0))/(mp2*tb);
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return (tb*(-1.0/2.0)+tp*(1.0/2.0))/(mp1*yb);
        if (index == 4) return (tb*(-1.0/2.0)+tp*(1.0/2.0))/(mp1*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return (tb*(-1.0/2.0)+tp*(1.0/2.0))/(mp1*tb);
        if (index == 4) return (tb*(-1.0/2.0)+tp*(1.0/2.0))/(mp1*tb);
      }
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return (tb*(tp*tp)*(-3.0/2.0)+(tb*tb)*tp-(tb*tb*tb)*(1.0/2.0)+tp*tp*tp+tb*(tp*tp)*log(fabs(tb))-tb*(tp*tp)*log(fabs(tp)))/(mp2*tb*yb*(tb-tp));
        if (index == 4) return (-(tp*tp)*log(fabs(tb))+(tp*tp)*log(fabs(tp))+tb*tp*2.0-(tb*tb)*(1.0/2.0)-(tp*tp)*(3.0/2.0))/(mp2*yb*(tb-tp));
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return (-(tm*tm)*log(fabs(tb))+(tm*tm)*log(fabs(tm))+tb*tm*2.0-(tb*tb)*(1.0/2.0)-(tm*tm)*(3.0/2.0))/(mp1*yb*(tb-tm));
        if (index == 2) return (tb*(tm*tm)*(-3.0/2.0)+(tb*tb)*tm-(tb*tb*tb)*(1.0/2.0)+tm*tm*tm+tb*(tm*tm)*log(fabs(tb))-tb*(tm*tm)*log(fabs(tm)))/(mp1*tb*yb*(tb-tm));
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return (1.0/(tb*tb)*pow(tb-tm,2.0)*(-1.0/2.0))/mp1;
        if (index == 2) return (1.0/(tb*tb*tb)*(tb+tm)*pow(tb-tm,2.0)*(-1.0/2.0))/mp1;
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    }
    // bdc 값의 에러 출력
    else {
      printf ("%s%d\n", "greens_integral_tau bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 차원의 에러 출력
  else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

/* 그린함수의 첫번째, 두번째변수로의 미분과 basis function의 적분 */
double greens_integral_ttau (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
  /*
  nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  2차원의 무한대에서의 점근적함수는 u (x, y) = log (x^2 + y^2)
  3차원축대칭의 무한대에서의 점근적함수는 u (r, z) = 1 / (r^2 + z^2)

  index: 그린함수와 basis function의 적분을 가리키는 변수
  1: 왼쪽(아래쪽) 점에서 1인 piecewise linear basis function
  2: 가운데 점에서 1인 tau보다 큰 부분에서는 0인 basis function
  3: 가운데 점에서 1인 tau보다 작은 부분에서는 0인 basis function
  4: 오른쪽(위쪽) 점에서 1인 piecewise linear basis function

  node: x축 성분인지 y축 성분인지를 나타내는 변수
  1: x축 성분
  2: y축 성분

  bdc: boundary condition
  0: Dirichlet - Dirichlet
  1: Neumann - Dirichlet
  2: Dirichlet - Neumann
  3: Infinity - Dirichlet
  4: Dirichlet - Infinity
  5: Singularity - Dirichlet
  6: Dirichlet - Singularity

  mp1: 왼쪽(아래쪽) 부분의 conductivity
  mp2: 오른쪽(위쪽) 부분의 conductivity
  */

  // 차원
  int nD = 2;

  // node 값의 에러 출력
  if (node > 2) {printf ("%s%d\n", "greens_integral_tau node error, node = ", node); exit (1);}

  // 2차원
  if (nD == 2) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // 양 쪽의 conductivity가 같은 경우
      if (IsEqualDouble (mp1, mp2)) {
        if (index == 1) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 2) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 3) return (tb*(-1.0/2.0)+tp*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 4) return (tb*(-1.0/2.0)+tp*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      }
      // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
      else {
        if (index == 1) return (tb-tm)/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return -(tb-tp)/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return -(1.0/(tb*tb)*tp*(tb-tp))/mp2;
      if (index == 4) return -(tb-tp)/(mp2*tb);
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      if (index == 1) return -(tb-tm)/(mp1*tb);
      if (index == 2) return -(1.0/(tb*tb)*tm*(tb-tm))/mp1;
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    }
    // bdc 값의 에러 출력
    else {
      printf ("%s%d\n", "greens_integral_tau bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 3차원 축대칭
  else if (nD == 3) {
    // Dirichlet - Dirichlet
    if (bdc == 0) {
      // x축 성분
      if (node == 1) {
        // 양 쪽의 conductivity가 같은 경우
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 3) return (tb*(-1.0/2.0)+tp*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 4) return (tb*(-1.0/2.0)+tp*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        }
        // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
        else {
          if (index == 1) return (tb-tm)/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return -(tb-tp)/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        }
      }
      // y축 성분
      else if (node == 2) {
        // 양 쪽의 conductivity가 같은 경우
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return (tb-tm-tb*log(fabs(tb))+tb*log(fabs(tm)))/(tb*(tb-tm)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 2) return -(tb-tm-tm*log(fabs(tb))+tm*log(fabs(tm)))/(tb*(tb-tm)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 3) return (tb-tp-tp*log(fabs(tb))+tp*log(fabs(tp)))/(tb*(tb-tp)*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 4) return -(tp-tb*(-log(fabs(tb))+log(fabs(tp))+1.0))/(tb*(tb-tp)*(mp1*log(fabs(tb))-mp2*log(fabs(tb))+mp2*log(fabs(tm))-mp1*log(fabs(tp))));
        }
        // 양 쪽의 conductivity가 다른 경우(Interface 위의 점)
        else {
          if (index == 1) return -(log(fabs(tb))-log(fabs(tm)))/(tb*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return (log(fabs(tb))-log(fabs(tp)))/(tb*(mp2*(log(fabs(tb))-log(fabs(tm)))-mp1*(log(fabs(tb))-log(fabs(tp)))));
        }
      }
    }
    // Neumann - Dirichlet
    else if (bdc == 1) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    }
    // Dirichlet - Neumann
    else if (bdc == 2) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    }
    // Infinity - Dirichlet
    else if (bdc == 3) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return -(1.0/(tb*tb)*tp*(tb-tp))/(mp2*yb);
        if (index == 4) return -(tb-tp)/(mp2*tb*yb);
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    }
    // Dirichlet - Infinity
    else if (bdc == 4) {
      // x축 성분
      if (node == 1) {
        if (index == 1) return -(tb-tm)/(mp1*tb*yb);
        if (index == 2) return -(1.0/(tb*tb)*tm*(tb-tm))/(mp1*yb);
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
      // y축 성분
      else if (node == 2) {
        if (index == 1) return (1.0/(tb*tb*tb)*(tb*tm-(tb*tb)*2.0+tm*tm)*(1.0/2.0))/mp1;
        if (index == 2) return (1.0/(tb*tb*tb*tb)*tm*(tb*tm+tb*tb-(tm*tm)*2.0)*(-1.0/2.0))/mp1;
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    }
    // bdc 값의 에러 출력
    else {
      printf ("%s%d\n", "greens_integral_tau bdc error, bdc = ", bdc);
      exit (1);
    }
  }
  // 차원의 에러 출력
  else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

/* Solution의 Representation formula의 경계값의 계수 */
double greens_coefficient_t (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
  /*
  nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  2차원의 무한대에서의 점근적함수는 u (x, y) = log (x^2 + y^2)
  3차원축대칭의 무한대에서의 점근적함수는 u (r, z) = 1 / (r^2 + z^2)

  eps: conductivity를 저장하는 변수

  node: x축 성분인지 y축 성분인지를 나타내는 변수
  1: x축 성분
  2: y축 성분

  bdc: boundary condition
  0: Dirichlet - Dirichlet
  1: Neumann - Dirichlet
  2: Dirichlet - Neumann
  3: Infinity - Dirichlet
  4: Dirichlet - Infinity
  5: Singularity - Dirichlet
  6: Dirichlet - Singularity

  mp1: 왼쪽(아래쪽) 부분의 conductivity
  mp2: 오른쪽(위쪽) 부분의 conductivity
  */

  // conductivity
  double eps = 0.0;

  // 차원
  int nD = 2;

  // node 값의 에러 출력
  if (node > 2) {printf ("%s%d\n", "greens_coefficient_t node error, node = ", node); exit (1);}

  // 2차원
  if (nD == 2) {
    if (t < tb) eps = mp1;
    else        eps = mp2;
  }
  // 3차원 축대칭
  else if (nD == 3) {
    // x축 성분
    if ( node == 1 ) {
      if (t < tb) eps = mp1 * yb;
      else        eps = mp2 * yb;
    }
    // y축 성분
    else if ( node == 2 ) {
      if (t < tb) eps = mp1 * t;
      else        eps = mp2 * t;
    }
  }
  // 차원의 에러 출력
  else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }

  // Dirichlet - Dirichlet
  if (bdc == 0) {
    if (IsEqualDouble (t, tm)) return   eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    if (t < tb) return   eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return - eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
  }
  // Neumann - Dirichlet
  else if (bdc == 1) {
    if (IsEqualDouble (t, tm)) return - eps * greens_function (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    if (t < tb) return - eps * greens_function (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   1.0E0;
  }
  // Dirichlet - Neumann
  else if (bdc == 2) {
    if (IsEqualDouble (t, tm)) return 1.0E0;
    if (t < tb) return   1.0E0;
    else        return   eps * greens_function (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
  }
  // Infinity - Dirichlet
  else if (bdc == 3) {
    if (IsEqualDouble (t, tm)) return   0.0;
    if (t < tb) return   0.0;
    else        return - eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
  }
  // Dirichlet - Infinity
  else if (bdc == 4) {
    if (IsEqualDouble (t, tm)) return   eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    if (t < tb) return   eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   0.0;
  }
  // bdc값의 에러 출력
  else {
    printf ("greens_coefficient_t bdc error, bdc = %d\n", bdc);
    exit (1);
  }
}

/* 미분의 Representation formula의 경계값의 계수 */
double greens_coefficient_ttau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
  /*
  nD: 2차원 혹운 3차원 축대칭 (2: 2차원, 3: 3차원 축대칭)
  2차원의 무한대에서의 점근적함수는 u (x, y) = log (x^2 + y^2)
  3차원축대칭의 무한대에서의 점근적함수는 u (r, z) = 1 / (r^2 + z^2)

  eps: conductivity를 저장하는 변수

  node: x축 성분인지 y축 성분인지를 나타내는 변수
  1: x축 성분
  2: y축 성분

  bdc: boundary condition
  0: Dirichlet - Dirichlet
  1: Neumann - Dirichlet
  2: Dirichlet - Neumann
  3: Infinity - Dirichlet
  4: Dirichlet - Infinity
  5: Singularity - Dirichlet
  6: Dirichlet - Singularity

  mp1: 왼쪽(아래쪽) 부분의 conductivity
  mp2: 오른쪽(위쪽) 부분의 conductivity
  */

  // conductivity
  double eps = 0.0;

  // 차원
  int nD = 2;

  // 2차원
  if (nD == 2) {
    if (t < tb) eps = mp1;
    else        eps = mp2;
  }
  // 3차원 축대칭
  else if (nD == 3) {
    // x축 성분
    if ( node == 1 ) {
      if (t < tb) eps = mp1 * yb;
      else        eps = mp2 * yb;
    }
    // y축 성분
    else if ( node == 2 ) {
      if (t < tb) eps = mp1 * t;
      else        eps = mp2 * t;
    }
  }
  // 차원의 에러 출력
  else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }

  // Dirichlet - Dirichlet
  if (bdc == 0) {
    if (IsEqualDouble (t, tm)) return   eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    if (t < tb) return   eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return - eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
  }
  // Neumann - Dirichlet
  else if (bdc == 1) {
    if (IsEqualDouble (t, tm)) return   eps * greens_function_tau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    if (t < tb) return   eps * greens_function_tau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return - eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
  }
  // Dirichlet - Neumann
  else if (bdc == 2) {
    if (IsEqualDouble (t, tm)) return   eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    if (t < tb) return   eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   eps * greens_function_tau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
  }
  // Infinity - Dirichlet
  else if (bdc == 3) {
    if (IsEqualDouble (t, tm)) return   0.0;
    if (t < tb) return   0.0;
    else        return - eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
  }
  // Dirichlet - Infinity
  else if (bdc == 4) {
    if (IsEqualDouble (t, tm)) return   eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    if (t < tb) return   eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   0.0;
  }
  // bdc값의 에러 출력
  else {
    printf ("greens_coefficient_ttau bdc error, bdc = %d\n", bdc);
    exit (1);
  }
}

#endif
