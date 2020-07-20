#ifndef FTNS_H
#define FTNS_H

#include "Read.hpp"

// exact solution u
double u_ftn(double x, double y) {
    double t = 1.1;
    return t + pow(x, 2) + pow(y, 2);

}

double u_ftn(double x, double y, double t) {
    return exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2));

}

double u_ftn(Coordinate src) {
    double x = src[0], y = src[1];

    return u_ftn(x, y);
}

double u_ftn(Point *pt) {
    char errorMassage[256];
    // double eps = pt->Pressure ()->MaterialProperty ();
    double x = pt->Pressure()->Coord().Value('x');
    double y = pt->Pressure()->Coord().Value('y');
    double t = pt->Pressure()->Time();
    double dt = pt->Pressure()->Dt();

    // double gf = exp (- (x * x + y * y) / (2.0E0 * t * t));
    // double spi = sqrt (2.0E0 * M_PI);

    if (pt->Mark() == "u") {
        if (pt->Condition() == 'N')
            return u_ftn(pt->Diff('x')) * pt->Pressure()->Normal()[0] +
                   u_ftn(pt->Diff('y')) * pt->Pressure()->Normal()[1];
        else return exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2));
        // else                         return sin (x * y);
    }
    if (pt->Mark() =="v") {
        if (pt->Condition() == 'N')
            return u_ftn(pt->Diff('x')) * pt->Pressure()->Normal()[0] +
                   u_ftn(pt->Diff('y')) * pt->Pressure()->Normal()[1];
        else return exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2));
    }
    if (pt->Mark() == "ux") return -x * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 2);
    if (pt->Mark() == "uy") return -y * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 2);
    if (pt->Mark() == "vx") return -x * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 2);
    if (pt->Mark() == "vy") return -y * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 2);
    if (pt->Mark() == "uxx")
        return (-pow(t, 2) + pow(x, 2))*exp(-1.0/2.0*(pow(x, 2) + pow(y, 2))/pow(t, 2))/pow(t, 4);
    if (pt->Mark() == "uyy")
        return (-pow(t, 2) + pow(y, 2))*exp(-1.0/2.0*(pow(x, 2) + pow(y, 2))/pow(t, 2))/pow(t, 4);
    if (pt->Mark() == "vxx")
        return -exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 2) +
               pow(x, 2) * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 4);
    if (pt->Mark() == "vyy")
        return -exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 2) +
               pow(y, 2) * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 4);
    if (pt->Mark() == "phi")
        return 0.25 * (-pow(x, 2) + pow(y, 2)) *
               (pow(t, 4) * exp(-1.0 / 2.0 * (pow(x, 2) + pow(y, 2)) / pow(dt - t, 2)) +
                pow(dt - t, 4) * exp(-1.0 / 2.0 * (pow(x, 2) + pow(y, 2)) / pow(t, 2))) / (pow(t, 4) * pow(dt - t, 4));
    if (pt->Mark() == "psi")
        return 0.25 * (-pow(x, 2) + pow(y, 2)) *
               (pow(t, 4) * exp(-1.0 / 2.0 * (pow(x, 2) + pow(y, 2)) / pow(dt - t, 2)) +
                pow(dt - t, 4) * exp(-1.0 / 2.0 * (pow(x, 2) + pow(y, 2)) / pow(t, 2))) / (pow(t, 4) * pow(dt - t, 4));
    if (pt->Mark() == "diff_u")
        return sqrt(u_ftn(pt->Velocity('u')->Diff('x')) * u_ftn(pt->Velocity('u')->Diff('x')) +
                    u_ftn(pt->Velocity('u')->Diff('y')) * u_ftn(pt->Velocity('u')->Diff('y')));
    if (pt->Mark() == "diff_v")
        return sqrt(u_ftn(pt->Velocity('v')->Diff('x')) * u_ftn(pt->Velocity('v')->Diff('x')) +
                    u_ftn(pt->Velocity('v')->Diff('y')) * u_ftn(pt->Velocity('v')->Diff('y')));
    if (pt->Mark() == "diff_diff_u")
        return sqrt(u_ftn(pt->Velocity('u')->Diff('x')->Diff('x')) * u_ftn(pt->Velocity('u')->Diff('x')->Diff('x')) +
                    u_ftn(pt->Velocity('u')->Diff('y')->Diff('y')) * u_ftn(pt->Velocity('u')->Diff('y')->Diff('y')));
    if (pt->Mark() == "diff_diff_v")
        return sqrt(u_ftn(pt->Velocity('v')->Diff('x')->Diff('x')) * u_ftn(pt->Velocity('v')->Diff('x')->Diff('x')) +
                    u_ftn(pt->Velocity('v')->Diff('y')->Diff('y')) * u_ftn(pt->Velocity('v')->Diff('y')->Diff('y')));

    sprintf(errorMassage, "u_ftn (Point *pt), pt->Mark () == %s, please check Mark.", pt->Mark().c_str());
    PrintError(errorMassage);

    exit(1);
}

double u_ftn_Dirichlet(Point *pt) {
    char errorMassage[256];
    double x = pt->Pressure()->Coord().Value('x');
    double y = pt->Pressure()->Coord().Value('y');
    double t = pt->Pressure()->Time();

    // double gf = exp (- (x * x + y * y) / (2.0E0 * t * t));
    // double spi = sqrt (2.0E0 * M_PI);

    // if (!pt->Mark ().compare ("u"))  return sin (x * y);
    if (pt->Mark() == "u") return exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2));
    if (pt->Mark() == "v") return exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2));

    if (pt->Mark() == "ux") return u_ftn(pt);
    if (pt->Mark() == "uy") return u_ftn(pt);
    if (pt->Mark() == "vx") return u_ftn(pt);
    if (pt->Mark() == "vy") return u_ftn(pt);
    if (pt->Mark() == "uxx") return u_ftn(pt);
    if (pt->Mark() == "uyy") return u_ftn(pt);
    if (pt->Mark() == "vxx") return u_ftn(pt);
    if (pt->Mark() == "vyy") return u_ftn(pt);
    if (pt->Mark() == "phi") return u_ftn(pt);
    if (pt->Mark() == "psi") return u_ftn(pt);
    if (pt->Mark() == "diff_u") return u_ftn(pt);
    if (pt->Mark() == "diff_v") return u_ftn(pt);
    if (pt->Mark() == "diff_diff_u") return u_ftn(pt);
    if (pt->Mark() == "diff_diff_v") return u_ftn(pt);

    sprintf(errorMassage, "u_ftn_Dirichlet (Point *pt), pt->Mark () == %s, please check Mark.", pt->Mark().c_str());
    PrintError(errorMassage);
    exit(1);
}

// differentiation with respect to x
double dudx_ftn(double x, double y) {

    return 2.0E0 * x;

    // ux (x, y) = 1/(x^2 + y^2) - (2*x^2)/(x^2 + y^2)^2
    // ux (x, y) = x/(x^2 + y^2)
    return x / (x * x + y * y);
    // return 1.0;

}

// differentiation with respect to y
double dudy_ftn(double x, double y) {

    return -2.0E0 * y;

    // uy (x, y) = 1/(x^2 + y^2) - (2*y^2)/(x^2 + y^2)^2
    // uy (x, y) = y/(x*x+y*y)
    return y / (x * x + y * y);
    // return -y*1.0/pow(x*x+y*y,3.0/2.0);
    // return 0.0;

}

// Dirichlet 경계값을 주기위한 함수
double b_u_ftn(double x, double y) {

    // 경계값으로 exact solution을 준다.
    return u_ftn(x, y);

}

double f_ftn(double x, double y, double t) {
    return 2.0 * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 2) -
           (-pow(x, 2) - pow(y, 2)) * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 3) -
           1.0 * pow(x, 2) * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 4) -
           1.0 * pow(y, 2) * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 4);

}

// External force f
double f_ftn(Point *pt) {
    double eps = pt->Pressure()->MaterialProperty();
    double x = pt->Pressure()->Coord().Value('x');
    double y = pt->Pressure()->Coord().Value('y');
    double t = pt->Pressure()->Time();
    // double gf = exp (- (x * x + y * y) / (2.0E0 * t * t));
    // double spi = sqrt (2.0E0 * M_PI);
    double dt = pt->Pressure()->Dt();

    return (1.0 * pow(t, 4) * pow(dt - t, 2) * exp(-1.0 / 2.0 * (pow(x, 2) + pow(y, 2)) / pow(dt - t, 2)) -
            0.5 * pow(t, 4) * (dt - t) * (pow(x, 2) + pow(y, 2)) *
            exp(-1.0 / 2.0 * (pow(x, 2) + pow(y, 2)) / pow(dt - t, 2)) -
            0.5 * pow(t, 4) * (pow(x, 2) + pow(y, 2)) * exp(-1.0 / 2.0 * (pow(x, 2) + pow(y, 2)) / pow(dt - t, 2)) +
            1.0 * pow(t, 2) * pow(dt - t, 4) * exp(-1.0 / 2.0 * (pow(x, 2) + pow(y, 2)) / pow(t, 2)) +
            0.5 * t * pow(dt - t, 4) * (pow(x, 2) + pow(y, 2)) * exp(-1.0 / 2.0 * (pow(x, 2) + pow(y, 2)) / pow(t, 2)) -
            0.5 * pow(dt - t, 4) * (pow(x, 2) + pow(y, 2)) * exp(-1.0 / 2.0 * (pow(x, 2) + pow(y, 2)) / pow(t, 2))) /
           (pow(t, 4) * pow(dt - t, 4));

    return 5.0E-1 * (f_ftn(x, y, t) + f_ftn(x, y, t - pt->Pressure()->Dt()));

    printf("Mark = %s\n", pt->Mark().c_str());
    exit(1);
}

double eps_ftn(double x, double y) {

    return 1.0;

}

double phi_ftn(double x, double y, double t) {

    return -0.5 * pow(x, 2) * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 4) +
           0.5 * pow(y, 2) * exp((1.0 / 2.0) * (-pow(x, 2) - pow(y, 2)) / pow(t, 2)) / pow(t, 4);
}

double phi_ftn(double z, double r) {
    double x = z, y = r;
    double R = sqrt(x * x + y * y);
    double Theta = atan2(y, x);

    return -6.0E0 * x;

    // return -1.0/pow(x*x+y*y,5.0/2.0)*((x*x)*2.0E2-(y*y)*1.0E2);
    return y * 1.0 / pow(x * x + y * y, 5.0 / 2.0) * ((x * x) * 2.0 - y * y) * -1.0E2;
    // return -1.0/pow(x*x+y*y,5.0/2.0)*((x*x)*y*2.0-y*y*y);

    return 1.0 / pow(x * x + y * y, 3.0) * ((x * x) * 3.0 - y * y) * -2.0;
    return y * 1.0 / pow(x * x + y * y, 3.0) * ((x * x) * 3.0 - y * y) * -2.0;

    if (IsEqualDouble(x, 0.0) && IsEqualDouble(y, 0.0)) return 0.0;

    if (Theta < 0.0) Theta += 2.0 * PI;

    return 2.0 / 9.0 * pow(R, -4.0 / 3.0) * sin(4.0 / 3.0 * Theta);

    // return -(r*1.0/pow(r*r+z*z,3.0/2.0)*-2.0+(r*r*r)*1.0/pow(r*r+z*z,5.0/2.0)*3.0);
    return 1.0 / pow(r * r + z * z, 5.0 / 2.0) * (r * (z * z) * 2.0 - r * r * r);

}

double sol_x_approx(const double xm, const double xb, const double xp, const double ym, const double yb, const double yp,
             const double T, const double Dt) {

    double value = ZeroValue;
    double dt = Dt;
    double VALUE[5] = {ZeroValue,};

    for (size_t i = 0; i < 5; i++) {
        value = ZeroValue;

        dt = Dt * pow(5.0E-1, i);

        value += greens_coefficient_t(xp, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(xp, yb, T);
        value += greens_coefficient_t(xm, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(xm, yb, T);

        value += greens_integral(1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn(xm, yb, T);
        value += greens_integral(2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn(xb, yb, T);
        value += greens_integral(3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn(xb, yb, T);
        value += greens_integral(4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn(xp, yb, T);

        value += greens_integral(1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(xm, yb, T) * 5.0E-1;
        value += greens_integral(2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(xb, yb, T) * 5.0E-1;
        value += greens_integral(3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(xb, yb, T) * 5.0E-1;
        value += greens_integral(4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(xp, yb, T) * 5.0E-1;

        value -= greens_integral(1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(xm, yb, T) / dt * 5.0E-1;
        value -= greens_integral(2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(xb, yb, T) / dt * 5.0E-1;
        value -= greens_integral(3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(xb, yb, T) / dt * 5.0E-1;
        value -= greens_integral(4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(xp, yb, T) / dt * 5.0E-1;

        value += greens_integral(1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(xm, yb, T - dt) / dt * 5.0E-1;
        value += greens_integral(2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(xb, yb, T - dt) / dt * 5.0E-1;
        value += greens_integral(3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(xb, yb, T - dt) / dt * 5.0E-1;
        value += greens_integral(4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(xp, yb, T - dt) / dt * 5.0E-1;

        printf("value = %23.16e; ", fabs(u_ftn(xb, yb, T) - value));

        VALUE[i] = fabs(u_ftn(xb, yb, T) - value);

    }
    printf("\n");

    double MinOrder = VALUE[0] / VALUE[1];
    double Order = VALUE[0] / VALUE[1];

    for (size_t i = 1; i < 4; i++) {
        Order = VALUE[i] / VALUE[i + 1];
        if (MinOrder > Order) MinOrder = Order;
    }



    // printf ("1 = %23.16e\n", greens_integral (1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0));
    // printf ("2 = %23.16e\n", greens_integral (2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0));
    // printf ("3 = %23.16e\n", greens_integral (3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0));
    // printf ("4 = %23.16e\n", greens_integral (4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0));


    // printf ("11 = %23.16e\n", greens_function_t (xp, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0));
    // printf ("22 = %23.16e\n", greens_function_t (xm, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0));
    //
    // printf ("33 = %23.16e\n", u_ftn (xm, yb, T));
    //
    //
    // printf ("greens_coefficient_t (xp, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xp, yb, T) = %23.16e\n", greens_coefficient_t (xp, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xp, yb, T));
    // printf ("greens_coefficient_t (xm, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xm, yb, T) = %23.16e\n", greens_coefficient_t (xm, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xm, yb, T));
    // printf ("greens_integral (1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn (xm, yb, T) = %23.16e\n", greens_integral (1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn (xm, yb, T));
    // printf ("greens_integral (2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn (xb, yb, T) = %23.16e\n", greens_integral (2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn (xb, yb, T));
    // printf ("greens_integral (3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn (xb, yb, T) = %23.16e\n", greens_integral (3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn (xb, yb, T));
    // printf ("greens_integral (4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn (xp, yb, T) = %23.16e\n", greens_integral (4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * phi_ftn (xp, yb, T));
    // printf ("greens_integral (1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn (xm, yb, T) * 5.0E-1 = %23.16e\n", greens_integral (1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn (xm, yb, T) * 5.0E-1);
    // printf ("greens_integral (2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn (xb, yb, T) * 5.0E-1 = %23.16e\n", greens_integral (2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn (xb, yb, T) * 5.0E-1);
    // printf ("greens_integral (3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn (xb, yb, T) * 5.0E-1 = %23.16e\n", greens_integral (3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn (xb, yb, T) * 5.0E-1);
    // printf ("greens_integral (4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn (xp, yb, T) * 5.0E-1 = %23.16e\n", greens_integral (4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn (xp, yb, T) * 5.0E-1);
    // printf ("greens_integral (1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xm, yb, T) / Dt * 5.0E-1 = %23.16e\n", greens_integral (1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xm, yb, T) / Dt * 5.0E-1);
    // printf ("greens_integral (2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xb, yb, T) / Dt * 5.0E-1 = %23.16e\n", greens_integral (2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xb, yb, T) / Dt * 5.0E-1);
    // printf ("greens_integral (3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xb, yb, T) / Dt * 5.0E-1 = %23.16e\n", greens_integral (3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xb, yb, T) / Dt * 5.0E-1);
    // printf ("greens_integral (4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xp, yb, T) / Dt * 5.0E-1 = %23.16e\n", greens_integral (4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xp, yb, T) / Dt * 5.0E-1);
    // printf ("greens_integral (1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xm, yb, T - Dt) / Dt * 5.0E-1 = %23.16e\n", greens_integral (1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xm, yb, T - Dt) / Dt * 5.0E-1);
    // printf ("greens_integral (2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xb, yb, T - Dt) / Dt * 5.0E-1 = %23.16e\n", greens_integral (2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xb, yb, T - Dt) / Dt * 5.0E-1);
    // printf ("greens_integral (3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xb, yb, T - Dt) / Dt * 5.0E-1 = %23.16e\n", greens_integral (3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xb, yb, T - Dt) / Dt * 5.0E-1);
    // printf ("greens_integral (4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xp, yb, T - Dt) / Dt * 5.0E-1 = %23.16e\n", greens_integral (4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn (xp, yb, T - Dt) / Dt * 5.0E-1);


    return MinOrder;

}

double sol_x_approx(Point *pt) {

    double value = -1.0E0 * u_ftn(pt);
    double xm = pt->MinMaxCoordinate('x', 'm'), xb = pt->Coord().Value('x'), xp = pt->MinMaxCoordinate('x', 'p');
    double yb = pt->Coord().Value('y');
    double dt = pt->Dt();
    double ent_x = ZeroValue;
    double uu = ZeroValue;

    value += greens_coefficient_t(xm, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('W', 'W'));
    value += greens_coefficient_t(xp, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('E', 'E'));

    value += greens_integral(1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('W', 'W')->Phi());
    value += greens_integral(2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(pt->Phi());
    value += greens_integral(3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(pt->Phi());
    value += greens_integral(4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('E', 'E')->Phi());

    value += greens_integral(1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(pt->EWNS('W', 'W')) * 5.0E-1;
    value += greens_integral(2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(pt) * 5.0E-1;
    value += greens_integral(3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(pt) * 5.0E-1;
    value += greens_integral(4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(pt->EWNS('E', 'E')) * 5.0E-1;

    value -= greens_integral(1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('W', 'W')) / dt * 5.0E-1;
    value -= greens_integral(2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(pt) / dt * 5.0E-1;
    value -= greens_integral(3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(pt) / dt * 5.0E-1;
    value -= greens_integral(4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('E', 'E')) / dt * 5.0E-1;

    value += greens_integral(1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * pt->EWNS('W', 'W')->Pre()->Value() / dt *
             5.0E-1;
    value += greens_integral(2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * pt->Pre()->Value() / dt * 5.0E-1;
    value += greens_integral(3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * pt->Pre()->Value() / dt * 5.0E-1;
    value += greens_integral(4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * pt->EWNS('E', 'E')->Pre()->Value() / dt *
             5.0E-1;

    ent_x = -1.0E0;
    ent_x -= greens_integral(2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) / dt * 5.0E-1;
    ent_x -= greens_integral(3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) / dt * 5.0E-1;
    uu = u_ftn(pt);
    printf("i = 00, "), printf("value = %23.16e, ", ent_x * uu), printf("ent_x = %23.16e, ", ent_x), printf(
            "u_ftn = %23.16e\n", uu);

    ent_x = greens_coefficient_t(xp, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0);
    ent_x -= greens_integral(4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) / dt * 5.0E-1;
    uu = u_ftn(pt->EWNS('E', 'E'));
    printf("i = 01, "), printf("value = %23.16e, ", ent_x * uu), printf("ent_x = %23.16e, ", ent_x), printf(
            "u_ftn = %23.16e\n", uu);

    ent_x = greens_coefficient_t(xm, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0);
    ent_x -= greens_integral(1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) / dt * 5.0E-1;
    uu = u_ftn(pt->EWNS('W', 'W'));
    printf("i = 02, "), printf("value = %23.16e, ", ent_x * uu), printf("ent_x = %23.16e, ", ent_x), printf(
            "u_ftn = %23.16e\n", uu);

    ent_x = greens_integral(2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) +
            greens_integral(3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0);
    uu = u_ftn(pt->Phi());
    printf("i = 07, "), printf("value = %23.16e, ", ent_x * uu), printf("ent_x = %23.16e, ", ent_x), printf(
            "u_ftn = %23.16e\n", uu);

    ent_x = greens_integral(4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0);
    uu = u_ftn(pt->EWNS('E', 'E')->Phi());
    printf("i = 08, "), printf("value = %23.16e, ", ent_x * uu), printf("ent_x = %23.16e, ", ent_x), printf(
            "u_ftn = %23.16e\n", uu);

    ent_x = greens_integral(1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0);
    uu = u_ftn(pt->EWNS('W', 'W')->Phi());
    printf("i = 09, "), printf("value = %23.16e, ", ent_x * uu), printf("ent_x = %23.16e, ", ent_x), printf(
            "u_ftn = %23.16e\n", uu);

    uu = ZeroValue;

    uu += greens_integral(1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(pt->EWNS('W', 'W')) * 5.0E-1;
    uu += greens_integral(2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(pt) * 5.0E-1;
    uu += greens_integral(3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(pt) * 5.0E-1;
    uu += greens_integral(4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * f_ftn(pt->EWNS('E', 'E')) * 5.0E-1;

    // uu += greens_integral (1, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * pt->EWNS ('W', 'W')->Pre ()->Value () / dt * 5.0E-1;
    // uu += greens_integral (2, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * pt                 ->Pre ()->Value () / dt * 5.0E-1;
    // uu += greens_integral (3, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * pt                 ->Pre ()->Value () / dt * 5.0E-1;
    // uu += greens_integral (4, xm, xb, xp, xb, yb, 1, 0, 1.0E0, 1.0E0) * pt->EWNS ('E', 'E')->Pre ()->Value () / dt * 5.0E-1;

    printf("xdat->F = %23.16e\n", uu);

    // pt->SetValue (value);
    return value;
}

double sol_diff_approx(Point *pt) {
    double value = -pt->Boundaryvalue();
    double xm = pt->MinMaxCoordinate('x', 'm'), xb = pt->Coord().Value('x'), xp = pt->MinMaxCoordinate('x', 'p');
    double ym = pt->MinMaxCoordinate('y', 'm'), yb = pt->Coord().Value('y'), yp = pt->MinMaxCoordinate('y', 'p');
    double dt = pt->Dt();
    double ent_x[14] = {ZeroValue,}, ent_y[14] = {ZeroValue,};
    double xf = ZeroValue, yf = ZeroValue;
    double yen, yes, ywn, yws, xne, xnw, xse, xsw;
    double epsilon = 5.0E-1;
    Point *func_x[14], *func_y[14];
    string azimuth_x[6] = {"EE", "WW", "EN", "ES", "WN", "WS"};
    string azimuth_y[6] = {"NN", "SS", "NE", "NW", "SE", "SW"};

    func_x[0] = pt;
    func_y[0] = pt;
    func_x[7] = pt->Phi();
    func_y[7] = pt->Phi();

    for (size_t i = 1; i < 7; i++) {
        if (pt->EWNS(azimuth_x[i - 1][0], azimuth_x[i - 1][1]))
            func_x[i] = pt->EWNS(azimuth_x[i - 1][0], azimuth_x[i - 1][1]);
        else func_x[i] = nullptr;
        if (pt->EWNS(azimuth_y[i - 1][0], azimuth_y[i - 1][1]))
            func_y[i] = pt->EWNS(azimuth_y[i - 1][0], azimuth_y[i - 1][1]);
        else func_y[i] = nullptr;
        if (pt->EWNS(azimuth_x[i - 1][0], azimuth_x[i - 1][1]))
            func_x[i + 7] = pt->EWNS(azimuth_x[i - 1][0], azimuth_x[i - 1][1])->Phi();
        else func_x[i + 7] = nullptr;
        if (pt->EWNS(azimuth_y[i - 1][0], azimuth_y[i - 1][1]))
            func_y[i + 7] = pt->EWNS(azimuth_y[i - 1][0], azimuth_y[i - 1][1])->Phi();
        else func_y[i + 7] = nullptr;
    }

    if (pt->EWNS('E', 'N')) yen = pt->EWNS('E', 'N')->Coord().Value('y');
    if (pt->EWNS('E', 'S')) yes = pt->EWNS('E', 'S')->Coord().Value('y');
    if (pt->EWNS('W', 'N')) ywn = pt->EWNS('W', 'N')->Coord().Value('y');
    if (pt->EWNS('W', 'S')) yws = pt->EWNS('W', 'S')->Coord().Value('y');
    if (pt->EWNS('N', 'E')) xne = pt->EWNS('N', 'E')->Coord().Value('x');
    if (pt->EWNS('N', 'W')) xnw = pt->EWNS('N', 'W')->Coord().Value('x');
    if (pt->EWNS('S', 'E')) xse = pt->EWNS('S', 'E')->Coord().Value('x');
    if (pt->EWNS('S', 'W')) xsw = pt->EWNS('S', 'W')->Coord().Value('x');


    //    ⌠
    //    ⎮         ∂
    //    ⎮ φ(x, y)⋅──(G(x, ξ)) dx
    //    ⎮         ∂ξ
    //    ⌡

    ent_x[7] = greens_integral_tau(2, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) +
               greens_integral_tau(3, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon);
    ent_x[8] = greens_integral_tau(4, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon);
    ent_x[9] = greens_integral_tau(1, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon);

    //    ⌠
    //    ⎮          ∂
    //    ⎮ -u(x, y)⋅──(G(x, ξ))
    //    ⎮          ∂ξ
    //    ⎮ ───────────────────── dx
    //    ⎮          2⋅dt
    //    ⌡

    ent_x[0] = -5.0E-1 * (greens_integral_tau(2, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) +
                          greens_integral_tau(3, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon)) / dt;
    ent_x[1] = -5.0E-1 * greens_integral_tau(4, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) / dt;
    ent_x[2] = -5.0E-1 * greens_integral_tau(1, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) / dt;

    //              ∂
    //    ε⋅u(x, y)⋅──(G(x, ξ))
    //              ∂x

    ent_x[1] += greens_coefficient_ttau(xp, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon);
    ent_x[2] += greens_coefficient_ttau(xm, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon);

    //    ⌠
    //    ⎮         ∂
    //    ⎮ f(x, y)⋅──(G(x, ξ))
    //    ⎮         ∂ξ
    //    ⎮ ─────────────────── dx
    //    ⎮          2
    //    ⌡

    xf = 5.0E-1 * (greens_integral_tau(2, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) +
                   greens_integral_tau(3, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon)) * f_ftn(pt);

    //    ⌠
    //    ⎮          ∂
    //    ⎮ u(x, y)⋅──(G(x, ξ))
    //    ⎮          ∂ξ
    //    ⎮ ───────────────────── dx
    //    ⎮          4⋅dt
    //    ⌡

    xf += 2.5E-1 * (greens_integral_tau(2, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) +
                    greens_integral_tau(3, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon)) * pt->Pre()->Value() / dt;

    //    ⌠
    //    ⎮      ⎛  2              2         ⎞
    //    ⎮      ⎜ ∂              ∂          ⎟ ∂
    //    ⎮ 0.25⋅⎜───(u(x, y)) + ───(u(x, y))⎟⋅──(G(x, ξ)) dx
    //    ⎮      ⎜  2              2         ⎟ ∂ξ
    //    ⎮      ⎝∂x             ∂y          ⎠
    //    ⌡

    xf += 2.5E-1 * (greens_integral_tau(2, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) +
                    greens_integral_tau(3, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon)) *
          (u_ftn(HeadVelocity(pt)->Diff('x')->Diff('x')) + u_ftn(HeadVelocity(pt)->Diff('y')->Diff('y')));

    if (pt->EWNS('E', 'E'))
        xf += 5.0E-1 * greens_integral_tau(4, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) *
              f_ftn(pt->EWNS('E', 'E'))
              + 5.0E-1 * greens_integral_tau(4, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) *
                pt->EWNS('E', 'E')->Pre()->Value() / dt;

    if (pt->EWNS('W', 'W'))
        xf += 5.0E-1 * greens_integral_tau(1, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) *
              f_ftn(pt->EWNS('W', 'W'))
              + 5.0E-1 * greens_integral_tau(1, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) *
                pt->EWNS('W', 'W')->Pre()->Value() / dt;

    if (pt->EWNS('E', 'N'))
        ent_x[3] = ent_x[1] * CalcVerticalUCoefficient('N', xp, yb, yen, yes, 0, epsilon)
                   - ent_x[1] * CalcVerticalFCoefficient('N', xp, yb, yen, yes, 0, epsilon) * 5.0E-1 / dt,
        ent_x[10] = ent_x[1] * CalcVerticalPHICoefficient('N', xp, yb, yen, yes, 0, epsilon) +
                    ent_x[8] * (yb - yes) / (yen - yes),
        xf += ent_x[1] * CalcVerticalFCoefficient('N', xp, yb, yen, yes, 0, epsilon) * 5.0E-1 *
              f_ftn(pt->EWNS('E', 'N'))
              + 5.0E-1 * greens_integral_tau(4, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) * (yb - yes) /
                (yen - yes) * f_ftn(pt->EWNS('E', 'N'))
              + 5.0E-1 * greens_integral_tau(4, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) * (yb - yes) /
                (yen - yes) * pt->EWNS('E', 'N')->Pre()->Value() / dt
              + ent_x[1] * CalcVerticalFCoefficient('N', xp, yb, yen, yes, 0, epsilon) * 5.0E-1 / dt *
                pt->EWNS('E', 'N')->Pre()->Value();
    if (pt->EWNS('E', 'S'))
        ent_x[4] = ent_x[1] * CalcVerticalUCoefficient('S', xp, yb, yen, yes, 0, epsilon)
                   - ent_x[1] * CalcVerticalFCoefficient('S', xp, yb, yen, yes, 0, epsilon) * 5.0E-1 / dt,
        ent_x[11] = ent_x[1] * CalcVerticalPHICoefficient('S', xp, yb, yen, yes, 0, epsilon) +
                    ent_x[8] * (yen - yb) / (yen - yes),
        xf += ent_x[1] * CalcVerticalFCoefficient('S', xp, xb, yen, yes, 0, epsilon) * 5.0E-1 *
              f_ftn(pt->EWNS('E', 'S'))
              + 5.0E-1 * greens_integral_tau(4, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) * (yen - yb) /
                (yen - yes) * f_ftn(pt->EWNS('E', 'S'))
              + 5.0E-1 * greens_integral_tau(4, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) * (yen - yb) /
                (yen - yes) * pt->EWNS('E', 'S')->Pre()->Value() / dt
              + ent_x[1] * CalcVerticalFCoefficient('S', xp, yb, yen, yes, 0, epsilon) * 5.0E-1 / dt *
                pt->EWNS('E', 'S')->Pre()->Value();
    if (pt->EWNS('W', 'N'))
        ent_x[5] = ent_x[2] * CalcVerticalUCoefficient('N', xm, yb, ywn, yws, 0, epsilon)
                   - ent_x[2] * CalcVerticalFCoefficient('N', xm, yb, ywn, yws, 0, epsilon) * 5.0E-1 / dt,
        ent_x[12] = ent_x[2] * CalcVerticalPHICoefficient('N', xm, yb, ywn, yws, 0, epsilon) +
                    ent_x[9] * (yb - yws) / (ywn - yws),
        xf += ent_x[2] * CalcVerticalFCoefficient('N', xm, yb, ywn, yws, 0, epsilon) * 5.0E-1 *
              f_ftn(pt->EWNS('W', 'N'))
              + 5.0E-1 * greens_integral_tau(1, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) * (yb - yws) /
                (ywn - yws) * f_ftn(pt->EWNS('W', 'N'))
              + 5.0E-1 * greens_integral_tau(1, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) * (yb - yws) /
                (ywn - yws) * pt->EWNS('W', 'N')->Pre()->Value() / dt
              + ent_x[2] * CalcVerticalFCoefficient('N', xm, yb, ywn, yws, 0, epsilon) * 5.0E-1 / dt *
                pt->EWNS('W', 'N')->Pre()->Value();
    if (pt->EWNS('W', 'S'))
        ent_x[5] = ent_x[2] * CalcVerticalUCoefficient('S', xm, yb, ywn, yws, 0, epsilon)
                   - ent_x[2] * CalcVerticalFCoefficient('S', xm, yb, ywn, yws, 0, epsilon) * 5.0E-1 / dt,
        ent_x[13] = ent_x[2] * CalcVerticalPHICoefficient('S', xm, yb, ywn, yws, 0, epsilon) +
                    ent_x[9] * (ywn - yb) / (ywn - yws),
        xf += ent_x[2] * CalcVerticalFCoefficient('S', xm, yb, ywn, yws, 0, epsilon) * 5.0E-1 *
              f_ftn(pt->EWNS('W', 'S'))
              + 5.0E-1 * greens_integral_tau(1, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) * (ywn - yb) /
                (ywn - yws) * f_ftn(pt->EWNS('W', 'S'))
              + 5.0E-1 * greens_integral_tau(1, xm, xb, xp, xb, yb, 1, 0, epsilon, epsilon) * (ywn - yb) /
                (ywn - yws) * pt->EWNS('W', 'S')->Pre()->Value() / dt
              + ent_x[2] * CalcVerticalFCoefficient('S', xm, yb, ywn, yws, 0, epsilon) * 5.0E-1 / dt *
                pt->EWNS('W', 'S')->Pre()->Value();


    ent_y[0] = -5.0E-1 * (greens_integral_tau(2, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) +
                          greens_integral_tau(3, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon)) / dt;
    ent_y[1] = greens_coefficient_ttau(yp, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) -
               5.0E-1 * greens_integral_tau(4, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) / dt;
    ent_y[2] = greens_coefficient_ttau(ym, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) -
               5.0E-1 * greens_integral_tau(1, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) / dt;

    ent_y[7] = -greens_integral_tau(2, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) -
               greens_integral_tau(3, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon);
    ent_y[8] = -greens_integral_tau(4, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon);
    ent_y[9] = -greens_integral_tau(1, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon);

    yf = 5.0E-1 * (greens_integral_tau(2, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) +
                   greens_integral_tau(3, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon)) * f_ftn(pt)
         + 5.0E-1 * (greens_integral_tau(2, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) +
                     greens_integral_tau(3, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon)) * pt->Pre()->Value() / dt;

    if (pt->EWNS('N', 'N'))
        yf += 5.0E-1 * greens_integral_tau(4, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) *
              f_ftn(pt->EWNS('N', 'N'))
              + 5.0E-1 * greens_integral_tau(4, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) *
                pt->EWNS('N', 'N')->Pre()->Value() / dt;

    if (pt->EWNS('S', 'S'))
        yf += 5.0E-1 * greens_integral_tau(1, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) *
              f_ftn(pt->EWNS('S', 'S'))
              + 5.0E-1 * greens_integral_tau(1, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) *
                pt->EWNS('S', 'S')->Pre()->Value() / dt;

    if (pt->EWNS('N', 'E'))
        ent_y[3] = ent_y[1] * CalcHorizontalUCoefficient('E', xb, yp, xne, xnw, 0, epsilon)
                   - ent_y[1] * CalcHorizontalFCoefficient('E', xb, yp, xne, xnw, 0, epsilon) * 5.0E-1 / dt,
        ent_y[10] = ent_y[1] * CalcHorizontalPHICoefficient('E', xb, yp, xne, xnw, 0, epsilon) +
                    ent_y[8] * (xb - xnw) / (xne - xnw),
        yf += ent_y[1] * CalcHorizontalFCoefficient('E', xb, yp, xne, xnw, 0, epsilon) * 5.0E-1 *
              f_ftn(pt->EWNS('N', 'E'))
              + ent_y[1] * CalcHorizontalFCoefficient('E', xb, yp, xne, xnw, 0, epsilon) * 5.0E-1 / dt *
                pt->EWNS('N', 'E')->Pre()->Value()
              + 5.0E-1 * greens_integral_tau(4, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) * (xb - xnw) /
                (xne - xnw) * f_ftn(pt->EWNS('N', 'E'))
              + 5.0E-1 * greens_integral_tau(4, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) * (xb - xnw) /
                (xne - xnw) * pt->EWNS('N', 'E')->Pre()->Value() / dt;

    if (pt->EWNS('N', 'W'))
        ent_y[4] = ent_y[1] * CalcHorizontalUCoefficient('W', xb, yp, xne, xnw, 0, epsilon)
                   - ent_y[1] * CalcHorizontalFCoefficient('W', xb, yp, xne, xnw, 0, epsilon) * 5.0E-1 / dt,
        ent_y[11] = ent_y[1] * CalcHorizontalPHICoefficient('W', xb, yp, xne, xnw, 0, epsilon) +
                    ent_y[8] * (xne - xb) / (xne - xnw),
        yf += ent_y[1] * CalcHorizontalFCoefficient('W', xb, xp, xne, xnw, 0, epsilon) * 5.0E-1 *
              f_ftn(pt->EWNS('N', 'W'))
              + ent_y[1] * CalcHorizontalFCoefficient('W', xb, yp, xne, xnw, 0, epsilon) * 5.0E-1 / dt *
                pt->EWNS('N', 'W')->Pre()->Value()
              + 5.0E-1 * greens_integral_tau(4, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) * (xne - xb) /
                (xne - xnw) * f_ftn(pt->EWNS('N', 'W'))
              + 5.0E-1 * greens_integral_tau(4, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) * (xne - xb) /
                (xne - xnw) * pt->EWNS('N', 'W')->Pre()->Value() / dt;

    if (pt->EWNS('S', 'E'))
        ent_y[5] = ent_y[2] * CalcHorizontalUCoefficient('E', xb, ym, xse, xsw, 0, epsilon)
                   - ent_y[2] * CalcHorizontalFCoefficient('E', xb, ym, xse, xsw, 0, epsilon) * 5.0E-1 / dt,
        ent_y[12] = ent_y[2] * CalcHorizontalPHICoefficient('E', xb, ym, xse, xsw, 0, epsilon) +
                    ent_y[9] * (xb - xsw) / (xse - xsw),
        yf += ent_y[2] * CalcHorizontalFCoefficient('E', xb, ym, xse, xsw, 0, epsilon) * 5.0E-1 *
              f_ftn(pt->EWNS('S', 'E'))
              + ent_y[2] * CalcHorizontalFCoefficient('E', xb, ym, xse, xsw, 0, epsilon) * 5.0E-1 / dt *
                pt->EWNS('S', 'E')->Pre()->Value()
              + 5.0E-1 * greens_integral_tau(1, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) * (xb - xsw) /
                (xse - xsw) * f_ftn(pt->EWNS('S', 'E'))
              + 5.0E-1 * greens_integral_tau(1, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) * (xb - xsw) /
                (xse - xsw) * pt->EWNS('S', 'E')->Pre()->Value() / dt;

    if (pt->EWNS('S', 'W'))
        ent_y[5] = ent_y[2] * CalcHorizontalUCoefficient('W', xb, ym, xse, xsw, 0, epsilon)
                   - ent_y[2] * CalcHorizontalFCoefficient('W', xb, ym, xse, xsw, 0, epsilon) * 5.0E-1 / dt,
        ent_y[13] = ent_y[2] * CalcHorizontalPHICoefficient('W', xb, ym, xse, xsw, 0, epsilon) +
                    ent_y[9] * (xse - xb) / (xse - xsw),
        yf += ent_y[2] * CalcHorizontalFCoefficient('W', xb, ym, xse, xsw, 0, epsilon) * 5.0E-1 *
              f_ftn(pt->EWNS('S', 'W'))
              + ent_y[2] * CalcHorizontalFCoefficient('W', xb, ym, xse, xsw, 0, epsilon) * 5.0E-1 / dt *
                pt->EWNS('S', 'W')->Pre()->Value()
              + 5.0E-1 * greens_integral_tau(1, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) * (xse - xb) /
                (xse - xsw) * f_ftn(pt->EWNS('S', 'W'))
              + 5.0E-1 * greens_integral_tau(1, ym, yb, yp, xb, yb, 2, 0, epsilon, epsilon) * (xse - xb) /
                (xse - xsw) * pt->EWNS('S', 'W')->Pre()->Value() / dt;


    for (size_t i = 0; i < 14; i++) {
        ent_x[i] *= pt->Normal().Value('x');
        ent_y[i] *= pt->Normal().Value('y');
    }

    for (size_t i = 0; i < 14; i++) {
        if (func_x[i]) value += ent_x[i] * u_ftn_Dirichlet(func_x[i]);
        if (func_y[i]) value += ent_y[i] * u_ftn_Dirichlet(func_y[i]);
    }

    for (size_t i = 0; i < 14; i++) {
        if (func_x[i])
            printf("i = %02zu, ", i), printf("value = %23.16e, ", ent_x[i] * u_ftn_Dirichlet(func_x[i])), printf(
                    "ent_x = %23.16e, ", ent_x[i]), printf("u_ftn = %23.16e\n", u_ftn_Dirichlet(func_x[i]));
    }

    for (size_t i = 0; i < 14; i++) {
        if (func_y[i])
            printf("i = %02zu, ", i), printf("value = %23.16e, ", ent_y[i] * u_ftn_Dirichlet(func_y[i])), printf(
                    "ent_y = %23.16e, ", ent_y[i]), printf("u_ftn = %23.16e\n", u_ftn_Dirichlet(func_y[i]));
    }

    printf("xdat->F = %23.16e\n", xf * pt->Normal().Value('x'));
    printf("ydat->F = %23.16e\n", yf * pt->Normal().Value('y'));

    value += xf * pt->Normal().Value('x') + yf * pt->Normal().Value('y');

    return value;
}

double sol_y_approx(Point *pt) {

    double value = -1.0E0 * u_ftn(pt);
    double xb = pt->Coord().Value('x');
    double ym = pt->MinMaxCoordinate('y', 'm'), yb = pt->Coord().Value('y'), yp = pt->MinMaxCoordinate('y', 'p');
    double dt = pt->Dt();

    value += greens_coefficient_t(ym, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('S', 'S'));
    value += greens_coefficient_t(yp, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('N', 'N'));

    value -= greens_integral(1, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('S', 'S')->Phi());
    value -= greens_integral(2, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(pt->Phi());
    value -= greens_integral(3, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(pt->Phi());
    value -= greens_integral(4, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('N', 'N')->Phi());

    value += greens_integral(1, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * f_ftn(pt->EWNS('S', 'S')) * 5.0E-1;
    value += greens_integral(2, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * f_ftn(pt) * 5.0E-1;
    value += greens_integral(3, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * f_ftn(pt) * 5.0E-1;
    value += greens_integral(4, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * f_ftn(pt->EWNS('N', 'N')) * 5.0E-1;

    value -= greens_integral(1, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('S', 'S')) / dt * 5.0E-1;
    value -= greens_integral(2, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(pt) / dt * 5.0E-1;
    value -= greens_integral(3, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(pt) / dt * 5.0E-1;
    value -= greens_integral(4, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(pt->EWNS('N', 'N')) / dt * 5.0E-1;

    value += greens_integral(1, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * pt->EWNS('S', 'S')->Pre()->Value() / dt *
             5.0E-1;
    value += greens_integral(2, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * pt->Pre()->Value() / dt * 5.0E-1;
    value += greens_integral(3, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * pt->Pre()->Value() / dt * 5.0E-1;
    value += greens_integral(4, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * pt->EWNS('N', 'N')->Pre()->Value() / dt *
             5.0E-1;

    return value;
}

double sol_x_diff_diff_approx(Point *pt) {
    double value = ZeroValue;
    double xm = pt->MinMaxCoordinate('x', 'm'), xb = pt->Coord().Value('x'), xp = pt->MinMaxCoordinate('x', 'p');
    double yb = pt->Coord().Value('y');
    double dt = pt->Dt();
    double ent_x[14] = {ZeroValue,};
    double xf = -(5.0E-1 * HeadVelocity(pt)->F() + pt->Phi()->Value()) / HeadVelocity(pt)->MaterialProperty() +
                (5.0E-1 * (HeadVelocity(pt)->Value() - HeadVelocity(pt)->Pre()->Value()) / HeadVelocity(pt)->Dt()) /
                HeadVelocity(pt)->MaterialProperty();
    double yen, yes, ywn, yws, xne, xnw, xse, xsw;
    Point *func_x[14], *func_y[14];
    string azimuth_x[6] = {"EE", "WW", "EN", "ES", "WN", "WS"};
    string azimuth_y[6] = {"NN", "SS", "NE", "NW", "SE", "SW"};

    func_x[0] = pt->Diff('x');
    func_x[7] = pt->Phi();

    for (size_t i = 1; i < 7; i++) {
        if (pt->EWNS(azimuth_x[i - 1][0], azimuth_x[i - 1][1]))
            func_x[i] = pt->EWNS(azimuth_x[i - 1][0], azimuth_x[i - 1][1])->Diff('x');
        else func_x[i] = nullptr;
        if (pt->EWNS(azimuth_y[i - 1][0], azimuth_y[i - 1][1]))
            func_y[i] = pt->EWNS(azimuth_y[i - 1][0], azimuth_y[i - 1][1])->Diff('x');
        else func_y[i] = nullptr;
        if (pt->EWNS(azimuth_x[i - 1][0], azimuth_x[i - 1][1]))
            func_x[i + 7] = pt->EWNS(azimuth_x[i - 1][0], azimuth_x[i - 1][1])->Phi();
        else func_x[i + 7] = nullptr;
        if (pt->EWNS(azimuth_y[i - 1][0], azimuth_y[i - 1][1]))
            func_y[i + 7] = pt->EWNS(azimuth_y[i - 1][0], azimuth_y[i - 1][1])->Phi();
        else func_y[i + 7] = nullptr;
    }

    if (pt->EWNS('E', 'N')) yen = pt->EWNS('E', 'N')->Coord().Value('y');
    if (pt->EWNS('E', 'S')) yes = pt->EWNS('E', 'S')->Coord().Value('y');
    if (pt->EWNS('W', 'N')) ywn = pt->EWNS('W', 'N')->Coord().Value('y');
    if (pt->EWNS('W', 'S')) yws = pt->EWNS('W', 'S')->Coord().Value('y');
    if (pt->EWNS('N', 'E')) xne = pt->EWNS('N', 'E')->Coord().Value('x');
    if (pt->EWNS('N', 'W')) xnw = pt->EWNS('N', 'W')->Coord().Value('x');
    if (pt->EWNS('S', 'E')) xse = pt->EWNS('S', 'E')->Coord().Value('x');
    if (pt->EWNS('S', 'W')) xsw = pt->EWNS('S', 'W')->Coord().Value('x');

    ent_x[0] = 5.0E-1 * (greens_integral_ttau(2, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) +
                         greens_integral_ttau(3, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue)) / dt;
    ent_x[1] = greens_coefficient_ttau(xp, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) +
               5.0E-1 * greens_integral_ttau(4, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) / dt;
    ent_x[2] = greens_coefficient_ttau(xm, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) +
               5.0E-1 * greens_integral_ttau(1, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) / dt;

    ent_x[7] = -greens_integral_ttau(2, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue)
               -greens_integral_ttau(3, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue);
    ent_x[8] = -greens_integral_ttau(4, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue);
    ent_x[9] = -greens_integral_ttau(1, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue);

    xf = -5.0E-1 * (greens_integral_ttau(2, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) -
                    greens_integral_ttau(3, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue)) * f_ftn(pt)
         -5.0E-1 * (greens_integral_ttau(2, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) -
                    greens_integral_ttau(3, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue)) * pt->Pre()->Value() / dt;

    if (pt->EWNS('E', 'E'))
        xf -= 5.0E-1 * greens_integral_ttau(4, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) *
              f_ftn(pt->EWNS('E', 'E'))
              + 5.0E-1 * greens_integral_ttau(4, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) *
                pt->EWNS('E', 'E')->Pre()->Value() / dt;

    if (pt->EWNS('W', 'W'))
        xf -= 5.0E-1 * greens_integral_ttau(1, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) *
              f_ftn(pt->EWNS('W', 'W'))
              + 5.0E-1 * greens_integral_ttau(1, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) *
                pt->EWNS('W', 'W')->Pre()->Value() / dt;

    if (pt->EWNS('E', 'N'))
        ent_x[3] = ent_x[1] * CalcVerticalUCoefficient('N', xp, yb, yen, yes, 0, UnitValue)
                   - ent_x[1] * CalcVerticalFCoefficient('N', xp, yb, yen, yes, 0, UnitValue) * 5.0E-1 / dt,
        ent_x[10] = ent_x[1] * CalcVerticalPHICoefficient('N', xp, yb, yen, yes, 0, UnitValue) +
                    ent_x[8] * (yb - yes) / (yen - yes),
        xf += ent_x[1] * CalcVerticalFCoefficient('N', xp, yb, yen, yes, 0, UnitValue) * 5.0E-1 *
              f_ftn(pt->EWNS('E', 'N'))
              - 5.0E-1 * greens_integral_ttau(4, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) * (yb - yes) /
                (yen - yes) * f_ftn(pt->EWNS('E', 'N'))
              - 5.0E-1 * greens_integral_ttau(4, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) * (yb - yes) /
                (yen - yes) * pt->EWNS('E', 'N')->Pre()->Value() / dt
              + ent_x[1] * CalcVerticalFCoefficient('N', xp, yb, yen, yes, 0, UnitValue) * 5.0E-1 / dt *
                pt->EWNS('E', 'N')->Pre()->Value();
    if (pt->EWNS('E', 'S'))
        ent_x[4] = ent_x[1] * CalcVerticalUCoefficient('S', xp, yb, yen, yes, 0, UnitValue)
                   - ent_x[1] * CalcVerticalFCoefficient('S', xp, yb, yen, yes, 0, UnitValue) * 5.0E-1 / dt,
        ent_x[11] = ent_x[1] * CalcVerticalPHICoefficient('S', xp, yb, yen, yes, 0, UnitValue) +
                    ent_x[8] * (yen - yb) / (yen - yes),
        xf += ent_x[1] * CalcVerticalFCoefficient('S', xp, xb, yen, yes, 0, UnitValue) * 5.0E-1 *
              f_ftn(pt->EWNS('E', 'S'))
              - 5.0E-1 * greens_integral_ttau(4, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) * (yen - yb) /
                (yen - yes) * f_ftn(pt->EWNS('E', 'S'))
              - 5.0E-1 * greens_integral_ttau(4, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) * (yen - yb) /
                (yen - yes) * pt->EWNS('E', 'S')->Pre()->Value() / dt
              + ent_x[1] * CalcVerticalFCoefficient('S', xp, yb, yen, yes, 0, UnitValue) * 5.0E-1 / dt *
                pt->EWNS('E', 'S')->Pre()->Value();
    if (pt->EWNS('W', 'N'))
        ent_x[5] = ent_x[2] * CalcVerticalUCoefficient('N', xm, yb, ywn, yws, 0, UnitValue)
                   - ent_x[2] * CalcVerticalFCoefficient('N', xm, yb, ywn, yws, 0, UnitValue) * 5.0E-1 / dt,
        ent_x[12] = ent_x[2] * CalcVerticalPHICoefficient('N', xm, yb, ywn, yws, 0, UnitValue) +
                    ent_x[9] * (yb - yws) / (ywn - yws),
        xf += ent_x[2] * CalcVerticalFCoefficient('N', xm, yb, ywn, yws, 0, UnitValue) * 5.0E-1 *
              f_ftn(pt->EWNS('W', 'N'))
              - 5.0E-1 * greens_integral_ttau(1, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) * (yb - yws) /
                (ywn - yws) * f_ftn(pt->EWNS('W', 'N'))
              - 5.0E-1 * greens_integral_ttau(1, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) * (yb - yws) /
                (ywn - yws) * pt->EWNS('W', 'N')->Pre()->Value() / dt
              + ent_x[2] * CalcVerticalFCoefficient('N', xm, yb, ywn, yws, 0, UnitValue) * 5.0E-1 / dt *
                pt->EWNS('W', 'N')->Pre()->Value();
    if (pt->EWNS('W', 'S'))
        ent_x[5] = ent_x[2] * CalcVerticalUCoefficient('S', xm, yb, ywn, yws, 0, UnitValue)
                   - ent_x[2] * CalcVerticalFCoefficient('S', xm, yb, ywn, yws, 0, UnitValue) * 5.0E-1 / dt,
        ent_x[13] = ent_x[2] * CalcVerticalPHICoefficient('S', xm, yb, ywn, yws, 0, UnitValue) +
                    ent_x[9] * (ywn - yb) / (ywn - yws),
        xf += ent_x[2] * CalcVerticalFCoefficient('S', xm, yb, ywn, yws, 0, UnitValue) * 5.0E-1 *
              f_ftn(pt->EWNS('W', 'S'))
              - 5.0E-1 * greens_integral_ttau(1, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) * (ywn - yb) /
                (ywn - yws) * f_ftn(pt->EWNS('W', 'S'))
              - 5.0E-1 * greens_integral_ttau(1, xm, xb, xp, xb, yb, 1, 0, UnitValue, UnitValue) * (ywn - yb) /
                (ywn - yws) * pt->EWNS('W', 'S')->Pre()->Value() / dt
              + ent_x[2] * CalcVerticalFCoefficient('S', xm, yb, ywn, yws, 0, UnitValue) * 5.0E-1 / dt *
                pt->EWNS('W', 'S')->Pre()->Value();


    for (double & i : ent_x) {
        i *= pt->Normal().Value('x');
    }

    for (size_t i = 0; i < 14; i++) {
        if (func_x[i]) value += ent_x[i] * u_ftn_Dirichlet(func_x[i]);
    }

    for (size_t i = 0; i < 14; i++) {
        if (func_x[i])
            printf("i = %02zu, ", i), printf("value = %23.16e, ", ent_x[i] * u_ftn_Dirichlet(func_x[i])), printf(
                    "ent_x = %23.16e, ", ent_x[i]), printf("u_ftn = %23.16e\n", u_ftn_Dirichlet(func_x[i]));
    }

    printf("xdat->F = %23.16e\n", xf * pt->Normal().Value('x'));

    value += xf;

    return value;
}

double
sol_y_approx(const double xm, const double xb, const double xp, const double ym, const double yb, const double yp,
             const double T, const double Dt) {

    double value = ZeroValue;
    double dt = Dt;
    double VALUE[5] = {ZeroValue,};

    for (size_t i = 0; i < 5; i++) {
        value = ZeroValue;

        dt = Dt * pow(5.0E-1, i);

        value += greens_coefficient_t(yp, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(xb, yp, T);
        value += greens_coefficient_t(ym, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(xb, ym, T);

        value -= greens_integral(1, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * phi_ftn(xb, ym, T);
        value -= greens_integral(2, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * phi_ftn(xb, yb, T);
        value -= greens_integral(3, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * phi_ftn(xb, yb, T);
        value -= greens_integral(4, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * phi_ftn(xb, yp, T);

        value += greens_integral(1, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * f_ftn(xb, ym, T) * 5.0E-1;
        value += greens_integral(2, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * f_ftn(xb, yb, T) * 5.0E-1;
        value += greens_integral(3, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * f_ftn(xb, yb, T) * 5.0E-1;
        value += greens_integral(4, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * f_ftn(xb, yp, T) * 5.0E-1;

        value -= greens_integral(1, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(xb, ym, T) / dt * 5.0E-1;
        value -= greens_integral(2, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(xb, yb, T) / dt * 5.0E-1;
        value -= greens_integral(3, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(xb, yb, T) / dt * 5.0E-1;
        value -= greens_integral(4, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(xb, yp, T) / dt * 5.0E-1;

        value += greens_integral(1, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(xb, ym, T - dt) / dt * 5.0E-1;
        value += greens_integral(2, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(xb, yb, T - dt) / dt * 5.0E-1;
        value += greens_integral(3, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(xb, yb, T - dt) / dt * 5.0E-1;
        value += greens_integral(4, ym, yb, yp, xb, yb, 2, 0, 1.0E0, 1.0E0) * u_ftn(xb, yp, T - dt) / dt * 5.0E-1;

        printf("value = %23.16e; ", fabs(u_ftn(xb, yb, T) - value));

        VALUE[i] = fabs(u_ftn(xb, yb, T) - value);

    }
    printf("\n");

    double MinOrder = VALUE[0] / VALUE[1];
    double Order = VALUE[0] / VALUE[1];

    for (size_t i = 1; i < 4; i++) {
        Order = VALUE[i] / VALUE[i + 1];
        if (MinOrder > Order) MinOrder = Order;
    }

    return MinOrder;

}

double
sol_y_approx(const double ym, const double yb, const double yp, const double xb, const double T, const double Dt) {
    const double value = (-0.25 * Dt * pow(yb, 2) * ym + 0.25 * Dt * pow(yb, 2) * yp + 0.25 * Dt * yb * pow(ym, 2) -
                          0.25 * Dt * yb * pow(yp, 2) - 0.25 * Dt * pow(ym, 2) * yp + 0.25 * Dt * ym * pow(yp, 2)
                          + 1.0 * pow(T, 2) * ym - 1.0 * pow(T, 2) * yp +
                          0.16666666666666666 * pow(xb, 2) * pow(yb, 3) * sin(xb * ym) -
                          0.16666666666666666 * pow(xb, 2) * pow(yb, 3) * sin(xb * yp)
                          - 0.33333333333333331 * pow(xb, 2) * pow(yb, 2) * ym * sin(xb * yb) -
                          0.33333333333333331 * pow(xb, 2) * pow(yb, 2) * ym * sin(xb * ym) +
                          0.16666666666666666 * pow(xb, 2) * pow(yb, 2) * ym * sin(xb * yp)
                          + 0.33333333333333331 * pow(xb, 2) * pow(yb, 2) * yp * sin(xb * yb) -
                          0.16666666666666666 * pow(xb, 2) * pow(yb, 2) * yp * sin(xb * ym) +
                          0.33333333333333331 * pow(xb, 2) * pow(yb, 2) * yp * sin(xb * yp)
                          + 0.33333333333333331 * pow(xb, 2) * yb * pow(ym, 2) * sin(xb * yb) +
                          0.16666666666666666 * pow(xb, 2) * yb * pow(ym, 2) * sin(xb * ym) +
                          0.33333333333333331 * pow(xb, 2) * yb * ym * yp * sin(xb * ym)
                          - 0.33333333333333331 * pow(xb, 2) * yb * ym * yp * sin(xb * yp) -
                          0.33333333333333331 * pow(xb, 2) * yb * pow(yp, 2) * sin(xb * yb) -
                          0.16666666666666666 * pow(xb, 2) * yb * pow(yp, 2) * sin(xb * yp)
                          - 0.33333333333333331 * pow(xb, 2) * pow(ym, 2) * yp * sin(xb * yb) -
                          0.16666666666666666 * pow(xb, 2) * pow(ym, 2) * yp * sin(xb * ym) +
                          0.33333333333333331 * pow(xb, 2) * ym * pow(yp, 2) * sin(xb * yb)
                          + 0.16666666666666666 * pow(xb, 2) * ym * pow(yp, 2) * sin(xb * yp) +
                          1.0 * yb * sin(xb * ym) - 1.0 * yb * sin(xb * yp) + 1.0 * ym * sin(xb * yp) -
                          1.0 * yp * sin(xb * ym)) / (ym - yp);

    return value;
}

#endif
