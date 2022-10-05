#include <iostream>
#include <fstream>
#include "ground_state.cpp"
#include "constants.h"

const double omega_L = (2*M_PI*c/lambda_L);
const double E_0 = sqrt(I_a/I);
const double tau = tau_p/sqrt(2*log(2.));

double dA_dt(double t)
{
    return -(c*E_0/omega_L)*(omega_L*cos(omega_L*t)*exp(-t*t/(tau*tau)) -
                  sin(omega_L*t)*(2*t/(tau*tau))*exp(-t*t/(tau*tau)));
}

//аналитически вычисленный градиент модельного потенциала
double grad_V(double x)
{
    return x/sqrt(pow(x*x+2, 3));
}

int main()
{
    const int N = 2048;
    const int M = 10000;
    const double Xmin = -15.;
    const double Xmax = 15.;
    const double dt = 0.01;
    const double dx = (Xmax - Xmin)/N;
    const double dp = 2.0*M_PI/(dx*N);
    const double dw = 2.0*M_PI/(dt*N);
    std::vector<double> coordinate(N);
    double *func = new double[N];
    fftw_complex *func_in = new fftw_complex[N];
    fftw_complex *func_out = new fftw_complex[N];
    fftw_plan plan_fwd, plan_bwd;
    plan_fwd = fftw_plan_dft_1d(N, func_in, func_out, FFTW_FORWARD, FFTW_MEASURE);
    plan_bwd = fftw_plan_dft_1d(N, func_out, func_in, FFTW_BACKWARD, FFTW_MEASURE);

    ground_state(func, N, Xmin, Xmax);
#if 0
    std::ofstream print_func("func.dat");
    double x = Xmin;
    for(int i = 0; i < N; ++i)
    {
        print_func << x << "\t" << func[i] << std::endl;
        x += dx;
    }
    print_func.close();
#endif
    for(int i = 0; i < N; ++i)
    {
        func_in[i][0] = func[i];
        func_in[i][1] = 0.0;
    }

    delete [] func;
    func = nullptr;

    double X = Xmin;
    for(int i = 0; i < N; ++i)
    {
        coordinate[i] = X;
        X += dx;
    }

    double *p = new double[N];
    for(int i = 0; i < N/2; ++i)
    {
        p[i] = dp*i;
    }
    for(int i = N/2; i < N; ++i)
    {
        p[i] = -dp*(N-i);
    }

    double *t = new double[M];
    for(int i = 0; i < M; ++i)
    {
       t[i] = i*dt;
    }

    fftw_complex *a_t = new fftw_complex[M];
    double re_integral;
    double im_integral;
    double re_norma;
    double im_norma;
    double Integral_sqrpsi_gradV = 0;

    for(int i = 0; i < M; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            func_in[j][0] = func_in[j][0] * cos(V(coordinate[j])*dt) +
                    func_in[j][1] * sin(V(coordinate[j])*dt);
            func_in[j][1] = -func_in[j][0] * sin(V(coordinate[j])*dt) +
                    func_in[j][1] * cos(V(coordinate[j])*dt);
        }

        fftw_execute(plan_fwd);

        for(int j = 0; j < N; ++j)
        {
            func_out[j][0] = func_out[j][0] * cos(p[j]*p[j]*dt/2.) +
                    func_out[j][1] * sin(p[j]*p[j]*dt/2.);
            func_out[j][1] = -func_out[j][0] * sin(p[j]*p[j]*dt/2.) +
                    func_out[j][1] * cos(p[j]*p[j]*dt/2.);
        }

        fftw_execute(plan_bwd);

        for(int j = 0; j < N; ++j)
        {
            func_in[j][0] = (1./N)*func_in[j][0];
            func_in[j][1] = (1./N)*func_in[j][1];
        }

        re_integral = 0;
        im_integral = 0;

        for(int j = 0; j < N; ++j)
        {
            re_integral = re_integral + (func_in[j][0]*func_in[j][0])*dx;
            im_integral = im_integral + (func_in[j][1]*func_in[j][1])*dx;
        }

        re_norma = 1./sqrt(re_integral);
        im_norma = 1./sqrt(im_integral);

        for(int j = 0; j < N; ++j)
        {
            func_in[j][0] = func_in[j][0]*re_norma;
            func_in[j][1] = func_in[j][1]*im_norma;
        }

        for(int j = 0; j < N; ++j)
        {
             Integral_sqrpsi_gradV += (func_in[j][0]*func_in[j][0] +
                func_in[j][1]*func_in[j][1])*grad_V(coordinate[j])*dx;
        }

        a_t[i][0] = -(1/c)*dA_dt(t[i]) - Integral_sqrpsi_gradV;
        a_t[i][1] = 0.;

    }
#if 0
    std::ofstream print_psi("psi.dat");
    for(int i = 0; i < N; ++i)
    {
        print_psi << func_in[i][0] << '\t' << func_in[i][1] << std::endl;
    }
    print_psi.close();
#endif

    delete [] func_out;
    func_out = nullptr;
    delete [] func_in;
    func_in = nullptr;

    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_bwd);

    std::ofstream _out_("a_t.dat");
    for(int i = 0; i < M; ++i)
    {
        _out_ << t[i] << "\t" << a_t[i][0] << std::endl;
    }
    _out_.close();

    fftw_complex *a_in = new fftw_complex[M];
    fftw_complex *a_omega = new fftw_complex[M];
    fftw_plan plan_fwd_a = fftw_plan_dft_1d(M, a_in, a_omega, FFTW_FORWARD, FFTW_MEASURE);

    for(int i = 0; i < M; ++i)
    {
        a_in[i][0] = a_t[i][0];
        a_in[i][1] = a_t[i][1];
    }

    fftw_execute(plan_fwd_a);

    std::ofstream print_omega("a_omega.dat");
    for(int i = 0; i < M; ++i)
    {
        print_omega << a_omega[i][0] << '\t' << a_omega[i][1] << std::endl;
    }
    print_omega.close();

    double *omega = new double[M];
    for(int i = 0; i < M; ++i)
    {
        omega[i] = dw*i;
    }
#if 0
    for(int i = M/2; i < M; ++i)
    {
        omega[i] = -dw*(M - i);
    }
#endif
    double *HHG = new double[M];
    for(int i = 0; i < M; ++i)
    {
        HHG[i] = a_omega[i][0]*a_omega[i][0];//+a_omega[i][1]*a_omega[i][1];
    }

    std::ofstream _out("HHG_spectrum.dat");
    //_out.precision(10);
    for(int i = 0; i < M; ++i)
    {
        //_out << std::fixed << HHG[i] << "\t" << omega[i]/omega_L << std::endl;
        _out << omega[i]/omega_L << "\t" << HHG[i] << std::endl;
    }
    _out.close();

    fftw_destroy_plan(plan_fwd_a);

    delete [] HHG;
    HHG = nullptr;
    delete [] omega;
    omega = nullptr;
    delete [] a_in;
    a_in = nullptr;
    delete [] a_omega;
    a_omega = nullptr;
    delete [] a_t;
    a_t = nullptr;
    delete [] t;
    t = nullptr;
    delete [] p;
    p = nullptr;

    return 0;
}
