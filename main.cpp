#include <iostream>
#include <fstream>
#include "ground_state.cpp"
#include "constants.h"

const double omega_L = (2*M_PI*c/(lambda_L/t_au));
const double E_0 = sqrt(I/I_a);
const double tau = (tau_p/t_au)/sqrt(2*log(2.));
const double r_osc = std::abs(E_0/(omega_L*omega_L));

double dA_dt(double t)
{
    return -(c*E_0/omega_L)*exp(-t*t/(tau*tau))*(omega_L*cos(omega_L*t) -
                  sin(omega_L*t)*(2*t/(tau*tau)));
}

//аналитически вычисленный градиент модельного потенциала
double grad_V(double x)
{
    return x/sqrt(pow(x*x+2, 3));
}

int main()
{
    const int N = 4096;
    const int M = 10000;
    //const double Xmin = -15.;
    //const double Xmax = 15.;
    const double Z_au_min = -4*r_osc;
    const double Z_au_max = 4*r_osc;
    std::cout << "E_0 = " << E_0 << std::endl;
    std::cout << "r_osc = " << r_osc << std::endl;
    std::cout << "omega_L = " << omega_L << std::endl;
    std::cout << "Zmax = " << Z_au_max<<std::endl;
    const double dz = (Z_au_max - Z_au_min)/N;
    std::cout << "dz = " << dz << std::endl;
    //const double t_min = (-4)*(tau_p/t_au);//это для расчёта интеграла остаточной плотности тока
    //const double t_max = 4*(tau_p/t_au);
    const double dt = 0.02;// в атомных единицах
    //const int time_steps = (t_max - t_min)/dt;
    //const double dx = (Xmax - Xmin)/N;
    const double dp = 2.0*M_PI/(dz*N);
    const double dw = 2.0*M_PI/(dt*N);
    std::vector<double> coordinate(N);
    double *func = new double[N];
    fftw_complex *func_in = new fftw_complex[N];
    fftw_complex *func_out = new fftw_complex[N];
    fftw_plan plan_fwd, plan_bwd;
    plan_fwd = fftw_plan_dft_1d(N, func_in, func_out, FFTW_FORWARD, FFTW_MEASURE);
    plan_bwd = fftw_plan_dft_1d(N, func_out, func_in, FFTW_BACKWARD, FFTW_MEASURE);

    //ground_state(func, N, Xmin, Xmax);
    ground_state(func, N, Z_au_min, Z_au_max);
    //return 0; //new
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

    double X = Z_au_min;
    for(int i = 0; i < N; ++i)
    {
        coordinate[i] = X;
        X += dz;
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

    double Integral_sqrpsi_gradV = 0;

    for(int i = 0; i < M; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            double re_part_psi_in = func_in[j][0];
            double im_part_psi_in = func_in[j][1];
            func_in[j][0] = re_part_psi_in * cos(V(coordinate[j])*dt) +
                    im_part_psi_in * sin(V(coordinate[j])*dt);
            func_in[j][1] = -re_part_psi_in * sin(V(coordinate[j])*dt) +
                    im_part_psi_in * cos(V(coordinate[j])*dt);
        }

        fftw_execute(plan_fwd);

        for(int j = 0; j < N; ++j)
        {
            double re_part_psi_out = func_out[j][0];
            double im_part_psi_out = func_out[j][1];
            func_out[j][0] = re_part_psi_out * cos(p[j]*p[j]*dt/2.) +
                    im_part_psi_out * sin(p[j]*p[j]*dt/2.);
            func_out[j][1] = -re_part_psi_out * sin(p[j]*p[j]*dt/2.) +
                    im_part_psi_out * cos(p[j]*p[j]*dt/2.);
        }

        fftw_execute(plan_bwd);

        for(int j = 0; j < N; ++j)
        {
            func_in[j][0] = (1./N)*func_in[j][0];
            func_in[j][1] = (1./N)*func_in[j][1];
        }

        for(int j = 0; j < N; ++j)
        {
             Integral_sqrpsi_gradV += (func_in[j][0]*func_in[j][0] +
                func_in[j][1]*func_in[j][1])*grad_V(coordinate[j])*dz;
        }

        a_t[i][0] = -(1/c)*dA_dt(t[i]) - Integral_sqrpsi_gradV;// на каждом шаге по времени находим дипольное ускорение
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

#if 0
    std::ofstream _out_("a_t.dat");
    for(int i = 0; i < M; ++i)
    {
        _out_ << t[i] << "\t" << a_t[i][0] << std::endl;
    }
    _out_.close();
#endif

    fftw_complex *a_in = new fftw_complex[M];
    fftw_complex *a_omega = new fftw_complex[M];
    fftw_plan plan_fwd_a = fftw_plan_dft_1d(M, a_in, a_omega, FFTW_FORWARD, FFTW_MEASURE);

    for(int i = 0; i < M; ++i)
    {
        a_in[i][0] = a_t[i][0];
        a_in[i][1] = a_t[i][1];
    }

    fftw_execute(plan_fwd_a);

#if 0
    std::ofstream print_omega("a_omega.dat");
    for(int i = 0; i < M; ++i)
    {
        print_omega << a_omega[i][0] << '\t' << a_omega[i][1] << std::endl;
    }
    print_omega.close();
#endif

    double *omega = new double[M];
    for(int i = 0; i < M/2; ++i)
    {
        omega[i] = dw*i;
    }
//#if 0
    for(int i = M/2; i < M; ++i)
    {
        omega[i] = -dw*(M - i);
    }
//#endif
    double *HHG = new double[M];
    for(int i = 0; i < M; ++i)
    {
        HHG[i] = a_omega[i][0]*a_omega[i][0]+a_omega[i][1]*a_omega[i][1];//по определению спектра
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
