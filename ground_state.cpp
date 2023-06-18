#include "ground_state.h"

double V(double x)
{
    return (-1.)/(sqrt(pow(x,2)+2));
}

double Psi_solution(double x)
{
    return exp(-sqrt(x*x+2))*(1+sqrt(x*x+2));
}

void ground_state(double func[], const int N,
                  const double dx, std::vector<double> coord)
{
    const int M = 20000;
    const double dt = 0.02;
    const double dp = 2.0*M_PI/(dx*N);
    fftw_complex *func_in = new fftw_complex[N];
    fftw_complex *func_out = new fftw_complex[N];
    fftw_plan plan_fwd, plan_bwd;
    plan_fwd = fftw_plan_dft_1d(N, func_in, func_out, FFTW_FORWARD, FFTW_MEASURE);
    plan_bwd = fftw_plan_dft_1d(N, func_out, func_in, FFTW_BACKWARD, FFTW_MEASURE);

    for (int i = 0; i < N; ++i)
    {
        func_in[i][0] = exp(-coord[i]*coord[i]);//exp(-x^2);
        func_in[i][1] = 0.;
    }

//#if 0
    std::ofstream print_an("ground_state_analitical.dat");
    for(int i = 0; i < N; ++i)
    {
        //print_an << coordinate[i] << '\t' << Psi_solution(coordinate[i]) << std::endl;
        print_an << coord[i] << '\t' << Psi_solution(coord[i]) << std::endl;
    }
    print_an.close();
//#endif

    std::ofstream print_func_in("real_part_func_in.dat");
    for(int i = 0; i < N; ++i)
    {
        //print_func_in << coordinate[i] << '\t' << func_in[i][0] << std::endl;
        print_func_in << coord[i] << '\t' << func_in[i][0] << std::endl;
    }
    print_func_in.close();

    double *p = new double[N];
    for(int i = 0; i < N/2; ++i)
    {
        p[i] = dp*i;
    }
    for(int i = N/2; i < N; ++i)
    {
        p[i] = -dp*(N-i);
    }

    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            //func_in[j][0] = exp(-V(coordinate[j])*dt)*func_in[j][0];
            func_in[j][0] = exp(-V(coord[j])*dt)*func_in[j][0];
        }

        fftw_execute(plan_fwd);

        for (int j = 0; j < N; ++j)
        {
            func_out[j][0] = exp(-p[j]*p[j]*dt/2.)*func_out[j][0];
            func_out[j][1] = exp(-p[j]*p[j]*dt/2.)*func_out[j][1];
        }
            fftw_execute(plan_bwd);

        for (int j = 0; j < N; ++j)
        {
            func_in[j][0] = (1./N)*func_in[j][0];
        }

        double psi_norma = 0;
        for (int j = 0; j < N; ++j)
        {
             psi_norma = psi_norma + (func_in[j][0]*func_in[j][0])*dx;
        }


        double re_member = 1./sqrt(psi_norma);

        for (int j = 0; j < N; ++j)
        {
             func_in[j][0] = func_in[j][0]*re_member;
        }

        for (int j = 0; j < N; ++j)
        {
             func_in[j][1] = 0;
        }
    }

    for(int i = 0; i < N; ++i)
    {
        func[i] = func_in[i][0];
    }
//#if 0
    std::ofstream print_gs("ground_state.dat");
    for(int i = 0; i < N; ++i)
    {
        //print_gs << coordinate[i] << '\t' << func_in[i][0] << std::endl;
        print_gs << coord[i] << '\t' << func_in[i][0] << std::endl;
    }
    print_gs.close();
//#endif

    fftw_destroy_plan(plan_bwd);
    fftw_destroy_plan(plan_fwd);

    delete [] p;
    p = nullptr;
    delete [] func_out;
    func_out = nullptr;
    delete [] func_in;
    func_in = nullptr;
}
