#include "ground_state.cpp"
#include "constants.h"

const double omega_L = (2*M_PI*c/(lambda_L/t_au));
const double E_0 = sqrt(I/I_a);// 0;
const double tau = (tau_p/t_au)/sqrt(2*log(2.));
const double r_osc = std::abs(E_0/(omega_L*omega_L));

double dA_dt(double t)
{
    return -(c*E_0/omega_L)*exp(-t*t/(tau*tau))*(omega_L*cos(omega_L*t) -
                  sin(omega_L*t)*(2*t/(tau*tau)));
}

double E_t(double t)
{
    return (-1/c)*dA_dt(t); // домножил на 0 23.03.2023, чтобы проверить дипольное ускорение в отсутствие внешнего электрического поля
}
//аналитически вычисленный градиент модельного потенциала
double grad_V(double x)
{
    return x/sqrt(pow(x*x+2, 3));
}

int main()
{
    const int M = 32500;
    const double dz = 0.1; // T = M*dz (a.u)
    const double Z_au_min = -4*r_osc;// -221.721
    const double Z_au_max = 4*r_osc;// 221.721
    const double N1 = (Z_au_max - Z_au_min)/dz;
    const int N = pow(2, (int(log(N1)/log(2)))+1);
    std::cout << "N = " << N << std::endl;
    //std::cout << "E_0 = " << E_0 << std::endl;
    //std::cout << "r_osc = " << r_osc << std::endl;
    //std::cout << "omega_L = " << omega_L << std::endl;
    //std::cout << "tau = " << tau << std::endl;
    //return 0;
    const double dt = 0.02;// в атомных единицах
    const double dp = 2.0*M_PI/(dz*N);
    const double dw = 2.0*M_PI/(dt*M);
    std::vector<double> coordinate(N);
    double *func = new double[N];
    fftw_complex *func_in = new fftw_complex[N];
    fftw_complex *func_out = new fftw_complex[N];
    fftw_plan plan_fwd, plan_bwd;
    plan_fwd = fftw_plan_dft_1d(N, func_in, func_out, FFTW_FORWARD, FFTW_MEASURE);
    plan_bwd = fftw_plan_dft_1d(N, func_out, func_in, FFTW_BACKWARD, FFTW_MEASURE);

    double Z = dz/2;
    for(int i = N/2; i < N; ++i)
    {
        coordinate[i] = Z;
        Z += dz;//dx;
    }
    Z = - dz/2;
    for(int i = N/2 - 1; i >= 0; --i)
    {
        coordinate[i] = Z;
        Z -= dz;//dx;
    }

    /*если не надо пересчитывать основное состояние, например, при изменении параметров,
    то ради ускорения вычислений будем считывать данные с ранее записанного файла*/

    //ground_state(func, N, dz, coordinate);
    //return 0; //new
//#if 0
    std::ifstream file_with_gs("ground_state.dat");
    int var = 0;
    if(file_with_gs)
    {
        double r = 0.;
        while(true)
        {
            file_with_gs >> r >> func[var];
            if(file_with_gs.eof()) break;
            var++;
        }
    }

    for(int i = 0; i < N; ++i)
    {
        func_in[i][0] = func[i];
        func_in[i][1] = 0.0;
    }
//#endif
    delete [] func;
    func = nullptr;

#if 0
    double X = Xmin;
    for(int i = 0; i < N; ++i)
    {
        coordinate[i] = X;
        X += dx;
    }
#endif
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
       t[i] = -3*tau + i*dt;
    }

    fftw_complex *a_t = new fftw_complex[M];

    for(int i = 0; i < M; ++i)
    {
        double Integral_sqrpsi_gradV = 0;
        for(int j = 0; j < N; ++j)
        {
            double re_part_psi_in = func_in[j][0];
            double im_part_psi_in = func_in[j][1];
            func_in[j][0] = re_part_psi_in * cos((V(coordinate[j])-coordinate[j]*E_t(t[i]))*dt) +
                    im_part_psi_in * sin((V(coordinate[j])-coordinate[j]*E_t(t[i]))*dt);
            func_in[j][1] = -re_part_psi_in * sin((V(coordinate[j])-coordinate[j]*E_t(t[i]))*dt) +
                    im_part_psi_in * cos((V(coordinate[j])-coordinate[j]*E_t(t[i]))*dt);
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
            Integral_sqrpsi_gradV += (func_in[j][0]*func_in[j][0]+
                    func_in[j][1]*func_in[j][1])*grad_V(coordinate[j])*dt;
        }

        a_t[i][0] = E_t(t[i]) - Integral_sqrpsi_gradV;// на каждом шаге по времени находим дипольное ускорение
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

//#if 0
    std::ofstream _out_("a_t.dat");
    for(int i = 0; i < M; ++i)
    {
        _out_ << t[i] << "\t" << a_t[i][0] << std::endl;
    }
    _out_.close();
    //return 0;   //временный выход из программы 22.03.2023
//#endif
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
    double omega_Nyq = 2*M_PI/(2*dt);
    for(int i = 0; i < M; ++i)
    {
        omega[i] = -omega_Nyq + dw*i;
    }
#if 0
    for(int i = M/2+1; i < M; ++i)
    {
        omega[i] = -dw*(M-i);//-omega_Nyq + dw*(i - M/2);
    }
#endif
    std::ofstream print_a_omega("a_omega.dat");
    for(int i = M/2; i < M; ++i)
    {
        print_a_omega << omega[i-M/2] <<'\t' << std::sqrt(a_omega[i][0]*a_omega[i][0]+
                a_omega[i][1]*a_omega[i][1]) << std::endl;
    }
    for(int i = 0; i < M/2; ++i)
    {
        print_a_omega << omega[i+M/2] <<'\t' << std::sqrt(a_omega[i][0]*a_omega[i][0]+
                a_omega[i][1]*a_omega[i][1]) << std::endl;
    }
    print_a_omega.close();
    //return 0; //new 3.12.2022
    double *HHG = new double[M];
    for(int i = 0; i < M; ++i)
    {
        HHG[i] = a_omega[i][0]*a_omega[i][0]+a_omega[i][1]*a_omega[i][1];//по определению спектра
    }

    std::ofstream _out("HHG_spectrum.dat");
    //_out.precision(10);
    for(int i = M/2; i < M; ++i)
    {
        _out << omega[i-M/2]/omega_L << "\t" << HHG[i] << std::endl;
    }
    for(int i = 0; i < M/2; ++i)
    {
        //_out << std::fixed << HHG[i] << "\t" << omega[i]/omega_L << std::endl;
        _out << omega[i+M/2]/omega_L << "\t" << HHG[i] << std::endl;
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
