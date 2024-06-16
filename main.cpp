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

//функция-маска для устранения шумов, связанных с остаточными осцилляциями дипольного ускорения
double Mask(double t, double t_max)
{
    return 0.5*(1.0 + tanh((-(t - (t_max - t_tau))/(t_tau/8))*(180/M_PI)));
}

//волновая функция основного состояния гармонического осциллятора, зависимость от х
double Psi_ground_state_osc(double x)
{
    return pow((omega_L/M_PI), 0.25) * exp((-omega_L*pow(x, 2))/2);
}

int main()
{
    const int M = 32500;    // было 32500
    const double dz = 0.1; // T = M*dz (a.u)
    const double Z_au_min = -221.721/2;// -221.721    -2*4*r_osc
    const double Z_au_max = 221.721/2;// 221.721     2*4*r_osc
    const double N1 = (Z_au_max - Z_au_min)/dz;
    const int N = pow(2, (int(log(N1)/log(2)))+1);
    std::cout << "N = " << N << std::endl;
    //return 0;
    const double dt = 0.02;// в атомных единицах
    const double dp = 2.0*M_PI/(dz*N);
    const double dw = 2.0*M_PI/(dt*M);
    std::vector<double> coordinate(N);      //массив координаты z
    std::vector<double> array_E_t(M);       //массив E(t)
    //std::vector<double> dipole_moment(M);   //массив диполного момента
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

    double Z_abs = coordinate[N-1] - 50;    //coordinate[N-1] - 50
    double CAP = 1.0;//exp(-0.1*dt);

#if 0
    std::ofstream print_grad_V("grad_V.dat");
    for(int i = 0; i < N; ++i)
    {
        print_grad_V << grad_V(coordinate[i]) << '\t' << coordinate[i] << std::endl;
    }
    print_grad_V.close();
    //return 0;
#endif

    /*если не надо пересчитывать основное состояние, например, при изменении параметров,
    то ради ускорения вычислений будем считывать данные с ранее записанного файла*/

    ground_state(func, N, dz, coordinate);
    //return 0; //new
#if 0
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
#endif
//#if 0
    for(int i = 0; i < N; ++i)
    {
        //внесём небольшое смещение в.ф. относительно центра
        /*
        if(i < N - 50) {
            func_in[i][0] = func[i+50];
        }
        else {
           func_in[i][0] = func[N-1];
        }*/

        func_in[i][0] = func[i];
        func_in[i][1] = 0.0;
    }
//#endif

    delete [] func;
    func = nullptr;

    //сдвинутое по координате начальное значение волновой функции (27.08.2023)
    std::ofstream print_ground_state("ground_state.dat");
    for(int i = 0; i < N; ++i)  {
        double x = coordinate[i];
        print_ground_state << x << '\t' << sqrt(pow(func_in[i][0],2) +
                           pow(func_in[i][1],2)) << std::endl;
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
       t[i] = -3*tau + i*dt;
    }

    double t_max = t[M-1];
    double *a_t = new double[M];

#if 0
    //записываем точки для построения графика функции маски
    std::ofstream print_acceleration_mult_mask("Mask_a_t.dat");
    for(int i = 0; i < M; ++i)
    {
        //print_acceleration_mult_mask << t[i] << '\t' << a_in[i][0] << std::endl;
        print_acceleration_mult_mask << t[i] << '\t' << Mask(t[i],t_max) << std::endl;
    }
    print_acceleration_mult_mask.close();
    return 0; // 26.03.2023
    //маска расчитана правильно
#endif

    //double Integral_sqrpsi_gradV;
    //double Integral_sqrpsi_z;
    //#if 0
    //рассчитаем данные для зависимости электрического поля от времени
    for(int i = 0; i < M; ++i)
    {
       array_E_t[i] = E_t(i);
    }

    std::ofstream print_E_t("E_t.dat");
    for(int i = 0; i < M; ++i)
    {
        print_E_t << t[i] << '\t' << array_E_t[i] << std::endl;
    }
    print_E_t.close();
    //#endif
    //return 0; // 20.07.2023 построение зависимости электрического поля от времени

    for(int i = 0; i < M; ++i)
    {
        double Integral_sqrpsi_gradV = 0;
        //Integral_sqrpsi_z = 0;           //для расчёта дипольного момента
        for(int j = 0; j < N; ++j)
        {
            double re_part_psi_in = func_in[j][0];
            double im_part_psi_in = func_in[j][1];
            const double U = V(coordinate[j])-coordinate[j]*E_t(t[i]);
            if(coordinate[i] < Z_abs) {
                func_in[j][0] = re_part_psi_in * cos(-U*dt) -
                        im_part_psi_in * sin(-U*dt);
                func_in[j][1] = re_part_psi_in * sin(-U*dt) +
                        im_part_psi_in * cos(-U*dt);
            } else {
                func_in[j][0] = (re_part_psi_in * cos(-U*dt) -
                        im_part_psi_in * sin(-U*dt))*CAP;
                func_in[j][1] = (re_part_psi_in * sin(-U*dt) +
                        im_part_psi_in * cos(-U*dt))*CAP;
            }
        }

        fftw_execute(plan_fwd);

        for(int j = 0; j < N; ++j)
        {
            double re_part_psi_out = func_out[j][0];
            double im_part_psi_out = func_out[j][1];
            func_out[j][0] = re_part_psi_out * cos(-p[j]*p[j]*dt/2.) -
                    im_part_psi_out * sin(-p[j]*p[j]*dt/2.);
            func_out[j][1] = re_part_psi_out * sin(-p[j]*p[j]*dt/2.) +
                    im_part_psi_out * cos(-p[j]*p[j]*dt/2.);
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
                    func_in[j][1]*func_in[j][1])*grad_V(coordinate[j])*dz;
        }

        a_t[i] = -E_t(t[i]) - Integral_sqrpsi_gradV;// на каждом шаге по времени находим дипольное ускорение
        array_E_t[i] = E_t(t[i]);

        //Записываем значение волновой функции в различные моменты времени
        if( i == M/4)
        {
            std::ofstream print_psi("psi_t_025max.dat");
            for(int j = 0; j < N; ++j)
            {
                double x = coordinate[j];
                print_psi << x << '\t' <<
                      pow(func_in[j][0],2) + pow(func_in[j][1],2) << std::endl;
            }
            print_psi.close();
        }
        else if ( i == (int)(0.125*M))
        {
            std::ofstream print_psi("psi_t_0125max.dat");
            for(int j = 0; j < N; ++j)
            {
                print_psi << coordinate[j] << '\t' <<
                             pow(func_in[j][0],2) + pow(func_in[j][1],2) << std::endl;
            }
            print_psi.close();
        }
        else if ( i == (int)(0.375*M))
        {
            std::ofstream print_psi("psi_t_0375max.dat");
            for(int j = 0; j < N; ++j)
            {
                print_psi << coordinate[j] << '\t' <<
                             pow(func_in[j][0],2) + pow(func_in[j][1],2) << std::endl;
            }
            print_psi.close();
        }
        else if( i == M/2 )
        {
            std::ofstream print_psi_mid("psi_t_0500max.dat");
            for(int j = 0; j < N; ++j)
            {
                print_psi_mid << coordinate[j] << '\t' <<
                      pow(func_in[j][0],2.) + pow(func_in[j][1],2) << '\t' <<
                      func_in[j][0] << '\t' << func_in[j][1] << std::endl;
            }
            print_psi_mid.close();
        }
        else if( i == (int)0.625*M)
        {
           std::ofstream print_psi("psi_t_0625max.dat");
           for(int j = 0; j < N; ++j)
           {
               print_psi << coordinate[j] << '\t' <<
                       pow(func_in[j][0],2) + pow(func_in[j][1],2) << std::endl;
           }
           print_psi.close();
        }
        else if( i == 0.750*M)
        {
            std::ofstream print_psi("psi_t_0750max.dat");
            for(int j = 0; j < N; ++j)
            {
                print_psi << coordinate[j] << '\t' <<
                        pow(func_in[j][0],2) + pow(func_in[j][1],2) << std::endl;
            }
            print_psi.close();
        }
        else if( i == (int)(0.875*M) )
        {
            std::ofstream print_psi("psi_t_0875max.dat");
            for(int j = 0; j < N; ++j)
            {
                print_psi << coordinate[j] << '\t' <<
                        pow(func_in[j][0],2) + pow(func_in[j][1],2) << std::endl;
            }
            print_psi.close();
        }
        else if( i == (M-1) )
        {
            std::ofstream print_psi_end("psi_t_end.dat");
            for(int j = 0; j < N; ++j)
            {
                print_psi_end << coordinate[j] << '\t' <<
                      pow(func_in[j][0],2.) + pow(func_in[j][1],2) << '\t' <<
                      func_in[j][0] << '\t' << func_in[j][1] << std::endl;
            }
            print_psi_end.close();

            double norma = 0;
            for(int j = 0; j < N; ++j)
            {
                norma += (pow(func_in[j][0],2) + pow(func_in[j][1],2))*dz;
            }
            std::cout << norma << std::endl;
        }
#if 0
        for(int j = 0; j < N; ++j)
        {
            Integral_sqrpsi_z += (func_in[j][0]*func_in[j][0]+
                    func_in[j][1]*func_in[j][1])*coordinate[j]*dz;
        }

        dipole_moment[i] = Integral_sqrpsi_z;
        //if(i%100 == 0) std::cout << dipole_moment[i] << std::endl; //1.04.2023
#endif
    }
    //return 0; // 1.04.2023
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

    //запись в файл аналитического выражения в.ф. гармонического осциллятора в момент времени t_max
    std::ofstream print_psi_osc("Analitical_Psi_osc.dat");
    for(int i = 0; i < N; ++i)
    {

        print_psi_osc << coordinate[i] << '\t' << Psi_ground_state_osc(coordinate[i])*cos((omega_L*t_max)/2) <<
                         '\t' << -Psi_ground_state_osc(coordinate[i])*sin((omega_L*t_max)/2);
    }
#if 0
    std::ofstream _out_("a_t.dat");
    for(int i = 0; i < M; ++i)
    {
        _out_ << t[i] << "\t" << a_t[i] << std::endl;
    }
    _out_.close();
    //return 0;   //временный выход из программы 22.03.2023
#endif

    //Считаем спектр дипольного ускорения
    fftw_complex *a_in = new fftw_complex[M];
    fftw_complex *a_omega = new fftw_complex[M];
    fftw_plan plan_fwd_a = fftw_plan_dft_1d(M, a_in, a_omega, FFTW_FORWARD, FFTW_MEASURE);

    for(int i = 0; i < M; ++i)
    {
        a_in[i][0] = a_t[i]*Mask(t[i],t_max);  //Домножаем дипольное ускорение на маску
        a_in[i][1] = 0.;
    }

    fftw_execute(plan_fwd_a);
#if 0
    std::ofstream print_acceleration_mult_mask("Mask_a_t.dat");
    for(int i = 0; i < M; ++i)
    {
        print_acceleration_mult_mask << t[i] << '\t' << a_in[i][0] << std::endl;
        //print_acceleration_mult_mask << t[i] << '\t' << Mask(t[i],t_max) << std::endl;
    }
#endif
    //return 0;   //26.03.2023

    delete [] a_t;
    a_t = nullptr;

    //Считаем спектр внешнего электрического поля
    fftw_complex *E_in = new fftw_complex[M];
    fftw_complex *E_omega = new fftw_complex[M];
    fftw_plan plan_fwd_E = fftw_plan_dft_1d(M, E_in, E_omega, FFTW_FORWARD, FFTW_MEASURE);

    for(int i = 0; i < M; ++i)
    {
        E_in[i][0] = array_E_t[i];
        E_in[i][1] = 0.;
    }

    fftw_execute(plan_fwd_E);

#if 0
    std::ofstream print_omega("a_omega.dat");
    for(int i = 0; i < M; ++i)
    {
        print_omega << a_omega[i][0] << '\t' << a_omega[i][1] << std::endl;
    }
    print_omega.close();
#endif

    //заполняем массив частот omega
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
//#endif
    std::ofstream print_a_omega("a_omega.dat");
    for(int i = M/2; i < M; ++i)
    {
        print_a_omega << omega[i-M/2] <<'\t' << std::sqrt(a_omega[i][0]*a_omega[i][0]+
                a_omega[i][1]*a_omega[i][1]) << std::endl;
        //print_a_omega << omega[i-M/2] <<'\t' << a_omega[i][0]*a_omega[i][0]+
                //a_omega[i][1]*a_omega[i][1] << std::endl;
    }
    for(int i = 0; i < M/2; ++i)
    {
        print_a_omega << omega[i+M/2] <<'\t' << std::sqrt(a_omega[i][0]*a_omega[i][0]+
                a_omega[i][1]*a_omega[i][1]) << std::endl;
        //print_a_omega << omega[i+M/2] <<'\t' << std::sqrt(a_omega[i][0]*a_omega[i][0]+
                //a_omega[i][1]*a_omega[i][1]) << std::endl;
    }
    print_a_omega.close();
#endif
    //Записываем в файл квадрат спектра внешнего электрического поля
    std::ofstream print_spectra_E("E_omega.dat");
    for(int i = M/2; i < M; ++i)
    {
        print_spectra_E << omega[i-M/2]/omega_L << '\t' << E_omega[i][0]*E_omega[i][0]+
                                           E_omega[i][1]*E_omega[i][1] << std::endl;
    }
    for(int i = 0; i < M/2; ++i)
    {
        print_spectra_E << omega[i+M/2]/omega_L << '\t' << E_omega[i][0]*E_omega[i][0]+
                                           E_omega[i][1]*E_omega[i][1] << std::endl;
    }
    print_spectra_E.close();

#if 0
    //Считаем спектр дипольного момента
    fftw_complex *dip_mom = new fftw_complex[M];
    fftw_complex *dm_omega = new fftw_complex[M];
    fftw_plan plan_fwd_dm = fftw_plan_dft_1d(M, dip_mom, dm_omega, FFTW_FORWARD, FFTW_MEASURE);

    for(int i = 0; i < M; ++i)
    {
        dip_mom[i][0] = dipole_moment[i];
        dip_mom[i][1] = 0.;
    }

    fftw_execute(plan_fwd_dm);

    //Записываем в файл квадрат спектра дипольного момента
    std::ofstream print_dm_omega("psi_x_psi_omega.dat");
    for(int i = M/2; i < M; ++i)
    {
        print_dm_omega << omega[i-M/2]/omega_L << "\t" << dm_omega[i][0]*dm_omega[i][0]+
                dm_omega[i][1]*dm_omega[i][1] << std::endl;
    }
    for(int i = 0; i < M/2; ++i)
    {
        print_dm_omega << omega[i+M/2]/omega_L << '\t' << dm_omega[i][0]*dm_omega[i][0]+
                                           dm_omega[i][1]*dm_omega[i][1] << std::endl;
    }
    print_dm_omega.close();
#endif

    //return 0; //new 3.12.2022
    double *HHG = new double[M];
    for(int i = 0; i < M; ++i)
    {
        HHG[i] = a_omega[i][0]*a_omega[i][0]+a_omega[i][1]*a_omega[i][1];//по определению спектра
    }

    //Записываем в файл квадрат спектра дипольного ускорения
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

    //fftw_destroy_plan(plan_fwd_dm);
    fftw_destroy_plan(plan_fwd_E);
    fftw_destroy_plan(plan_fwd_a);

    delete [] HHG;
    HHG = nullptr;
    //delete [] dip_mom;
   // dip_mom = nullptr;
    //delete [] dm_omega;
    //dm_omega = nullptr;
    delete [] E_in;
    E_in = nullptr;
    delete [] E_omega;
    E_omega = nullptr;
    delete [] omega;
    omega = nullptr;
    delete [] a_in;
    a_in = nullptr;
    delete [] a_omega;
    a_omega = nullptr;
    delete [] t;
    t = nullptr;
    delete [] p;
    p = nullptr;

    return 0;
}
