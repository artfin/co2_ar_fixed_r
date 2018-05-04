#include <mpi/mpi.h>

#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include "constants.hpp"
#include "trajectory.hpp"
#include "matrix_euler.hpp"
#include "co2_ar_dipole.hpp"

#include "awp.hpp"
#include "vmblock.hpp"
#include "gear.hpp"
#include "basis.hpp"

const int DIM = 10;
const double RDIST = 40.1;
const double sampling_time = 200.0;
const int MaxTrajectoryLength = 16384;

void syst( REAL t, REAL * y, REAL * f )
{
    (void)(t);

    double * out = new double [10];

    rhs(out, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9] );

    f[0] = out[0]; // dR / dt
    f[1] = out[1]; // dpR / dt
    f[2] = out[2]; // d Theta / dt
    f[3] = out[3]; // d pTheta / dt
    f[4] = out[4]; // d phi / dt
    f[5] = out[5]; // d p_phi / dt
    f[6] = out[6]; // d theta / dt
    f[7] = out[7]; // d p_theta / dt
    f[8] = out[8]; // d psi / dt
    f[9] = out[9]; // d p_psi / dt

    delete [] out;
}


void read_initial_conditions( std::vector<std::vector<double>> & contents, const std::string filename,
                              int lines_to_read = -1 )
{
    std::ifstream inFile(filename);
    if ( !inFile )
        throw std::invalid_argument( "Can't find file with initial conditions." );

    const int MAXLINE = 300;
    char buf[MAXLINE];

    std::stringstream ss;
    std::vector<double> temp;
    double value;

    int lines_counter = 0;
    while ( inFile.getline(buf, MAXLINE) )
    {
        // воспринимаем строчку с # как комментарий
        std::string line(buf);
        size_t found = line.find('#');
        if ( found != std::string::npos )
            continue;

        ss << buf;
        for ( int k = 0; k < DIM; ++k )
        {
            ss >> value;
            temp.push_back(value);
        }

        if ( lines_counter == lines_to_read )
            break;

        contents.push_back(temp);
        temp.clear();
        ss.clear();
        ss.str("");

        ++lines_counter;
    }

    inFile.close();
}

void create_frequencies_vector( std::vector<double> & freqs, const double sampling_time )
{
    const double FREQ_MAX = 1000.0;
    double FREQ_STEP = 1.0 / (sampling_time * constants::ATU) / constants::LIGHTSPEED_CM / MaxTrajectoryLength; // cm^-1
    int FREQ_SIZE = (int) FREQ_MAX / FREQ_STEP + 1;

    freqs.resize( FREQ_SIZE );
    for( int k = 0; k <  FREQ_SIZE; ++k )
        freqs[k] = k * FREQ_STEP;
}

//structure of array after Fourier transform used in GSL
void compl_mod( std::vector<double> & farr )
{
    int N = (int) farr.size();

    farr[0] = farr[0] * farr[0];
    farr[N/2] = farr[N/2] * farr[N/2];

    for (int i = 1; i < N/2; ++i)
        farr[i] = farr[i] * farr[i] + farr[N-i] * farr[N-i];
}

void slave_code ( int world_rank )
{
    std::vector<double> freqs;
    create_frequencies_vector( freqs, sampling_time );
    std::vector<double> spectralFunction( freqs.size() );

    Trajectory trajectory( DIM, sampling_time, MaxTrajectoryLength, RDIST );

    bool exit_status;
    std::vector<double> dipx(MaxTrajectoryLength);
    std::vector<double> dipy(MaxTrajectoryLength);
    std::vector<double> dipz(MaxTrajectoryLength);

    while ( true )
    {
        exit_status = trajectory.receive_initial_conditions();
        if ( exit_status )
        {
            std::cout << "(" << world_rank << ") Received exit tag." << std::endl;
            break;
        }

        double pr_mu = trajectory.get_pr_mu();

        trajectory.run_trajectory( syst, dipx, dipy, dipz );

        gsl_fft_real_radix2_transform( &dipx[0], 1, MaxTrajectoryLength);
        gsl_fft_real_radix2_transform( &dipy[0], 1, MaxTrajectoryLength);
        gsl_fft_real_radix2_transform( &dipz[0], 1, MaxTrajectoryLength);

        compl_mod( dipx );
        compl_mod( dipy );
        compl_mod( dipz );

        for ( size_t k = 0; k < freqs.size(); ++k )
            spectralFunction[k] = pr_mu * (dipx[k] + dipy[k] + dipz[k]) * sampling_time * sampling_time / 2.0 / M_PI;

        MPI_Send( &spectralFunction[0], (int) spectralFunction.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );

        std::fill( dipx.begin(), dipx.end(), 0.0 );
        std::fill( dipy.begin(), dipy.end(), 0.0 );
        std::fill( dipz.begin(), dipz.end(), 0.0 );

    }

    std::cout << "(" << world_rank << ") Message after breaking out of cycle" << std::endl;
}

void master_code( int world_size )
{
    MPI_Status status;
    int source;

    std::vector<double> sample;
    std::vector<std::vector<double>> samples;
    const int samples_to_read = 100;

    read_initial_conditions(samples, "./input/starting_points_-+.txt", samples_to_read);
    std::cout << "Samples len: " << samples.size() << std::endl;

    int sent = 0;
    int received = 0;

    // status of calculation
    bool is_finished = false;

    std::vector<double> freqs;
    create_frequencies_vector( freqs, sampling_time );

    // sending first trajectory
    for ( int i = 1; i < world_size; i++ )
    {
        sample = samples[sent];

        MPI_Send( &sample[0], DIM, MPI_DOUBLE, i, 0, MPI_COMM_WORLD );
        //std::cout << "(master) sent point " << std::endl;

        ++sent;
    }

    std::vector<double> specfunc_total( freqs.size() );
    std::vector<double> specfunc_package( freqs.size() );

    //double Volume = 4.0 / 3.0 * M_PI * pow( RDIST, 3);

    while( true )
    {
        if ( is_finished )
        {
            for ( int i = 1; i < world_size; i++ )
                MPI_Send( &is_finished, 1, MPI_INT, i, tags::EXIT_TAG, MPI_COMM_WORLD );

            break;
        }

        std::fill( specfunc_package.begin(), specfunc_package.end(), 0.0 );
        MPI_Recv( &specfunc_package[0], (int) specfunc_package.size(), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
        source = status.MPI_SOURCE;

        for ( size_t i = 0; i < specfunc_package.size(); ++i )
            specfunc_total[i] += specfunc_package[i];

        ++received;

        if ( received % 10 == 0 && received != 0 )
            std::cout << "(master) received = " << received << std::endl;

        if ( received == (int) samples.size() )
        {
            std::string filename = "./spectral_function.txt";
            std::cout << "(final) Saving spectral function to " << filename << std::endl;

            std::ofstream out(filename);
            for ( size_t k = 0; k < specfunc_total.size(); ++k )
            {
                //20106.2 is the ratio of  V * \frac{\int exp(-H/kt) dpR dpTheta dpPhi dTheta dPhi}
                //                                  {\int exp(-H/kT) dR dpR dpTheta dTheta dPhi}
                specfunc_total[k] = specfunc_total[k] * constants::HTOJ * std::pow(constants::ALU, 6.0) * constants::ATU * \
                        20106.2 / received * 1.0E19;
                out << freqs[k] << " " << specfunc_total[k] << std::endl;
            }

            out.close();


            std::cout << "(final) Exiting..." << std::endl;
            is_finished = true;
        }

        if ( sent < (int) samples.size() )
        {
            sample = samples[sent];
            MPI_Send( &sample[0], DIM, MPI_DOUBLE, source, 0, MPI_COMM_WORLD );

            ++sent;
        }
    }

    std::cout << "(master) out of main cycle" << std::endl;
}

int main( int argc, char * argv[] )
{
    //Initialize the MPI environment
    MPI_Init( &argc, &argv );

    //getting id of the current process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //getting number of running processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if ( world_rank == 0 )
    {
        std::clock_t start = clock();

        master_code( world_size );

        std::cout << "Time elapsed: " << (clock() - start) / (double) CLOCKS_PER_SEC << "s" << std::endl;
    }
    else
    {
        slave_code( world_rank );
    }

    MPI_Finalize();

    return 0;
}

