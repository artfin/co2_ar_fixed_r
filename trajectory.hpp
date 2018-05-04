#pragma once

#include <mpi.h>

#include <cmath>
#include <vector>

#include "matrix_euler.hpp"

#include "tags.hpp"
#include "basis.hpp"
#include "vmblock.hpp"
#include "gear.hpp"

class Trajectory
{
public:
    Trajectory( const int N, const double sampling_time, const int MaxTrajectoryLength, const double RDIST ) :
        N(N), sampling_time(sampling_time), MaxTrajectoryLength(MaxTrajectoryLength), RDIST(RDIST)
    {
        vmblock = vminit();
        y0 = (REAL*) vmalloc( vmblock, VEKTOR, N, 0 );

        if ( ! vmcomplete(vmblock) )
        {
            vmfree( vmblock ); // free memory in list
        }
    }

    ~Trajectory( )
    {
        std::cout << "Entered Trajectory destructor" << std::endl;
        vmfree( vmblock );
        std::cout << "Freed gear data" << std::endl;
    }

    bool receive_initial_conditions( void );
    void show_initial_conditions( void );

    void zero_out_dipoles( void );

    void run_trajectory( dglsysfnk syst, std::vector<double> & dipx, std::vector<double> & dipy, std::vector<double> & dipz );

    double get_pr_mu() const { return -y0[1] / 36440.0; }

private:
    REAL * y0 = 0; // pointer to [0..n-1]-vector; initial value

    int N;
    double sampling_time;
    int MaxTrajectoryLength;
    double RDIST;

    void * vmblock; // start of the vector/matrix list
};
