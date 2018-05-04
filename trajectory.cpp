#include "trajectory.hpp"

bool Trajectory::receive_initial_conditions( void )
{
    MPI_Status status;
    MPI_Recv( y0, N, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
    //std::cout << "(slave) Received initial points" << std::endl;

    if ( status.MPI_TAG == tags::EXIT_TAG )
        return true;

    return false;
}

void Trajectory::show_initial_conditions( void )
{
    std::cout << "y0: ";
    for ( int i = 0; i != N; i++ )
        std::cout << y0[i] << " ";
    std::cout << std::endl;
}

void Trajectory::run_trajectory( dglsysfnk syst, std::vector<double> & dipx,
                                                 std::vector<double> & dipy,
                                                 std::vector<double> & dipz )
{
    REAL epsabs = 1.0E-13; // absolute error bound
    REAL epsrel = 1.0E-13; // relative error bound

    REAL t0 =  0.0;
    REAL h = 0.1; // initial, final step size
    REAL xend = sampling_time;         /* right edge of integration interval        */

    long fmax = 1e9; // maximum number of calls of right side in gear4()
    long aufrufe; // actual number of function calls
    int fehler;  // error code from gear4()

    std::vector<double> dip_(3);

    for ( int counter = 0; y0[0] < RDIST; ++counter )
    {
        if ( counter == MaxTrajectoryLength )
        {
            std::cout << "Trajectory cut" << std::endl;
            break;
        }

        fehler = gear4( &t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe );
        if ( fehler != 0 )
        {
            std::cout << "Gear4: error n = " << 10 + fehler << std::endl;
            break;
        }

        transform_dipole( dip_, y0[0], y0[2], y0[4], y0[6], y0[8] );
        dipx[counter] = dip_[0];
        dipy[counter] = dip_[1];
        dipz[counter] = dip_[2];

        xend = sampling_time * (counter + 2);
        aufrufe = 0;
    }
}

