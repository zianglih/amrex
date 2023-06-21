#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

// This routine sets the initial condition for the specific problem.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    // Extract the geometry data.
    int coord_type = geomdata.Coord();
    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    // Get the x-location for cell i
    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);
    Real y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
    Real z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);

    // Our state variables to set:
    // -- density: C_RHO
    // -- momentum: C_RHOUX,Y,Z
    // -- total energy: C_RHOE
    // -- temperature: C_T
    // -- species: C_RHOY
    // -- passive scalar: C_PHI

    // ~~~~ declare variables.
    

    
    state(i,j,k,C_PHI) = 0.0;
}




// This routine specifies a custom boundary condition. Leave it empty if not applicable.
// Make sure input file is configured with "inflow" option to activate this function.
// ::RAL:: Need to make another version so we can a custom in inflow and a custom injecting wall.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void bc_custom_func
(
    const amrex::Real* x,       // spatial coordinate
    const amrex::Real* s_int,   // internal cell (size NUM_STATE)
    const amrex::Real* s_int_reflect, // internal reflected cell data 
    amrex::Real* s_ext,         // ghost cell (size NUM_STATE, the output)
    const int idir,             // component direction (0 is x, 1 is y, 2 is z)
    const int sgn,              // sign (+1 means ghost cell is on left/bottom/back, -1 is right/top/front)
    const amrex::Real time,     // time
    amrex::GeometryData const& geomdata // geometry
)
{



}


#endif 
