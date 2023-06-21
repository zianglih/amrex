#include <AMReX_BLFort.H>
#include <HMM.H>
#include <BC_funcs.H>

// problem-specific function definition
#include "problem_setup.H"

using namespace amrex;

// Dummy 
struct GenericFill
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& /*iv*/, 
                        Array4<Real> const& /*data*/,
                        const int /*dcomp*/, 
                        const int /*numcomp*/,
                        GeometryData const& /*geom*/, 
                        const Real /*time*/,
                        const BCRec* /*bcr*/,
                        const int /*bcomp*/,
                        const int /*orig_comp*/) const
        {   
            // Dummy routine that does nothing for inflow boundaries.
            // We assume that there are no inflow boundaries for a
            // generic fill, since we overwrote them with first-order
            // extrapolation in variableSetUp().
        }   
};



// Takes care of the ext_dir BC fill.
// -- this is how amrex documentation tells us to fill BC
struct StateFill
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, // cell location in integer space 
                        Array4<Real> const& data, // output data 
                        const int /*dcomp*/,  
                        const int /*numcomp*/,
                        GeometryData const& geom, // geometry 
                        const Real time, // time 
                        const BCRec* bcr, // bcrec object 
                        const int /*bcomp*/,
                        const int /*orig_comp*/) const
    {    
        // THIS WILL SET EXTERNAL DIRICHLET CONDITIONS FOR CELL "iv"
        // --- if external BCs are not specified, this function does nothing.

        // Get domain params
        const int* domlo = geom.Domain().loVect();
        const int* domhi = geom.Domain().hiVect();
        const amrex::Real* prob_lo = geom.ProbLo();

        // get cell spacing
        const amrex::Real* dx = geom.CellSize();
                    
        // get bcrec data: 
        const int* bc = bcr->data(); 
        
        // get spatial coordinate
        const amrex::Real x[AMREX_SPACEDIM] = 
        {AMREX_D_DECL(  prob_lo[0] + static_cast<amrex::Real>(iv[0] + 0.5) * dx[0],
                        prob_lo[1] + static_cast<amrex::Real>(iv[1] + 0.5) * dx[1],
                        prob_lo[2] + static_cast<amrex::Real>(iv[2] + 0.5) * dx[2])}; 

        // initialize the boundary data 
        amrex::Real s_int[NUM_STATE] = {0.0}; // internal cell 
        amrex::Real s_int_reflect[NUM_STATE] = {0.0}; // internal reflected cell 
        amrex::Real s_ext[NUM_STATE] = {0.0}; // external (ghost) cell -- this is the output 
        amrex::IntVect loc_reflect(AMREX_D_DECL(0,0,0)); 

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // ~~~~ x-direction 
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        int idir = 0; 
        // if the BC is external type and cell "iv" is a ghost cell on left part of domain: 
        if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir]))
        {
            // Get the location of interior boundary cell in x-direction  
            amrex::IntVect loc(AMREX_D_DECL(domlo[idir], iv[1], iv[2]));
            
            // fill s_int with interior cell  
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int[n] = data(loc, n); // set the internal variable 
            }

            // Get the location of interior reflected cell in x-direction 
            // if currrent cell is -1, reflected cell is 0
            if (iv[idir] == domlo[idir]-1)
            {
                // loc_reflect = {AMREX_D_DECL(domlo[idir], iv[1], iv[2])};
                loc_reflect[0] = domlo[idir];
                loc_reflect[1] = iv[1];
                loc_reflect[2] = iv[2];

            }
            // if current cell is -2, reflected cell is 1 
            else if (iv[idir] == domlo[idir]-2)
            {
                // loc_reflect = {AMREX_D_DECL(domlo[idir]+1, iv[1], iv[2])};
                loc_reflect[0] = domlo[idir]+1;
                loc_reflect[1] = iv[1];
                loc_reflect[2] = iv[2];

            }
            // if current cell is -3, reflected cell is 2 
            else if (iv[idir] == domlo[idir]-3)
            {
                // loc_reflect = {AMREX_D_DECL(domlo[idir]+2, iv[1], iv[2])};
                loc_reflect[0] = domlo[idir]+2;
                loc_reflect[1] = iv[1];
                loc_reflect[2] = iv[2];

            }
            else
            {
                amrex::Error("SB: Error in StateFill function (BC_funcs.cpp). Update stencil to accomodate NUM_GHOST > 3");
            }

            // populate the reflected internal state: 
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int_reflect[n] = data(loc_reflect, n); // set the internal variable 
            }

            // Fill s_ext (ghost cell) 
            // bc_custom_func(/* spatial coordinate */         x,
            //                /* internal cell data */     s_int,
            //                /* ghost cell data (out) */  s_ext,
            //                /* component direction */     idir,
            //                /* sgn */                      +1,
            //                /* time */                    time,
            //                /* geometry */                geom);

            bc_custom_func(/* spatial coordinate */         x,
                           /* internal cell data */     s_int,
                           /* internal reflected cell data */     s_int_reflect,
                           /* ghost cell data (out) */  s_ext,
                           /* component direction */     idir,
                           /* sgn */                      +1,
                           /* time */                    time,
                           /* geometry */                geom);

            // copy into array4  
            for (int n = 0; n < NUM_STATE; n++)
            {
                data(iv, n) = s_ext[n];
            }

            // amrex::Print() << "SB: in BC func for cell " << iv << std::endl;

        } 
        // else, if the BC is external type and cell "iv" is a ghost cell on the right part of domain 
        else if ( (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) && (iv[idir] > domhi[idir]))
        {
            // Get the location of interior boundary cell in x-direction  
            amrex::IntVect loc(AMREX_D_DECL(domhi[idir], iv[1], iv[2]));

            // fill s_int with interior cell
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int[n] = data(loc, n);
            }

            // Get the location of interior reflected cell in x-direction 
            if (iv[idir] == domhi[idir]+1)
            {
                // loc_reflect = {AMREX_D_DECL(domhi[idir], iv[1], iv[2])};
                loc_reflect[0] = domhi[idir];
                loc_reflect[1] = iv[1];
                loc_reflect[2] = iv[2];
            }
            else if (iv[idir] == domhi[idir]+2)
            {
                // loc_reflect = {AMREX_D_DECL(domhi[idir]-1, iv[1], iv[2])};
                loc_reflect[0] = domhi[idir]-1;
                loc_reflect[1] = iv[1];
                loc_reflect[2] = iv[2];
            }
            else if (iv[idir] == domhi[idir]+3)
            {
                loc_reflect[0] = domhi[idir]-2;
                loc_reflect[1] = iv[1];
                loc_reflect[2] = iv[2];
            }
            else
            {
                amrex::Error("SB: Error in StateFill function (BC_funcs.cpp). Update stencil to accomodate NUM_GHOST > 3");
            }

            // populate the reflected internal state: 
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int_reflect[n] = data(loc_reflect, n); // set the internal variable 
            }



            // Fill s_ext (ghost cell)
            bc_custom_func(/* spatial coordinate */         x,
                           /* internal cell data */     s_int,
                           /* internal reflected cell data */     s_int_reflect,
                           /* ghost cell data (out) */  s_ext,
                           /* component direction */     idir,
                           /* sgn */                      -1,
                           /* time */                    time,
                           /* geometry */                geom);

            // copy into array4  
            for (int n = 0; n < NUM_STATE; n++) 
            {
                data(iv, n) = s_ext[n];
            }
        }
        
        #if AMREX_SPACEDIM >= 2 
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // ~~~~ y-direction 
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        idir = 1; 
        // if the BC is external type and cell "iv" is a ghost cell on left part of domain: 
        if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir]))
        {
            // Get the location of interior boundary cell in x-direction  
            amrex::IntVect loc(AMREX_D_DECL(iv[0], domlo[idir], iv[2]));

            // fill s_int with interior cell  
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int[n] = data(loc, n); // set the internal variable 
            }


            // Get the location of interior reflected cell in y-direction 
            // if currrent cell is -1, reflected cell is 0
            if (iv[idir] == domlo[idir]-1)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domlo[idir], iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = domlo[idir];
                loc_reflect[2] = iv[2];
            }
            // if current cell is -2, reflected cell is 1 
            else if (iv[idir] == domlo[idir]-2)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domlo[idir]+1, iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = domlo[idir]+1;
                loc_reflect[2] = iv[2];
            }
            // if current cell is -3, reflected cell is 2 
            else if (iv[idir] == domlo[idir]-3)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domlo[idir]+2, iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = domlo[idir]+2;
                loc_reflect[2] = iv[2];
            }
            else
            {
                amrex::Error("SB: Error in StateFill function (BC_funcs.cpp). Update stencil to accomodate NUM_GHOST > 3");
            }

            // populate the reflected internal state: 
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int_reflect[n] = data(loc_reflect, n); // set the internal variable 
            }

            // Fill s_ext (ghost cell) 
            bc_custom_func(/* spatial coordinate */         x,
                           /* internal cell data */     s_int,
                           /* internal reflected cell data */     s_int_reflect,
                           /* ghost cell data (out) */  s_ext,
                           /* component direction */     idir,
                           /* sgn */                      +1,
                           /* time */                    time,
                           /* geometry */                geom);

            // copy into array4  
            for (int n = 0; n < NUM_STATE; n++)
            {
                data(iv, n) = s_ext[n];
            }
        } 
        // else, if the BC is external type and cell "iv" is a ghost cell on the right part of domain 
        else if ( (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) && (iv[idir] > domhi[idir]))
        {
            // Get the location of interior boundary cell in x-direction  
            amrex::IntVect loc(AMREX_D_DECL(iv[0], domhi[idir], iv[2]));

            // fill s_int with interior cell
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int[n] = data(loc, n);
            }

            
            // Get the location of interior reflected cell in x-direction 
            if (iv[idir] == domhi[idir]+1)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domhi[idir], iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = domhi[idir];
                loc_reflect[2] = iv[2];
            }
            else if (iv[idir] == domhi[idir]+2)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domhi[idir]-1, iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = domhi[idir]-1;
                loc_reflect[2] = iv[2];
            }
            else if (iv[idir] == domhi[idir]+3)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domhi[idir]-2, iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = domhi[idir]-2;
                loc_reflect[2] = iv[2];
            }
            else
            {
                amrex::Error("SB: Error in StateFill function (BC_funcs.cpp). Update stencil to accomodate NUM_GHOST > 3");
            }

            // populate the reflected internal state: 
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int_reflect[n] = data(loc_reflect, n); // set the internal variable 
            }

            // Fill s_ext (ghost cell)
            bc_custom_func(/* spatial coordinate */         x,
                           /* internal cell data */     s_int,
                           /* internal reflected cell data */     s_int_reflect,
                           /* ghost cell data (out) */  s_ext,
                           /* component direction */     idir,
                           /* sgn */                      -1,
                           /* time */                    time,
                           /* geometry */                geom);

            // copy into array4  
            for (int n = 0; n < NUM_STATE; n++) 
            {
                data(iv, n) = s_ext[n];
            }
        }
        #endif

        #if AMREX_SPACEDIM == 3
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // ~~~~ z-direction 
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        idir = 2; 
        // if the BC is external type and cell "iv" is a ghost cell on left part of domain: 
        if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir]))
        {
            // Get the location of interior boundary cell in x-direction  
            amrex::IntVect loc(AMREX_D_DECL(iv[0], iv[1], domlo[idir]));

            // fill s_int with interior cell  
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int[n] = data(loc, n); // set the internal variable 
            }


            // Get the location of interior reflected cell in y-direction 
            // if currrent cell is -1, reflected cell is 0
            if (iv[idir] == domlo[idir]-1)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domlo[idir], iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = iv[1];
                loc_reflect[2] = domlo[idir];
            }
            // if current cell is -2, reflected cell is 1 
            else if (iv[idir] == domlo[idir]-2)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domlo[idir]+1, iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = iv[1];
                loc_reflect[2] = domlo[idir]+1;
            }
            // if current cell is -3, reflected cell is 2 
            else if (iv[idir] == domlo[idir]-3)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domlo[idir]+2, iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = iv[1];
                loc_reflect[2] = domlo[idir]+2;
            }
            else
            {
                amrex::Error("SB: Error in StateFill function (BC_funcs.cpp). Update stencil to accomodate NUM_GHOST > 3");
            }

            // populate the reflected internal state: 
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int_reflect[n] = data(loc_reflect, n); // set the internal variable 
            }

            // Fill s_ext (ghost cell) 
            bc_custom_func(/* spatial coordinate */         x,
                           /* internal cell data */     s_int,
                           /* internal reflected cell data */     s_int_reflect,
                           /* ghost cell data (out) */  s_ext,
                           /* component direction */     idir,
                           /* sgn */                      +1,
                           /* time */                    time,
                           /* geometry */                geom);

            // copy into array4  
            for (int n = 0; n < NUM_STATE; n++)
            {
                data(iv, n) = s_ext[n];
            }
        } 
        // else, if the BC is external type and cell "iv" is a ghost cell on the right part of domain 
        else if ( (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) && (iv[idir] > domhi[idir]))
        {
            // Get the location of interior boundary cell in x-direction  
            amrex::IntVect loc(AMREX_D_DECL(iv[0], iv[1], domhi[idir]));

            // fill s_int with interior cell
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int[n] = data(loc, n);
            }

            
            // Get the location of interior reflected cell in x-direction 
            if (iv[idir] == domhi[idir]+1)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domhi[idir], iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = iv[1];
                loc_reflect[2] = domhi[idir];
            }
            else if (iv[idir] == domhi[idir]+2)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domhi[idir]-1, iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = iv[1];
                loc_reflect[2] = domhi[idir]-1;
            }
            else if (iv[idir] == domhi[idir]+3)
            {
                // loc_reflect = {AMREX_D_DECL(iv[0], domhi[idir]-2, iv[2])};
                loc_reflect[0] = iv[0];
                loc_reflect[1] = iv[1];
                loc_reflect[2] = domhi[idir]-2;
            }
            else
            {
                amrex::Error("SB: Error in StateFill function (BC_funcs.cpp). Update stencil to accomodate NUM_GHOST > 3");
            }

            // populate the reflected internal state: 
            for (int n = 0; n < NUM_STATE; n++) 
            {
                s_int_reflect[n] = data(loc_reflect, n); // set the internal variable 
            }

            // Fill s_ext (ghost cell)
            bc_custom_func(/* spatial coordinate */         x,
                           /* internal cell data */     s_int,
                           /* internal reflected cell data */     s_int_reflect,
                           /* ghost cell data (out) */  s_ext,
                           /* component direction */     idir,
                           /* sgn */                      -1,
                           /* time */                    time,
                           /* geometry */                geom);

            // copy into array4  
            for (int n = 0; n < NUM_STATE; n++) 
            {
                data(iv, n) = s_ext[n];
            }
        }


        #endif
    }   
};


void inflow_fill(const Box& bx, Array4<Real> const& state,
                        Geometry const& geom, const Vector<BCRec>& bcr)
{
    // Make a copy of the BC that can be passed by value to the ParallelFor.

    BCRec bc = bcr[0];

    const auto domlo = geom.Domain().loVect3d();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        bool inflow_x_lo = (bc.lo(0) == EXT_DIR);

        // If the left side is an inflow: 
        if (inflow_x_lo && i < domlo[0])
        {
        
            state(i,j,k) = 1.0;
        
        
        }
        else
        {
            amrex::Error("BC function not ready");
        }
    });
}


void hmm_statefill( Box const& bx, 
                    FArrayBox& data,
                    const int dcomp, 
                    const int numcomp,
                    Geometry const& geom, 
                    const Real time,
                    const Vector<BCRec>& bcr, 
                    const int bcomp,
                    const int scomp )
{
    // Here dcomp is the component in the destination array that we
    // are filling and bcr is a vector of length ncomp which are the
    // BC values corresponding to components dcomp to dcomp + ncomp -
    // 1
    
    // First, fill all the BC data using the default routines.
    // We replace inflow with outflow in the generic fill to ensure that
    // valid data is always present.
    Vector<BCRec> bcr_noinflow{bcr};
    for (int i = 0; i < bcr_noinflow.size(); ++i) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (bcr_noinflow[i].lo(dir) == EXT_DIR) {
                bcr_noinflow[i].setLo(dir, FOEXTRAP);
            }   
            if (bcr_noinflow[i].hi(dir) == EXT_DIR) {
                bcr_noinflow[i].setHi(dir, FOEXTRAP);
            }   
        }   
    }   

    // Apply a generic fill to the data.
    GpuBndryFuncFab<GenericFill> gpu_bndry_func(GenericFill{});
    gpu_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr_noinflow, bcomp, scomp);

    // As of now, dirichlet ghost cells have been filled with outflow values. Overwrite if applicable. 
    inflow_fill(bx, data.array(dcomp), geom, bcr);

    // ~~~~ // // Override with problem-specific BC fill operation 
    // ~~~~ // // the input is a FAB
    // ~~~~ // const auto state = data.array();

    // ~~~~ // // Copy BCs to an Array1D so they can be passed by value to the ParallelFor.
    // ~~~~ // Array1D<BCRec, 0, NUM_STATE - 1> bcs;
    // ~~~~ // for (int n = 0; n < numcomp; ++n) {
    // ~~~~ //     bcs(n) = bcr[n];
    // ~~~~ // }   

    // ~~~~ // const auto geomdata = geom.data();
    // ~~~~ // amrex::ParallelFor(bx,
    // ~~~~ // [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    // ~~~~ // {   
    // ~~~~ //     problem_bc_fill(i, j, k, state, time, bcs, geomdata);
    // ~~~~ // }); 
}


void hmm_statefill_v2(  Box const& bx, 
                        FArrayBox& data,
                        const int dcomp, 
                        const int numcomp,
                        Geometry const& geom, 
                        const Real time,
                        const Vector<BCRec>& bcr, 
                        const int bcomp,
                        const int scomp )

{
    GpuBndryFuncFab<StateFill> gpu_bndry_func(StateFill{});
    gpu_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}
