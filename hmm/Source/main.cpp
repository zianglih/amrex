#include <new>
#include <iostream>
#include <iomanip>

// AMR includes: 
#include <AMReX_Amr.H>
#include <HMMAmr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include <HMM.H>
#include <AMReX.H>

using namespace amrex;

amrex::LevelBld* getLevelBld ();

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    auto dRunTime1 = amrex::second();

    int  max_step;
    Real strt_time;
    Real stop_time;

    {
        ParmParse pp;

        max_step  = -1;
        strt_time =  0.0;
        stop_time = -1.0;

        pp.query("max_step",max_step);
        pp.query("strt_time",strt_time);
        pp.query("stop_time",stop_time);
    }

    if (strt_time < 0.0) {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0) {
	amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {
        
        HMMAmr amr(getLevelBld());
        amr.init(strt_time, stop_time);       

        // ~~~~ From castro: 
        // If we set the regrid_on_restart flag and if we are *not* going to take
        //    a time step then we want to go ahead and regrid here.
        if (amr.RegridOnRestart() && ( (amr.levelSteps(0) >= max_step) || (amr.cumTime() >= stop_time) ))  
        {   
            //  
            // Regrid only!
            //  
            //

            // amr.RegridOnly(amr.cumTime());
            amr.RegridOnlyMod(amr.cumTime());
            amr.writePlotFile();
            amr.checkPoint();
        }   

        // ~~~~ Running the simulation
        while ( amr.okToContinue() &&
               (amr.levelSteps(0) < max_step || max_step < 0) &&
               (amr.cumTime() < stop_time || stop_time < 0.0) )

        {
            // amr.coarseTimeStepMod(stop_time); // SB: for reaction v-cycle implementation.
            amr.coarseTimeStep(stop_time);
        }

        // Write final checkpoint and plotfile
        if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
            amr.checkPoint();
        }

        if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
            amr.writePlotFile();
        }



        // ~~~~ Sandbox: 
        // Initialize a multifab with the level 0 boxArray, and two ghost cells. Basically Sborder.
        // Do a fillPatch on this multifab at time 0. 
        // Check the boundary values. 
    
        // // Get the state data for level 0. 
        // AmrLevel& amrlev = amr.getLevel(0);        
        // MultiFab& mf = amrlev.get_new_data(0);


        // // Initialize a separate multifab:
        // // things I need: Sborder.define(grids,dmap,NUM_STATE,NUM_GROW);
        // BoxArray ba = amrlev.boxArray(); 
        // DistributionMapping dmap = amrlev.DistributionMap();
        // int nghost = 2;
        // int ncomp = 1;
        // amrex::MultiFab data_ghost(ba, dmap, ncomp, nghost);
    
        // // Do a fillpatch: 
        // amrlev.FillPatch(amrlev, data_ghost, nghost, 0.0, 0, 0, 1); 

        // // Loop through this data: 
        // for (MFIter mfi(data_ghost); mfi.isValid(); ++mfi)
        // {    
        //     const Box& bx     = mfi.fabbox(); // includes ghost cells 
        //     
        //     const amrex::Dim3 lo = amrex::lbound(bx);
        //     const amrex::Dim3 hi = amrex::ubound(bx);  
        //     amrex::FArrayBox& phi_ghost = data_ghost[mfi]; 
        //     amrex::Array4<amrex::Real> const& phi_temp_ghost = phi_ghost.array(); 

        //     for (int i = lo.x; i <= hi.x; ++i)
        //     {   
        //         printf("data[%d] = %g\n", i, phi_temp_ghost(i,0,0));
        //     }   
        // }   








    }

    Print() << "\n\n\n\n" << std::endl;
    auto dRunTime2 = amrex::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run time = " << dRunTime2 << std::endl;

    amrex::Finalize();

    return 0;
}
