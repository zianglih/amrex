#include <AMReX_Geometry.H>
#include <AMReX_TagBox.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_CoordSys.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_Cluster.H>
#include <AMReX_LevelBld.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_PROB_AMR_F.H>
#include <AMReX_Amr.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabSet.H>
#include <AMReX_StateData.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <HMMAmr.H>

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#ifdef BL_USE_ARRAYVIEW
#include <DatasetClient.H>
#endif

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrInSituBridge.H>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <limits>
#include <list>
#include <sstream>

namespace amrex {

namespace
{
    const std::string CheckPointVersion("CheckPointVersion_1.0");

    //   
    // These are ParmParse'd in. Defaults set in base class
    //   
    int  regrid_on_restart;
    int  plotfile_on_restart;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// constructors/destructors: 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HMMAmr::HMMAmr(LevelBld* a_levelbld) 
    : Amr (a_levelbld)
{
}

HMMAmr::HMMAmr (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord,
          LevelBld* a_levelbld)
    : Amr(rb, max_level_in, n_cell_in, coord, a_levelbld)
{
}

HMMAmr::~HMMAmr ()
{
    levelbld->variableCleanUp();
    Amr::Finalize();
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// functions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void 
HMMAmr::timeStep_init(int level, 
                        Real time, 
                        Real stop_time)
{
    BL_PROFILE("HMMAmr::timeStep_init()");
    BL_COMM_PROFILE_NAMETAG("HMMAmr::timeStep_init TOP");

    if (verbose > 0)
    {
        amrex::Print() << "[Level " << level << " step " << level_steps[level]+1 << "] "
                       << "timeStep_init at time = " << time << "\n";
    }

    // This is used so that the AmrLevel functions can know which level is being advanced
    //      when regridding is called with possible lbase > level.
    which_level_being_advanced = level;

    // Update so that by default, we don't force a post-step regrid.
    amr_level[level]->setPostStepRegrid(0);

    //
    // Allow regridding of level 0 calculation on restart.
    //
    if (max_level == 0 && regrid_on_restart)
    {
        regrid_level_0_on_restart();
    }
    else
    {
        int lev_top = std::min(finest_level, max_level-1);
        for (int i(level); i <= lev_top; ++i)
        {
            const int old_finest = finest_level;

            if (okToRegrid(i))
            {
                regrid(i,time);

                //
                // Compute new dt after regrid if at level 0 and compute_new_dt_on_regrid.
                //
                if ( compute_new_dt_on_regrid && (i == 0) )
                {
                    int post_regrid_flag = 1;
                    amr_level[0]->computeNewDt(finest_level,
                                               sub_cycle,
                                               n_cycle,
                                               ref_ratio,
                                               dt_min,
                                               dt_level,
                                               stop_time,
                                               post_regrid_flag);
                }

                for (int k(i); k <= finest_level; ++k) {
                    level_count[k] = 0;
                }

                if (old_finest < finest_level)
                {
                    //
                    // The new levels will not have valid time steps
                    // and iteration counts.
                    //
                    for (int k(old_finest + 1); k <= finest_level; ++k)
                    {
                        dt_level[k]    = dt_level[k-1]/n_cycle[k];
                    }
                }
            }
            if (old_finest > finest_level) {
                lev_top = std::min(finest_level, max_level - 1);
            }
        }

        if (max_level == 0 && loadbalance_level0_int > 0 && loadbalance_with_workestimates)
        {
            if (level_steps[0] == 1 || level_count[0] >= loadbalance_level0_int) {
                LoadBalanceLevel0(time);
                level_count[0] = 0;
            }
        }
    }
    //
    // Check to see if should write plotfile.
    // This routine is here so it is done after the restart regrid.
    //
    if (plotfile_on_restart && ! (restart_chkfile.empty()) )
    {
        plotfile_on_restart = 0;
        writePlotFile();
    }

    //
    // Run for grids at higher level.
    //
    if (level < finest_level)
    {
        const int lev_fine = level+1;
        timeStep_init(lev_fine,time,stop_time);
    }
    
    // Set this back to negative so we know whether we are in fact in this routine
    which_level_being_advanced = -1;

}

void
HMMAmr::timeStep_BC (int  level,
               Real time,
               int  iteration,
               int  niter)
{
    BL_PROFILE("HMMAmr::timeStep_BC()");
    BL_COMM_PROFILE_NAMETAG("HMMAmr::timeStep_BC TOP");

    // This is used so that the AmrLevel functions can know which level is being advanced
    //      when regridding is called with possible lbase > level.
    which_level_being_advanced = level;

    //
    // Advance grids at this level.
    //
    if (verbose > 0)
    {
        amrex::Print() << "[Level " << level << " step " << level_steps[level]+1 << "] "
                       << "Berger-Collela ADVANCE with dt = " << dt_level[level] << "\n";
    }

    Real dt_new = amr_level[level]->advance(time,dt_level[level],iteration,niter);
    BL_PROFILE_REGION_STOP("amr_level.advance");

    dt_min[level] = iteration == 1 ? dt_new : std::min(dt_min[level],dt_new);

    level_steps[level]++;
    level_count[level]++;

    if (verbose > 0)
    {
        amrex::Print() << "[Level " << level << " step " << level_steps[level] << "] "
                       << "Advanced " << amr_level[level]->countCells() << " cells\n";
    }

    // If the level signified that it wants a regrid after the advance has
    // occurred, do that now.
    if (amr_level[level]->postStepRegrid()) {
        
        amrex::Error("Quitting because we should not do postStepRegrid. Can remove this error flag if needed.");
        int old_finest = finest_level;

        regrid(level, time);

        if (old_finest < finest_level)
        {
            //
            // The new levels will not have valid time steps.
            //
            for (int k = old_finest + 1; k <= finest_level; ++k)
            {
                dt_level[k] = dt_level[k-1] / n_cycle[k];
            }
        }
    }

    //
    // Advance grids at higher level.
    //
    if (level < finest_level)
    {
        const int lev_fine = level+1;

        if (sub_cycle)
        {
            const int ncycle = n_cycle[lev_fine];

            BL_COMM_PROFILE_NAMETAG("HMMAmr::timeStep_BC timeStep_BC subcycle");
            for (int i = 1; i <= ncycle; i++)
                timeStep_BC(lev_fine,time+(i-1)*dt_level[lev_fine],i,ncycle);
        }
        else
        {
            BL_COMM_PROFILE_NAMETAG("HMMAmr::timeStep_BC timeStep_BC nosubcycle");
            timeStep_BC(lev_fine,time,1,1);
        }
    }

    // 
    // Call the post_timestep function for the BC algorithm
    // -- this will invoke refluxing and volume averaging 
    //
    amr_level[level]->post_timestep(iteration);

    // Set this back to negative so we know whether we are in fact in this routine
    which_level_being_advanced = -1;
}

void
HMMAmr::timeStep_react (int  level,
                       Real time)
{
    BL_PROFILE("HMMAmr::timeStep_react()");
    BL_COMM_PROFILE_NAMETAG("HMMAmr::timeStep_react TOP");

    // This is used so that the AmrLevel functions can know which level is being advanced
    //      when regridding is called with possible lbase > level.
    which_level_being_advanced = level;

    //
    // Advance grids at this level.
    //
    if (verbose > 0)
    {
        amrex::Print() << "[Level " << level << " step " << level_steps[level]+1 << "] "
                       << "Reaction ADVANCE with dt = " << dt_level[0]/2.0 << "\n";
    }

    // we use half the coarse (level 0) timestep for all levels here. 
    Real dt_new = amr_level[level]->advance_react(time,dt_level[0]/2.0);

    if (verbose > 0)
    {
        amrex::Print() << "[Level " << level << " step " << level_steps[level]+1 << "] "
                       << "Reacted " << amr_level[level]->countCells() << " cells\n";
    }

    //
    // Advance grids at higher level.
    //
    if (level < finest_level)
    {
        const int lev_fine = level+1;
        BL_COMM_PROFILE_NAMETAG("HMMAmr::timeStep_react");
        timeStep_react(lev_fine,time);
    }

    // 
    // Call the post_timestep function
    //
    amr_level[level]->post_timestep_react(1);

    // Set this back to negative so we know whether we are in fact in this routine
    which_level_being_advanced = -1;
}
void
HMMAmr::coarseTimeStepMod (Real stop_time)
{
    double run_stop;
    double run_strt;
    BL_PROFILE_REGION_START("HMMAmr::coarseTimeStepMod()");
    BL_PROFILE("HMMAmr::coarseTimeStepMod()");
    std::stringstream stepName;
    stepName << "timeStep STEP " << level_steps[0];

    run_strt = amrex::second() ;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Compute new dt.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (levelSteps(0) > 0)
    {
        int post_regrid_flag = 0;
        amr_level[0]->computeNewDt(finest_level,
                                   sub_cycle,
                                   n_cycle,
                                   ref_ratio,
                                   dt_min,
                                   dt_level,
                                   stop_time,
                                   post_regrid_flag);
    }
    else
    {
        amr_level[0]->computeInitialDt(finest_level,
                                       sub_cycle,
                                       n_cycle,
                                       ref_ratio,
                                       dt_level,
                                       stop_time);
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Perform the timestep: 
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BL_PROFILE_REGION_START(stepName.str());

    // ~~~~ initialization
    // a) Allow regridding of level calculation on restart.
    // b) Check to see if should write plotfile.
    // c) Set default postStepRegrid to false 
    timeStep_init(0, cumtime, stop_time);

    // ~~~~ advance reactions (1/2) 
    timeStep_react(/* level = */      0,
                   /* time = */       cumtime);

    // ~~~~ advance with Berger-Collela subcycling: 
    timeStep_BC(/* level = */      0,
                /* time = */       cumtime,
                /* iteration = */  1,
                /* niter = */      1);


    // ~~~~ advance reactions (2/2)
    timeStep_react(/* level = */      0,
                   /* time = */       cumtime);

    // ~~~~ finalization 
    BL_PROFILE_REGION_STOP(stepName.str());

    cumtime += dt_level[0];

    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Perform the post coarsetimestep routines: 
    // -- this is where we print summary diagnostics 
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    amr_level[0]->postCoarseTimeStep(cumtime);


    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        run_stop = amrex::second() - run_strt;
        const int istep    = level_steps[0];

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelDescriptor::ReduceRealMax(run_stop,IOProc);
            amrex::Print() << "\n[STEP " << istep << "] Coarse TimeStep time: " << run_stop << '\n';
#ifdef BL_LAZY
        });
#endif

#ifndef AMREX_MEM_PROFILING
        Long min_fab_kilobytes  = amrex::TotalBytesAllocatedInFabsHWM()/1024;
        Long max_fab_kilobytes  = min_fab_kilobytes;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelDescriptor::ReduceLongMin(min_fab_kilobytes, IOProc);
            ParallelDescriptor::ReduceLongMax(max_fab_kilobytes, IOProc);

            amrex::Print() << "[STEP " << istep << "] FAB kilobyte spread across MPI nodes: ["
                           << min_fab_kilobytes << " ... " << max_fab_kilobytes << "]\n";
#ifdef BL_LAZY
            amrex::Print() << "\n";
            });
#endif
#endif
    }

#ifdef AMREX_MEM_PROFILING
    {
        std::ostringstream ss;
        ss << "[STEP " << level_steps[0] << "]";
        MemProfiler::report(ss.str());
    }
#endif

    BL_PROFILE_ADD_STEP(level_steps[0]);
    BL_PROFILE_REGION_STOP("HMMAmr::coarseTimeStepMod()");
    BL_COMM_PROFILE_NAMETAG(stepName.str());
    //BL_PROFILE_FLUSH();
    BL_TRACE_PROFILE_FLUSH();
    BL_COMM_PROFILE_FLUSH();

    if (verbose > 0)
    {
        amrex::Print()
            << "\nSTEP = " << level_steps[0]
            << " TIME = "  << cumtime
            << " DT = "    << dt_level[0] << "\n\n";
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "STEP = "  << level_steps[0]
               << " TIME = " << cumtime
               << " DT = "   << dt_level[0] << '\n';
    }
    if (record_run_info_terse && ParallelDescriptor::IOProcessor())
        runlog_terse << level_steps[0] << " " << cumtime << " " << dt_level[0] << '\n';

    int check_test = 0;

    if (check_per > 0.0)
    {

        // Check to see if we've crossed a check_per interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = static_cast<int>((cumtime-dt_level[0]) / check_per);
        int num_per_new = static_cast<int>((cumtime            ) / check_per);

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next check_per interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0_rt * std::abs(cumtime);
        const Real next_chk_time = (num_per_old + 1) * check_per;

        if ((num_per_new == num_per_old) && std::abs(cumtime - next_chk_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && std::abs((cumtime - dt_level[0]) - next_chk_time) <= eps)
        {
            num_per_old += 1;
        }

        if (num_per_old != num_per_new)
        {
            check_test = 1;
        }

    }

    int to_stop       = 0;
    int to_checkpoint = 0;
    int to_plot       = 0;
    int to_small_plot = 0;
    if (message_int > 0 && level_steps[0] % message_int == 0) {
        if (ParallelDescriptor::IOProcessor())
        {
            FILE *fp;
            if ((fp=fopen("dump_and_continue","r")) != 0)
            {
                remove("dump_and_continue");
                to_checkpoint = 1;
                fclose(fp);
            }
            else if ((fp=fopen("stop_run","r")) != 0)
            {
                remove("stop_run");
                to_stop = 1;
                fclose(fp);
            }
            else if ((fp=fopen("dump_and_stop","r")) != 0)
            {
                remove("dump_and_stop");
                to_checkpoint = 1;
                to_stop = 1;
                fclose(fp);
            }

            if ((fp=fopen("plot_and_continue","r")) != 0)
            {
                remove("plot_and_continue");
                to_plot = 1;
                fclose(fp);
            }

            if ((fp=fopen("small_plot_and_continue","r")) != 0)
            {
                remove("small_plot_and_continue");
                to_small_plot = 1;
                fclose(fp);
            }
        }
        int packed_data[4];
        packed_data[0] = to_stop;
        packed_data[1] = to_checkpoint;
        packed_data[2] = to_plot;
        packed_data[3] = to_small_plot;
        ParallelDescriptor::Bcast(packed_data, 4, ParallelDescriptor::IOProcessorNumber());
        to_stop = packed_data[0];
        to_checkpoint = packed_data[1];
        to_plot = packed_data[2];
        to_small_plot = packed_data[3];

    }

    if(to_stop == 1 && to_checkpoint == 0) {  // prevent main from writing files
        last_checkpoint = level_steps[0];
        last_plotfile   = level_steps[0];
    }

    if (to_checkpoint && write_plotfile_with_checkpoint) {
        to_plot = 1;
        to_small_plot = 1;
    }

    if ((check_int > 0 && level_steps[0] % check_int == 0) || check_test == 1
        || to_checkpoint)
    {
        checkPoint();
    }


    if (writePlotNow() || to_plot)
    {
        writePlotFile();
    }

    if (writeSmallPlotNow() || to_small_plot)
    {
        writeSmallPlotFile();
    }

    updateInSitu();

    bUserStopRequest = to_stop;
    if (to_stop)
    {
        ParallelDescriptor::Barrier("HMMAmr::coarseTimeStepMod::to_stop");
        if(ParallelDescriptor::IOProcessor()) {
            if (to_checkpoint)
            {
                amrex::ErrorStream() << "Stopped by user w/ checkpoint" << std::endl;
            }
            else
            {
                amrex::ErrorStream() << "Stopped by user w/o checkpoint" << std::endl;
            }
        }
    }
}


void
HMMAmr::RegridOnlyMod (Real time, bool do_io)
{
    // BL_ASSERT(regrid_on_restart == 1);

    if (max_level == 0)
    {    
        regrid_level_0_on_restart();
    }    
    else 
    {    
        int lev_top = std::min(finest_level, max_level-1);
        
        lev_top = max_level-1;    

        for (int i = 0; i <= lev_top; i++) 
        {
           regrid(i,time);
        }
    }    

    // if (do_io) {

    //     if (plotfile_on_restart)
    //         writePlotFile();

    //     if (checkpoint_on_restart)
    //         checkPoint();

    //     if (insitu_on_restart)
    //         updateInSitu();

    // }    
}







}
