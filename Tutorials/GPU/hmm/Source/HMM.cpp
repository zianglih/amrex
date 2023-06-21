// AMR includes: 
#include <HMM.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
// #include <AMReX.H>
#include <unistd.h>

// Include run-specific functions:
#include "problem_setup.H" 

// Include header files defining inlines
// #include "HMM_IL_timestep.H" 
using namespace amrex;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize the static members and other variables.   
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Protected variables: 
int      HMM::verbose               = 0;
Real     HMM::cfl                   = 0.9;
int      HMM::use_viscosity         = 0;
Real     HMM::dt_max                = 0.01; 
Real     HMM::dt_fixed              = -1.0; 
int      HMM::do_reflux             = 1;
int      HMM::do_subcycle           = 1;
int      HMM::num_state_type        = 0;   // this is reset in read_params 
BCRec    HMM::phys_bc;
std::vector<std::string> HMM::err_list_names;
std::vector<int> HMM::err_list_ng;
// std::vector<int> HMM::convFluxMethod; 
// std::vector<int> HMM::interpMethod; 
// std::vector<int> HMM::advanceMethod; 
// std::vector<int> HMM::do_reaction;
// int HMM::do_reaction_vcycle; 
// std::vector<int> HMM::do_convection;
// std::string HMM::reaction_subcycle_scheme = "explicit_forwardEuler"; 
// std::string HMM::reaction_algorithm = "strang"; 



// tagging: 
Vector<std::string> HMM::tag_field_list; 
Vector<int> HMM::tag_max_level_list;            // max level to tag
Vector<int> HMM::tag_ng_list;                   // ghost cells for derived mf 
Vector<int> HMM::tag_value_greater_flag;        // 1 if value_greater, else 0
Vector<int> HMM::tag_value_less_flag;           // 1 if value_less, else 0
Vector<int> HMM::tag_gradient_flag;             // 1 if gradient, else 0
Vector<amrex::Real> HMM::tag_value_greater;          // actual values
Vector<amrex::Real> HMM::tag_value_less;             // actual values
Vector<amrex::Real> HMM::tag_gradient;               // actual values
Vector<Vector<amrex::Real>> HMM::tag_domain_lo;               // actual values
Vector<Vector<amrex::Real>> HMM::tag_domain_hi;               // actual values




std::string    HMM::limiter = "minmod"; 
int      HMM::NUM_GROW        = 2;  // number of ghost cells used in Sborder

// Box mapping functions -- used by derived fields. 
static Box the_same_box (const Box& b) { return b; } 
static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Create arrays for boundary condition indexing.  
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//  BC labels: 
//  0: Interior, 
//  1: Inflow, 
//  2: Outflow,  
//  3: Symmetry,     
//  4: SlipWall,     
//  5: NoSlipWall
//  6: SlipWall, zeroGrad instead of reflect_even 


//      Interior, Inflow, Outflow,   Symmetry,      SlipWall,   NoSlipWall,     SlipWall zeroGrad
static int scalar_bc[] = 
    {
        INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, FOEXTRAP
    };

static int norm_vel_bc[] = 
    {
        INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD,  REFLECT_ODD,  REFLECT_ODD, REFLECT_ODD
    };

static int tang_vel_bc[] = 
    {
        INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_ODD, FOEXTRAP
    };

static
void
set_scalar_bc (BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++) 
    {    
      bc.setLo(i,scalar_bc[lo_bc[i]]);
      bc.setHi(i,scalar_bc[hi_bc[i]]);
    }    
}

static
void
set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,norm_vel_bc[lo_bc[0]]);
  bc.setHi(0,norm_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)
  bc.setLo(1,tang_vel_bc[lo_bc[1]]);
  bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
  bc.setLo(2,tang_vel_bc[lo_bc[2]]);
  bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,tang_vel_bc[lo_bc[0]]);
  bc.setHi(0,tang_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)
  bc.setLo(1,norm_vel_bc[lo_bc[1]]);
  bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
  bc.setLo(2,tang_vel_bc[lo_bc[2]]);
  bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,tang_vel_bc[lo_bc[0]]);
  bc.setHi(0,tang_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)
  bc.setLo(1,tang_vel_bc[lo_bc[1]]);
  bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
  bc.setLo(2,norm_vel_bc[lo_bc[2]]);
  bc.setHi(2,norm_vel_bc[hi_bc[2]]);
#endif
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Begin class function definitions.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//Default constructor.  Builds invalid object.
//
HMM::HMM ()
{
}

//
//The basic constructor.
//
HMM::HMM (Amr&            papa,
     	                  int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    BL_PROFILE("HMM::HMM()");
    if (verbose)
    {
        Print() << "~~~~ Creating a new HMM object for level " << level << " ~~~~" << std::endl; 
    } 

    // build metrics: 
    init_mfs(); 
     
    // ~~~~ Flux register initialization: 
    if (HMM::verbose) Print() << "initializing flux register on level " << level << std::endl;
    flux_reg = 0;
    if (level > 0 && do_reflux)
    {
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
        flux_reg->setVal(0.0);
    }    

    // // ~~~~ Initialize class multifabs: 
    // Sborder.define(grids,dmap,NUM_STATE,NUM_GROW);
    // source_reac.define(grids, dmap, NUM_STATE, 0);
    // source_conv.define(grids, dmap, NUM_STATE, 0);
    // for (int j = 0; j < BL_SPACEDIM; j++)
    // {
    //     BoxArray ba = get_new_data(State_Type).boxArray();
    //     ba.surroundingNodes(j);
    //     fluxes[j].define(ba, dmap, NUM_STATE, 0); 
    // }
    
    // ~~~ Some sanity checks: 
    if (NUM_GROW == 0)
    {
        Print() << "NUM_GROW is zero. Exiting." << std::endl;
        amrex::Error();
    }

    // ~~~~ Set all time level data to zero: 
    for (int k = 0; k < num_state_type; k++) 
    {
        MultiFab& data = get_new_data(k);
        data.setVal(0.0, data.nGrow());
    }    
    
    // Allocate space for old data only for State_Type! 
    if (state[State_Type].hasOldData()) 
    {
        if (HMM::verbose) Print() << "Space for old data already exists." << std::endl; 
    }    
    else
    {
        if (HMM::verbose) Print() << "Allocating space for old data within AmrLevel constructor... " << std::endl;
        state[State_Type].allocOldData();        // allocates space for old data
    }
    
    if (verbose)
    {
        Print() << "~~~~ Successfully created HMM object for level " << level << " ~~~~\n" << std::endl;
    }

}

//
//The destructor.
//
HMM::~HMM () 
{
   delete flux_reg;
}
// HMM::~HMM() = default; 



//
//Restart from a checkpoint file.
//
void
HMM::restart (Amr&          papa,
                      std::istream& is,
                      bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    // Initialize additional multifabs: 
    init_mfs(); 

    flux_reg = 0;
    if (level > 0 && do_reflux)
    {
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
        flux_reg->setVal(0.0);
    }

    if (HMM::verbose) Print() << "Done with restart on level " << level << std::endl;
}

void 
HMM::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old) 
{
}

//
//Configure the variables to be plotted
//
void
HMM::setPlotVariables ()
{
    
    // Set the defaults: 
    // --- this looks for the following in the input file:
    //      1) "amr.plot_vars" 
    //      2) "amr.derive_plot_vars"
    AmrLevel::setPlotVariables();

    // // Test this by removing all the state variables from the plot file: 
    // for (int i = 0; i < desc_lst[State_Type].nComp(); i++) 
    // {
    //     if (desc_lst[State_Type].name(i) == "rho") 
    //     {
    //         parent->deleteStatePlotVar(desc_lst[State_Type].name(i));
    //     }
    // }
}

//
//Write a plotfile to specified directory.
//
void
HMM::writePlotFile (const std::string& dir,
	 	            std::ostream&      os,
                    VisMF::How         how)
{
    
    // Base version: This overwrites the state variable output with a derived quantity. 
    // AmrLevel::writePlotFile (dir,os,how);


    // Castro version: this does not overwrite  
    plotFileOutput(dir, os, how, 0);
}

//
//Define data descriptors.
//
void
HMM::variableSetUp ()
{
    BL_PROFILE("HMM::variableSetUp()");
    BL_ASSERT(desc_lst.size() == 0);

    // ~~~~ Get parameters 
    read_params();

    // ~~~ Configure cantera: 
    // configure_cantera(); 


    // ~~~~ Add the state descriptors:
    desc_lst.addDescriptor( /*indx = */             State_Type,
                            /*IndexType = */        IndexType::TheCellType(),
                            /*StateDescriptor::TimeCenter = */ StateDescriptor::Point,
                            /*nextra (ghost?) = */  0,
                            /*num_comp = */         NUM_STATE,
			                /*Interpolater = */     &pc_interp);// pc_interp 


    // ~~~~ Set the boundary conditions:
    Vector<BCRec>       bcs_state(NUM_STATE);
    Vector<std::string> name_state(NUM_STATE);
   


    BCRec bc;




    // passive scalar: 
    set_scalar_bc(bc, phys_bc);
    bcs_state[C_PHI]      = bc; 
    name_state[C_PHI]     = "phi"; 


    // ~~~~ Declare the state boundary condition function. 
    if (HMM::verbose) Print() << "Declaring state boundary condition function... " << std::endl;
    StateDescriptor::BndryFunc stateBndryFunc(hmm_statefill_v2);


    // ~~~~ Set components: 
    desc_lst.setComponent(  /*indx = */ State_Type, 
                            /*comp = */ C_PHI, 
                            /*nm = */   name_state, 
                            /*bc = */   bcs_state, 
			                /*func = */ stateBndryFunc);



    num_state_type = desc_lst.size();

}

//
//Cleanup data descriptors at end of run.
//
void
HMM::variableCleanUp () 
{
    desc_lst.clear();
}

//
//Initialize grid data at problem start-up.
//
void
HMM::initData ()
{
    BL_PROFILE("HMM::initData()");

    const Real* dx  = geom.CellSize(); // this returns pointer. geom.CellSizeArray() returns GpuArray object.
    const Real* prob_lo = geom.ProbLo(); // geom.ProbLoArray() returns GpuArray object
    MultiFab& S_new = get_new_data(State_Type); // state[state_indx].newData() -- returns pointer to new time data
    Real cur_time   = state[State_Type].curTime();
    
    if (verbose) 
    {
        amrex::Print() << "~~~~ Initializing the data at level " << level << " ~~~~" << std::endl;
        amrex::Print() << "dx = " << dx[0] << std::endl; 
        amrex::Print() << "prob_lo = " << prob_lo[0] << std::endl; 
        amrex::Print() << "cur_time = " << cur_time << std::endl; 
    }

    

    // Set initial condition 
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        // Get the box and multifab data 
        const Box& box     = mfi.validbox();
        amrex::Array4<amrex::Real> const& data = S_new[mfi].array(); 

        #ifdef HMM_USE_DF
        amrex::Array4<amrex::Real> const& DF_data = DF_multifab[mfi].array(); 
        #endif
        
        // Get the geometry data: 
        GeometryData geomdata = geom.data();
        
        // Initialize problem data: 
        amrex::ParallelFor(box,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {    
            problem_initialize_state_data(i, j, k, 
                                data, 
                                #ifdef HMM_USE_DF
                                DF_data,
                                #endif
                                geomdata);
        });  
    }

    if (verbose) 
    {
        amrex::Print() << "~~~~ Done initializing the level " << level 
                       << " data ~~~~\n" << std::endl;
    }
}

//
//Initialize data on this level from another HMM (during regrid).
//
void
HMM::init (AmrLevel &old)
{
    BL_PROFILE("HMM::init(old)");
    
    HMM* oldlev = (HMM*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[State_Type].curTime();
    Real prev_time = oldlev->state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);


    // cycle through the different state types: 
    for (int s = 0; s < num_state_type; ++s) 
    {
        MultiFab& state_MF = get_new_data(s);
        
        // Fill the current time data: 
        FillPatch(old, state_MF, state_MF.nGrow(), cur_time, s, 0, state_MF.nComp());

        // Fill the previous time data, only if the old object had old data
        // --- this might be unnecessary 
        if (oldlev->state[s].hasOldData()) 
        {
            if (!state[s].hasOldData()) 
            {
                state[s].allocOldData();
            }
            MultiFab& old_state_MF = get_old_data(s);
            FillPatch(old, old_state_MF, old_state_MF.nGrow(), prev_time, s, 0, old_state_MF.nComp());
        }
    }
    
    if (HMM::verbose) Print() << "Done with HMM::init (AmrLevel &old) on level " << level << std::endl;

}

//
//Initialize data on this level after regridding if old level did not previously exist
//
void
HMM::init ()
{
    BL_PROFILE("HMM::init()");
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();
    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);
    setTimeLevel(cur_time,dt_old,dt);

    if (HMM::verbose) Print() << "Calling init FillCoarsePatch on level " << level << std::endl;
    for (int s = 0; s < num_state_type; ++s)
    {
        MultiFab& state_MF = get_new_data(s);
        FillCoarsePatch(state_MF, 0, cur_time, s, 0, state_MF.nComp(), state_MF.nGrow());
    }

    if (HMM::verbose) Print() << "Done with HMM::init() on level " << level << std::endl;
}

//
//Estimate time step.
//

Real
HMM::estTimeStep (Real)
{
/*    
    BL_PROFILE("HMM::estTimeStep(Real)");

    if (dt_fixed > 0)
    {
        
        // Estimate the CFL number: 
        MultiFab& S_new = get_new_data(State_Type);
        const auto dx = geom.CellSizeArray(); 
        amrex::Real cfl_est = 
        amrex::ReduceMax
        (
            S_new,
            0,   
            [=] AMREX_GPU_HOST_DEVICE
            (
                amrex::Box const& bx, 
                const amrex::Array4<const amrex::Real>& fab_arr
            ) noexcept -> 
                    amrex::Real{ return hmm_cfl_from_dt( bx, fab_arr, AMREX_D_DECL(dx[0], dx[1], dx[2])); }
        );  
        amrex::ParallelDescriptor::ReduceRealMax(cfl_est);
        
        if (HMM::verbose)
        {    
            Print() << "Using fixed dt of " << dt_fixed << "\t CFL = " << cfl_est << std::endl;
        }
        else
        {
            if (level == 0)
            {
                Print() << "Using fixed dt of " << dt_fixed << "\t CFL = " << cfl_est << std::endl;
            }

        }
        return dt_fixed;
    
    }
    else
    {
        // dt, time, dx 
        Real time = state[State_Type].curTime();
        const auto dx = geom.CellSizeArray(); 

        // get state multifab
        MultiFab& S_new = get_new_data(State_Type);
        if (HMM::verbose)
        {
            Print() << "Entered estTimeStep. Using cfl =  " << cfl << std::endl;
        }
        else
        {
            if (level == 0)
            {
                Print() << "Entered estTimeStep. Using cfl =  " << cfl << std::endl;
            }
        }

        // ~~~~ Hyperbolic timestep: 
        amrex::Real dt_est_hyp = 
        amrex::ReduceMin
        (
            S_new,
            0,   
            [=] AMREX_GPU_HOST_DEVICE
            (
                amrex::Box const& bx, 
                const amrex::Array4<const amrex::Real>& fab_arr
            ) noexcept -> 
                    amrex::Real{ return hmm_dt_hyperbolic( bx, fab_arr, AMREX_D_DECL(dx[0], dx[1], dx[2])); }
        );  
        dt_est_hyp = amrex::min<amrex::Real>(dt_max, dt_est_hyp); 

        // ~~~~ Viscosity timestep 

        // ~~~~ Thermal conductivity timestep

        // ~~~~ Species diffusion timestep 

        // // ~~~~ Chemical timestep: 
        // amrex::Real dt_est_chem;
        // if (use_chemical_timestep == 1)
        // {
        //     dt_est_chem = 
        //     amrex::ReduceMin
        //     (
        //         S_new,
        //         0,   
        //         [=] AMREX_GPU_HOST_DEVICE
        //         (
        //             amrex::Box const& bx, 
        //             const amrex::Array4<const amrex::Real>& fab_arr
        //         ) noexcept -> 
        //                 amrex::Real{ return hmm_dt_chemical( bx, fab_arr, AMREX_D_DECL(dx[0], dx[1], dx[2])); }
        //     );  
        //     dt_est_chem = amrex::min<amrex::Real>(dt_max, dt_est_chem); 
        // } 
        // else
        // {
        //     dt_est_chem = std::numeric_limits<amrex::Real>::max();
        // }
        // 

        // // ~~~~ Get min of hyperbolic dt and chemical dt:  
        // dt_est_hyp = amrex::min<amrex::Real>(dt_est_hyp, dt_est_chem);

        // ~~~~ Reduction across MPI ranks 
        amrex::ParallelDescriptor::ReduceRealMin(dt_est_hyp);

        if (HMM::verbose) 
        {
            amrex::Print() << "Done with HMM::estTimeStep at level " << level 
                           << ":  dt_est = " << dt_est_hyp << std::endl;
            amrex::Print() << "\n";
        }
        else
        {
            if (level == 0)
            {
                amrex::Print() << "Done with HMM::estTimeStep at level " << level 
                           << ":  dt_est = " << dt_est_hyp << std::endl;
                amrex::Print() << "\n";
            }
        }

        return dt_est_hyp;
    }
*/    
    return 0.0;
}

//
//Compute initial time step.
//
Real
HMM::initialTimeStep ()
{
    return estTimeStep(0.0);
}

//
//Compute initial `dt'.
//
void
HMM::computeInitialDt (int                   finest_level,
                               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
    BL_PROFILE("HMM::computeInitialDt()");
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    if (HMM::verbose) Print() << "~~~~ Entering ComputeInitialDt() ~~~~" << std::endl;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    //n_factor = 1;
    //for (int i = 0; i <= finest_level; i++)
    //{
    //    n_factor *= n_cycle[i];
    //    dt_level[i] = dt_0/n_factor;
    //}


    if (do_subcycle)
    {
        if (HMM::verbose) Print() << "Do subcycle: dt_coarse = dt_fine * rf" << std::endl; 
        n_factor = 1;
        for (int i = 0; i <= finest_level; i++)
        {   
            n_factor *= n_cycle[i];
            dt_level[i] = dt_0/n_factor;
        }   
    }
    else
    {
        if (HMM::verbose) Print() << "No subcycle: dt_coarse = dt_fine" << std::endl; 
        n_factor = 1;
        for (int i = 0; i <= finest_level; i++)
        {   
            n_factor *= n_cycle[i];
            dt_level[i] = dt_0/n_factor;
        }
        
        Real dt_finest_level = dt_level[finest_level];  
        for (int i = 0; i <= finest_level; i++)
        {
            dt_level[i] = dt_finest_level; 
            n_cycle[i] = 1;
        }
    }




    if (HMM::verbose) Print() << "~~~~ Done with ComputeInitialDt() ~~~~\n" << std::endl;
}

//
//Compute new `dt'.
//
void
HMM::computeNewDt (int                   finest_level,
                           int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{

    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++)
    {   
        HMM& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }   

    if (post_regrid_flag == 1)  
    {   
        //  
        // Limit dt's by pre-regrid dt
        //  
        for (int i = 0; i <= finest_level; i++)
        {   
            dt_min[i] = std::min(dt_min[i],dt_level[i]);
        }   
    }   
    else 
    {   
        //  
        // Limit dt's by change_max * old dt
        //  
        static Real change_max = 1.1;
        for (int i = 0; i <= finest_level; i++)
        {   
            dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
        }   
    }   
    
    //  
    // Find the minimum over all levels
    //  
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {   
        n_factor *= n_cycle[i]; // n_cycle[i] is the number of subcycled timesteps at level i. 
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }   

    //  
    // Limit dt's by the value of stop_time.
    //  
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }   

    // 
    // Re-compute dt at each level to conform to the subcycle criterion
    //
    if (do_subcycle)
    {
        if (HMM::verbose) Print() << "Do subcycle: dt_coarse = dt_fine * rf" << std::endl; 
        n_factor = 1;
        for (int i = 0; i <= finest_level; i++)
        {   
            n_factor *= n_cycle[i];
            dt_level[i] = dt_0/n_factor;
        }   
    }
    else
    {
        if (HMM::verbose) Print() << "No subcycle: dt_coarse = dt_fine" << std::endl; 
        n_factor = 1;
        for (int i = 0; i <= finest_level; i++)
        {   
            n_factor *= n_cycle[i];
            dt_level[i] = dt_0/n_factor;
        }
        
        Real dt_finest_level = dt_level[finest_level];  
        for (int i = 0; i <= finest_level; i++)
        {
            dt_level[i] = dt_finest_level; 
            n_cycle[i] = 1;
        }
    }
    



}



//
//Do work after regrid().
//At this point, we have data on this level that could have been interpolated from a coarser level.
//Should check to make sure that this data satisfies physical properties. 
//
void
HMM::post_regrid (int lbase, int new_finest) 
{

    // synchronize grids
    if (level > 0)
    // Only perform this on the coarse level.       
        return;
    else
    {

        if (HMM::verbose) Print() << "Doing post-initialization step for synchronizing grids." << std::endl;
        //
        // Average data down from finer levels
        // so that conserved data is consistent between levels.
        //
        int finest_level = parent->finestLevel();
        for (int k = finest_level-1; k>= 0; k--)
        {
            getLevel(k).avgDown();
        }

    }

     

    if (HMM::verbose) amrex::Print() << "~~~~ Successfully completed HMM::post_regrid for level " << level << " ~~~~" << std::endl; 
}

//
//Do work after a restart().
//
void
HMM::post_restart() 
{
    Print() << "Nothing to do after restart..." << std::endl;
}

amrex::Real HMM::advance (amrex::Real time,
                                 amrex::Real dt,
                                 int  iteration,
                                 int  ncycle)
{
    return dt;
} 

amrex::Real HMM::advance_react (amrex::Real time,
                                 amrex::Real dt)
{
    return dt;
}    


//
//Do work after berger-collela timestep.
//
void
HMM::post_timestep (int iteration)
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
}

void
HMM::post_timestep_react (int iteration)
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
}

//
//Do work after init().
//
void
HMM::post_init (Real stop_time)
{   
    // Only perform this on the coarse level. 
    if (level > 0)
    {
        return;
    }
    else
    {
        if (HMM::verbose) Print() << "Doing post-initialization step for synchronizing grids." << std::endl;
        //
        // Average data down from finer levels
        // so that conserved data is consistent between levels.
        //
        int finest_level = parent->finestLevel();
        for (int k = finest_level-1; k>= 0; k--)
        {
            getLevel(k).avgDown(State_Type);
        }
                        
    }
}

//
//Error estimation for regridding.
//
void
HMM::errorEst (TagBoxArray& tags,
	                    int          clearval,
                        int          tagval,
                        Real         time,
                        int          n_error_buf,
                        int          ngrow)
{
    BL_PROFILE("HMM::errorEst()"); 
    if (verbose) 
    {
        amrex::Print() << "~~~~ Entering HMM::errorEst for level " << level << " ~~~~" << std::endl; 
    }
    
    // Get the grown state data: 
    int tag_ng = 1;
    amrex::MultiFab S_grown( /*const BoxArray& = */ get_new_data(State_Type).boxArray(),
                            /*const DistributionMapping& dm = */ get_new_data(State_Type).DistributionMap(), 
                            /*int ncomp = */ NUM_STATE, 
                            /*int ngrow = */ tag_ng);

    if (HMM::verbose) Print() << "Fill patch in errorEst... " << std::endl;
    FillPatch(  /* AmrLevel& amrlevel = */ *this, 
            /* MultiFab& leveldata = */ S_grown, 
            /* int boxGrow = */ S_grown.nGrow(), 
            /* Real time = */ state[State_Type].curTime(),
            /* int index = */ State_Type, 
            /* int scomp = */ C_PHI,  
            /* int ncomp = */ NUM_STATE, 
            /* int dcomp = */ 0);

    // Declare temporary fab for derived data: 
    FArrayBox derived_data;

    // dummy bcrec required for derived functions: 
    const int *bcrec; 

    // Begin tagging routine: 
    for (MFIter mfi(tags); mfi.isValid(); ++mfi)
    {
        // Get boxes: 
        const Box& box_valid = mfi.validbox();
        const Box& box_grown = amrex::grow(box_valid, tag_ng);

        // Get array4 for tag: 
        auto tag = tags[mfi].array(); 

        // Get FAB for grown state data: 
        FArrayBox& state_data = S_grown[mfi]; 
        
        // Get FAB for derived data: 
        derived_data.resize(box_grown, 1);
       
        // zero the tags box (check to see if this is done by default -- may not be necessary)  
        // Print() << "Zero-ing tag box " << std::endl;
        amrex::ParallelFor(box_valid, 
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {   
            tag(i,j,k) = TagBox::CLEAR;
        }); 

        // Run the tagging routine: this loops through all the tagging vars specified in input file. 
        //#include "HMM_tagging.H"
    } // end MFIter

    if (verbose) 
    {
        amrex::Print() << "~~~~ Done with HMM::errorEst for level " << level << " ~~~~" << std::endl; 
    }
}

void
HMM::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("hmm");   

 
    pp.query("cfl",cfl);
    pp.query("dt_max",dt_max); 

    pp.query("do_reflux",do_reflux);
    pp.query("do_subcycle", do_subcycle); 
    pp.query("dt_fixed", dt_fixed);



    const Geometry& dgeom = DefaultGeometry(); 
    // This tutorial code only supports Cartesian coordinates.
    if (! dgeom.IsCartesian()) 
    {
        amrex::Abort("Please set geom.coord_sys = 0");
    }

    // Get boundary conditions:
    Vector<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    ParmParse pp_amr("amr");
    pp_amr.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
    pp_amr.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++) 
    {    
        // Set phys_bc, which is a BCRec object.
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }    

    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (dgeom.isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir<BL_SPACEDIM; dir++)
        {
            if (dgeom.isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "read_params:periodic in direction "
                              << dir
                              << " but low BC is not Interior\n";
                    amrex::Error();
                }    
                if (hi_bc[dir] != Interior)
                {    
                    std::cerr << "read_params:periodic in direction "
                              << dir
                              << " but high BC is not Interior\n";
                    amrex::Error();
                }
            }
        }
    }
    else
    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir=0; dir<BL_SPACEDIM; dir++)
        {
            if (lo_bc[dir] == Interior)
            {
                std::cerr << "read_params:interior bc in direction "
                          << dir
                          << " but not periodic\n";
                amrex::Error();
            }
            if (hi_bc[dir] == Interior)
            {
                std::cerr << "read_params:interior bc in direction "
                          << dir
                          << " but not periodic\n";
                amrex::Error();
            }
        }
    }

    // read the list of Euler fluxes to use on each level 
    int max_level; 
    pp_amr.query("max_level", max_level); 


    Print() << "MAXIMUM LEVEL IS " << max_level << std::endl;

    

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // read the field tagging params. 
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {   
        // ~~~~ Get the pp obj 
        ParmParse pptag("tagging"); 

        // ~~~~ list containing names of refinement indicators as supplied in input file 
        Vector<std::string> refinement_indicators; 
        pptag.queryarr("refinement_indicators", refinement_indicators, 0, pptag.countval("refinement_indicators"));

        // ~~~~ below lists are of size refinement_indicators.size()  
        // Vector<std::string> tag_field_list;                // field to derive from  
        // Vector<int> tag_max_level_list;            // max level to tag 
        // Vector<int> tag_value_greater_flag;        // 1 if value_greater, else 0
        // Vector<int> tag_value_less_flag;           // 1 if value_less, else 0 
        // Vector<int> tag_gradient_flag;             // 1 if gradient, else 0 
        // Vector<amrex::Real> tag_value_greater;          // actual values 
        // Vector<amrex::Real> tag_value_less;             // actual values
        // Vector<amrex::Real> tag_gradient;               // actual values 
        // Vector<amrex::Real> tag_domain_lo;           // domain boundary for fixed domain tagging 
        // Vector<amrex::Real> tag_domain_hi; 

        // ~~~~ Extract tagging parameters from input file 
        if (HMM::verbose) Print() << "\n beginning tag params loop" << std::endl;
        Vector<amrex::Real> value_vec; 
        value_vec.resize(3);
        for (int i = 0; i < refinement_indicators.size(); i++)
        {
            
            // Set the tagging param: 
            std::string ref_prefix = "tagging." + refinement_indicators[i]; 
            Print() << "tagging param: " << ref_prefix << std::endl; 
    
            // get pp object
            Print() << "get pp object." << std::endl;
            ParmParse ppr(ref_prefix); 

            // get the field:
            Print() << "get the field_name. " << std::endl;
            std::string field;
            ppr.get("field_name", field);
            tag_field_list.push_back(field);
            
            // get the max level: 
            Print() << "get the max level. " << std::endl; 
            int max_level; 
            if (ppr.countval("max_level") > 0)
            {
                ppr.get("max_level", max_level); 
            }
            else 
            {
                max_level = 0; 
            }
            tag_max_level_list.push_back(max_level);


            // get the tagging type:
            Print() << "get the tagging types. " << std::endl; 

            // value greater: 
            Print() << "\tvalue_greater... " << std::endl;
            int flag;  
            amrex::Real value;
            if (ppr.countval("value_greater") > 0) 
            {
                ppr.get("value_greater", value);
                flag = 1; 
            }
            else 
            {
                value = 1.0e30; 
                flag = 0; 
            }
            tag_value_greater_flag.push_back(flag); 
            tag_value_greater.push_back(value);

            // value less: 
            Print() << "\tvalue_less... " << std::endl;
            if (ppr.countval("value_less") > 0) 
            {
                ppr.get("value_less", value);
                flag = 1; 
            }
            else 
            {
                value = -1.0e30; 
                flag = 0; 
            }
            tag_value_less_flag.push_back(flag); 
            tag_value_less.push_back(value);

            // gradient 
            Print() << "\tgradient... " << std::endl;
            if (ppr.countval("gradient") > 0) 
            {
                ppr.get("gradient", value);
                flag = 1; 
            }
            else 
            {
                value = 1.0e30; 
                flag = 0; 
            }
            tag_gradient_flag.push_back(flag); 
            tag_gradient.push_back(value);


            // fixed domain: 
            Print() << "\tfixed domain..." << std::endl;
            if (ppr.countval("tag_domain_lo") > 0)
            {
                ppr.getarr("tag_domain_lo", value_vec, 0, 3);
            }
            else
            {
                value_vec[0] = -1.0e30;
                value_vec[1] = -1.0e30;
                value_vec[2] = -1.0e30;
            }
            tag_domain_lo.push_back(value_vec);

            if (ppr.countval("tag_domain_hi") > 0)
            {
                ppr.getarr("tag_domain_hi", value_vec, 0, 3);
            }
            else
            {
                value_vec[0] = -1.0e30;
                value_vec[1] = -1.0e30;
                value_vec[2] = -1.0e30;
            }
            tag_domain_hi.push_back(value_vec);

            // Summary: 
            Print() << "Tagging summary for " << refinement_indicators[i] << std::endl;
            Print() << "\ttag_field_list[" << i << "] = " << tag_field_list[i] << std::endl;  
            Print() << "\ttag_max_level_list[" << i << "] = " << tag_max_level_list[i] << std::endl;  
            Print() << "\ttag_value_greater_flag[" << i << "] = " << tag_value_greater_flag[i] << std::endl; 
            Print() << "\ttag_value_greater[" << i << "] = " << tag_value_greater[i] << std::endl; 
            Print() << "\ttag_value_less_flag[" << i << "] = " << tag_value_less_flag[i] << std::endl; 
            Print() << "\ttag_value_less[" << i << "] = " << tag_value_less[i] << std::endl; 
            Print() << "\ttag_gradient_flag[" << i << "] = " << tag_gradient_flag[i] << std::endl; 
            Print() << "\ttag_gradient[" << i << "] = " << tag_gradient[i] << std::endl; 
            Print() << "\ttag_domain_lo[" << i << "][0] = " << tag_domain_lo[i][0] << std::endl;
            Print() << "\ttag_domain_lo[" << i << "][1] = " << tag_domain_lo[i][1] << std::endl;
            Print() << "\ttag_domain_lo[" << i << "][2] = " << tag_domain_lo[i][2] << std::endl;
            Print() << "\ttag_domain_hi[" << i << "][0] = " << tag_domain_hi[i][0] << std::endl;
            Print() << "\ttag_domain_hi[" << i << "][1] = " << tag_domain_hi[i][1] << std::endl;
            Print() << "\ttag_domain_hi[" << i << "][2] = " << tag_domain_hi[i][2] << std::endl;

            Print() << "\n" << std::endl;
        }
    }

#ifdef AMREX_PARTICLES
    Print() << "Particles not implemented. Recompile without particles. Exiting... " << "\n"; 
    amrex::Error();
#endif 

}

// void
// HMM::configure_cantera ()
// {

// }

void
HMM::reflux ()
{

    BL_ASSERT(level < parent->finestLevel());

}

void
HMM::avgDown ()
{
    // Do not do for finest level.
    if (level == parent->finestLevel()) return;
    
    
    const auto strt = amrex::second(); 
    
    // Only do for state_type
    avgDown(State_Type);

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        auto      end    = amrex::second() - strt;

        ParallelDescriptor::ReduceRealMax(end,IOProc);

        amrex::Print() << "HMM::avgDown() at level " << level
                       << " : time = " << end << std::endl;
    }

}

void
HMM::avgDown (int state_indx)
{
    // Do not do for finest level. 
    if (level == parent->finestLevel()) return;

    // Get fine and coarse multifabs: 
    HMM& fine_lev = getLevel(level+1); // get fine AmrLevel object  
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx); // Get fine multifab
    MultiFab&  S_crse   = get_new_data(state_indx); // Get coarse (current) multifab 
   
    // Do the volume average on the coarse (current) multifab 
    amrex::average_down(S_fine,
                        S_crse,
                        fine_lev.geom,
                        geom,
                        0,
                        S_fine.nComp(),
                        parent->refRatio(level));

}



void HMM::plotFileOutput(const std::string& dir, 
                       std::ostream& os,
                       VisMF::How how, 
                       const int is_small)
{
    //   
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //   
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++) 
    {
        for (int comp = 0; comp < desc_lst[typ].nComp(); comp++) 
        {
            if (( (parent->isStatePlotVar(desc_lst[typ].name(comp)) && is_small == 0) ||
                    (parent->isStateSmallPlotVar(desc_lst[typ].name(comp)) && is_small == 1) ) &&
                    desc_lst[typ].getType() == IndexType::TheCellType()   ) 
            {
                plot_var_map.push_back(std::pair<int,int>(typ,comp));
            }
        }
    }

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (auto it = dlist.begin(); it != dlist.end(); ++it)
    {
        if ((parent->isDerivePlotVar(it->name()) && is_small == 0) ||
            (parent->isDeriveSmallPlotVar(it->name()) && is_small == 1))
        {
               derive_names.push_back(it->name());
               num_derive = num_derive + it->numDerive();
        }
    }

    int n_data_items = plot_var_map.size() + num_derive;

    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0) {
          amrex::Error("Must specify at least one valid data item to plot");
        }

        os << n_data_items << '\n';

        //
        // Names of variables -- first state, then derived
        //
        for (int i =0; i < plot_var_map.size(); i++)
        {
            int typ = plot_var_map[i].first;
            int comp = plot_var_map[i].second;
            os << desc_lst[typ].name(comp) << '\n';
        }

        for (auto it = derive_names.begin(); it != derive_names.end(); ++it)
        {
            const DeriveRec* rec = derive_lst.get(*it);
            if (rec->numDerive() > 1) {
                for (int i = 0; i < rec->numDerive(); ++i) {
                    os << rec->variableName(0) + '_' + std::to_string(i) + '\n';
                }
            }
            else {
                os << rec->variableName(0) << '\n';
            }
        }


        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++) {
            os << geom.ProbLo(i) << ' ';
        }
        os << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++) {
            os << geom.ProbHi(i) << ' ';
        }
        os << '\n';
        for (int i = 0; i < f_lev; i++) {
          os << parent->refRatio(i)[0] << ' ';
        }
        os << '\n';
        for (int i = 0; i <= f_lev; i++) {
          os << parent->Geom(i).Domain() << ' ';
        }
        os << '\n';
        for (int i = 0; i <= f_lev; i++) {
          os << parent->levelSteps(i) << ' ';
        }
        os << '\n';
        for (int i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++) {
              os << parent->Geom(i).CellSize()[k] << ' ';
            }
            os << '\n';
        }
        os << (int) geom.Coord() << '\n';
        os << "0\n"; // Write bndry data.

    }

    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/') {
      FullPath += '/';
    }
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor()) {
        if (!amrex::UtilCreateDirectory(FullPath, 0755)) {
            amrex::CreateDirectoryFailed(FullPath);
        }
    }
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (int i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (int n = 0; n < BL_SPACEDIM; n++) {
              os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
            }
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    const int nGrow = 0;
    MultiFab  plotMF(grids,dmap,n_data_items,nGrow);
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (int i = 0; i < plot_var_map.size(); i++)
    {
        int typ  = plot_var_map[i].first;
        int comp = plot_var_map[i].second;
        this_dat = &state[typ].newData();
        MultiFab::Copy(plotMF,*this_dat,comp,cnt,1,nGrow);
        cnt++;
    }
    //
    // Cull data from derived variables.
    //
    if (dlist.size() > 0)
    {
        for (auto it = dlist.begin(); it != dlist.end(); ++it)
        {
            if ((parent->isDerivePlotVar(it->name()) && is_small == 0) ||
                (parent->isDeriveSmallPlotVar(it->name()) && is_small == 1)) {

                auto derive_dat = derive(it->variableName(0), cur_time, nGrow);
                MultiFab::Copy(plotMF, *derive_dat, 0, cnt, it->numDerive(), nGrow);
                cnt = cnt + it->numDerive();
            }
        }
    }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;

    const Real io_start_time = ParallelDescriptor::second();

    if (amrex::AsyncOut::UseAsyncOut()) {
        VisMF::AsyncWrite(std::move(plotMF),TheFullPath);
    } else {
        VisMF::Write(plotMF,TheFullPath,how,true);
    }

    const Real io_time = ParallelDescriptor::second() - io_start_time;
    if (level == 0 && ParallelDescriptor::IOProcessor()) {
        writeJobInfo(dir, io_time);
    }

}




void
HMM::writeJobInfo(const std::string& dir, const Real io_time)
{
    // job_info file with details about the run
    std::ofstream jobInfoFile;
    std::string FullPathJobInfoFile = dir; 
    FullPathJobInfoFile += "/job_info";
    jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

    std::string PrettyLine = std::string(78, '=') + "\n";
    std::string OtherLine = std::string(78, '-') + "\n";
    std::string SkipSpace = std::string(8, ' ');

    // job information
    jobInfoFile << PrettyLine;
    jobInfoFile << " Job Information\n";
    jobInfoFile << PrettyLine;

    jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
    jobInfoFile << "\n";

    time_t now = time(0);

    // Convert now to tm struct for local timezone
    tm* localtm = localtime(&now);
    jobInfoFile   << "output date / time: " << asctime(localtm);

    char currentDir[FILENAME_MAX];
    if (getcwd(currentDir, FILENAME_MAX)) {
        jobInfoFile << "output dir:         " << currentDir << "\n";
    }    

    jobInfoFile << "I/O time (s):       " << io_time << "\n";

    jobInfoFile << "\n\n";

    jobInfoFile.close();
}




void 
HMM::init_mfs()
{

    if (HMM::verbose) Print() << "Entering init_mfs() function on level " << level << std::endl;

    // Volume 
    if (HMM::verbose) Print() << "\tDefining volume multifab..." << std::endl;
    volume.clear();
    volume.define(grids,dmap,1,NUM_GROW);
    geom.GetVolume(volume);


    if (HMM::verbose) Print() << "Done with init_mfs() function on level " << level << std::endl;
}
