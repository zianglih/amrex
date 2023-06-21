#ifndef _HMM_H_
#define _HMM_H_

// C++ standard:
#include <memory>
#include <iostream>
#include <random>
#include <unistd.h>
// OpenMP:
#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

// AMREX: 
#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_PlotFileUtil.H> 
// #include <AMReX.H>

using namespace amrex;

// Main headers:
#include <BC_funcs.H> 

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// BOOK KEEPING 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// StateType enum. Will be useful when extending to different physical models. 
enum StateType  { 
                    State_Type = 0
                    #ifdef HMM_USE_POSTPROC
                    ,Postproc_Type 
                    #endif
                    #ifdef HMM_USE_DF
                    ,DF_Type
                    #endif 
                    #ifdef HMM_USE_KMEANS
                    ,Kmeans_Type
                    #endif
                };

// These index the "State_Type" and "flux" multifabs. 

constexpr int C_PHI     = 0;        // passive scalar
constexpr int NUM_STATE = 1;        // number of state variables.  

// These index the "primitive" multifab. 
constexpr int P_PHI     = 0;
constexpr int NUM_STATE_PRIM = 1;   // number of primitive variables

// These index the pos/neg "face" multifabs. 
constexpr int F_PHI     = 0; 
constexpr int NUM_STATE_FACE = 1; // number of face vars 


// Useful index values for dimensionality 
#if AMREX_SPACEDIM == 1
constexpr int dg0 = 1; 
constexpr int dg1 = 0; 
constexpr int dg2 = 0; 

#elif AMREX_SPACEDIM == 2
constexpr int dg0 = 1; 
constexpr int dg1 = 1; 
constexpr int dg2 = 0; 

#else
constexpr int dg0 = 1; 
constexpr int dg1 = 1; 
constexpr int dg2 = 1; 

#endif


// Small number for avoiding divide by zeros
constexpr amrex::Real SMALL = 1.0e-16; 

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Class for implementing HMM ideas 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class HMM
    :
    public amrex::AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    HMM ();
    //
    //The basic constructor.
    //
    HMM (amrex::Amr&     papa,
                 int             lev,
                 const amrex::Geometry& level_geom,
                 const amrex::BoxArray& bl,
                 const amrex::DistributionMapping& dm,
                 amrex::Real            time);
    //
    //The destructor.
    //
    virtual ~HMM () override;
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (amrex::Amr&   papa,
                          std::istream& is,
			  bool          bReadSpecial = false) override;

    virtual void checkPoint (const std::string& dir,
			     std::ostream&      os,
			     amrex::VisMF::How  how = amrex::VisMF::NFiles,
			     bool               dump_old = true) override;

    //
    //Configure the variables that are generated in the plotfile
    //
    void setPlotVariables () override; 


    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                amrex::VisMF::How  how) override;
   



///
    void plotFileOutput(const std::string& dir,
                        std::ostream& os,
                        amrex::VisMF::How how,
                        const int is_small);


    /// 
    /// From castro: Write job info to file
    /// 
    /// @param dir       directory to write to
    /// @param io_time   how long it took to generate the plotfile/checkpoint
    /// 
    void writeJobInfo (const std::string& dir, const amrex::Real io_time);



    //Define data descriptors.
    //
    static void variableSetUp ();
    //
    //Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    //Initialize grid data at problem start-up.
    //
    virtual void initData () override;
    //
    //Initialize data on this level from another HMM (during regrid).
    //
    virtual void init (amrex::AmrLevel& old) override;
    //
    //Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () override;
    //
    //Initialize extra multifabs and build metrics 
    //
    void init_mfs(); 

    // GENERAL MULTIFAB FUNCTIONS 
    //
    int numberOfElementsMF(amrex::MultiFab& mf);

    void sampleRandomMF(
        amrex::MultiFab& mf, 
        int N_ghost,
        int seed);

    void standardizeMF(
        amrex::MultiFab& mf,
        int idx_comp, 
        int N_ghost, 
        int N_elements);

    void getMeanStdMF(
        amrex::Real *mean_std,
        const amrex::MultiFab& mf, 
        int idx_comp, 
        int N_ghost,
        int N_elements);
   
    void normalizeMF(
        amrex::MultiFab& mf,
        amrex::Real shift,
        amrex::Real scale,
        int idx_comp,
        int N_ghost);
      


    //
    // Do an integration step on this level.
    // -- this advance function is configured in the context of Berger-Collela subcycling. 
    // -- if reactions are included here, the reactions are treated in the same multigrid cycle as convection
    //
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int  iteration,
                                 int  ncycle) override;



    //
    // Do a reaction advance on this level.
    // -- this is called when we use Strang splitting with independent V-cycles for the reaction step. 
    // -- Used in the HMMAmr class, which is derived from Amr. 
    //
    virtual amrex::Real advance_react (Real time,
                                          Real dt) override;



    void swap_data_and_time(amrex::Real dt); 
    void swap_data(); 

    //
    //Estimate time step.
    //
    amrex::Real estTimeStep (amrex::Real dt_old);
    //
    //Compute initial time step.
    //
    amrex::Real initialTimeStep ();
    //
    //Compute initial `dt'.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   amrex::Vector<int>&           n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                  stop_time) override;
    //
    //Compute new `dt'.
    //
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               amrex::Vector<int>&           n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                  stop_time,
                               int                   post_regrid_flag) override;
    //
    //Do work after advance.
    //
    virtual void post_timestep (int iteration) override;

    //
    //Do work after advance_react.
    //
    virtual void post_timestep_react (int iteration) override;

    //
    //Do work after regrid().
    //
    virtual void post_regrid (int lbase, int new_finest) override;
    //
    //Do work after a restart().
    //
    virtual void post_restart () override;
    //
    //Do work after init().
    //
    virtual void post_init (amrex::Real stop_time) override;
    //
    //Error estimation for regridding.
    //
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           amrex::Real         time,
			                int          n_error_buf = 0, int ngrow = 0) override;


  


    static int numGrow();



    //
    // The public data.
    //
    // static data members: 
    static int  NUM_GROW;
    
    static std::string limiter; 

    
    static amrex::Real  cfl;
    static amrex::Real  dt_fixed; 

    // chemistry/thermo: 
    

protected:

    static void read_params ();

    // inline functions implemented below
    HMM& getLevel (int lev);
    amrex::FluxRegister& getFluxReg ();
    amrex::FluxRegister& getFluxReg (int lev);
    amrex::FluxRegister& getFluxReg_y ();
    amrex::FluxRegister& getFluxReg_y (int lev);
    // amrex::MultiFab& getDFVelocity(); 

    void reflux ();

    void avgDown ();

    void avgDown (int state_indx);

    //
    // The protected data.
    //

    // Flux register object. 
    amrex::FluxRegister        *flux_reg;
    
    // Vector of multifabs for faces: this seems redundant for 1d, but will generalize to higher dimensions.  
    // amrex::Vector<std::unique_ptr<amrex::MultiFab> > fluxes;
    // amrex::Vector<std::unique_ptr<amrex::MultiFab> > Spos; 
    // amrex::Vector<std::unique_ptr<amrex::MultiFab> > Sneg;  
    // amrex::MultiFab Spos[BL_SPACEDIM];
    // amrex::MultiFab Sneg[BL_SPACEDIM]; 
    amrex::MultiFab             volume; // volume of cell 
    amrex::Array<amrex::MultiFab,3> area; // area of face
    amrex::MultiFab             Sborder; // Multifab for state, but with ghost cells.  
    amrex::MultiFab             Sprim; // Multifab for primitive state, with ghost cells. 
    amrex::MultiFab fluxes[BL_SPACEDIM];
    amrex::MultiFab             grad_scalar; // generic multifab for scalar gradient. 



    //
    // Static data members. 
    //
    static int          verbose;
    static int          use_viscosity; 
    static amrex::Real  dt_max; 
    static int          do_reflux;
    static int          do_subcycle; 
    static int num_state_type;    
  

 
    
    // variable names: 
    static std::vector<std::string> C_NAMES; // name id

    // for error estimation / tagging  
    static std::vector<std::string> err_list_names; // name id
    static std::vector<int> err_list_ng; // number of ghost cells 

    static Vector<std::string> tag_field_list; 
    static Vector<int> tag_max_level_list;            // max level to tag
    static Vector<int> tag_value_greater_flag;        // 1 if value_greater, else 0
    static Vector<int> tag_value_less_flag;           // 1 if value_less, else 0
    static Vector<int> tag_gradient_flag;             // 1 if gradient, else 0
    static Vector<int> tag_ng_list;                   // number of ghost cells for derived multifab 
    static Vector<amrex::Real> tag_value_greater;          // actual values
    static Vector<amrex::Real> tag_value_less;             // actual values
    static Vector<amrex::Real> tag_gradient;               // actual values
    static Vector<Vector<amrex::Real>> tag_domain_lo;               // actual values
    static Vector<Vector<amrex::Real>> tag_domain_hi;               // actual values
    
    // for boundary conditions: 
    static amrex::BCRec phys_bc; // this is read-in from the input file. 

};    

//
// Inlines.
//

// Returns the HMM object from any level. "parent" is the AMR class object.  
inline
HMM& HMM::getLevel (int lev)
{
    return *(HMM *) &parent->getLevel(lev);
}

// Returns the flux register object at the current level. 
inline
amrex::FluxRegister& HMM::getFluxReg ()
{
    //BL_ASSERT(flux_reg);
    return *flux_reg;
}

// Returns the flux register object at another level. 
inline
amrex::FluxRegister& HMM::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}

// returns the number of ghost cells used on the state. 
inline
int HMM::numGrow()
{
    return NUM_GROW;
}

// // returns the DF_velocity multifab 
// inline
// amrex::MultiFab& HMM::getDFVelocity()
// {
//     return DF_velocity; 
// }



#endif /*_HMM_H_ header guard*/
