#ifndef _HMMAmr_H_
#define _HMMAmr_H_
#include <AMReX_Config.H>

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BCRec.H>
#include <AMReX_AmrCore.H>
#include <AMReX_Amr.H> 

#include <iosfwd>
#include <list>
#include <memory>

namespace amrex {

// class AmrLevel;
// class LevelBld;
// class BoxDomain;
// #if defined(BL_USE_SENSEI_INSITU)
// class AmrInSituBridge;
// #endif


/**
HMMAmr is a derived class for Amr. 
Main added functionality:
    - redesign of the coarse timestepping routine to allow for Strang 
    - integration of in-situ neural network training  
**/


class HMMAmr 
    : public Amr 
{

public:
    //! The constructor.
    HMMAmr (LevelBld* a_levelbld /* One must pass LevelBld* as an argument now*/);

    HMMAmr (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord,
        LevelBld* a_levelbld /* One must pass LevelBld* as an argument now*/);


    //! The destructor.
    virtual ~HMMAmr ();


    // Define a new coarseTimeStep 
    void coarseTimeStepMod (Real stop_time);

    void RegridOnlyMod( Real time, bool do_io = true);
    
protected:

    // Run some pre-timeStep routines on level L
    void timeStep_init (int  level,
                       Real time,
                       Real stop_time);

    // Run the Berger-Collela advance (multigrid cycle): 
    void timeStep_BC (int  level,
                       Real time,
                       int  iteration,
                       int  niter);


    // Run the reaction multigrid cycle: 
    void timeStep_react (int  level,
                       Real time);


};



}


#endif /*_Amr_H_*/
