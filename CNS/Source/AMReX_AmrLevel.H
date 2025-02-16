
#ifndef AMREX_AmrLevel_H_
#define AMREX_AmrLevel_H_
#include <AMReX_Config.H>

#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_LayoutData.H>
#include <AMReX_Derive.H>
#include <AMReX_BCRec.H>
#include <AMReX_Amr.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_StateDescriptor.H>
#include <AMReX_StateData.H>
#include <AMReX_VisMF.H>
#include <AMReX_RungeKutta.H>
#include <AMReX_FillPatcher.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBSupport.H>
#endif

#include <memory>
#include <map>

namespace amrex {

class TagBox;
class TagBoxArray;

/**
* \brief Virtual base class for managing individual levels.
* AmrLevel functions both as a container for state data on a level
* and also manages the advancement of data in time.
*/

class AmrLevel
{
    friend class Amr;
    friend class FillPatchIterator;
    friend class FillPatchIteratorHelper;

public:
    //! What time are we at?
    enum TimeLevel { AmrOldTime,
                     AmrHalfTime,
                     AmrNewTime,
                     Amr1QtrTime,
                     Amr3QtrTime,
                     AmrOtherTime };
    //! The destructor.
    virtual ~AmrLevel ();

    AmrLevel (const AmrLevel&) = delete;
    AmrLevel (AmrLevel&&) = delete;
    AmrLevel& operator= (const AmrLevel&) = delete;
    AmrLevel& operator= (AmrLevel&&) = delete;

    //! Get the level directory names
    void LevelDirectoryNames (const std::string &dir,
                              std::string &LevelDir,
                              std::string &FullPath) const;
    //! Create the Level_ directory for checkpoint and plot files
    virtual void CreateLevelDirectory (const std::string &dir);
    /**
    * \brief Set if the Level_ directory was created or to clear the value.
    * CreateLevelDirectory sets levelDirectoryCreated = true
    */
    void SetLevelDirectoryCreated(bool ldc) noexcept { levelDirectoryCreated = ldc; }
    /**
    * \brief A string written as the first item in writePlotFile() at
    * level zero.
    * It is so we can distinguish between different types of
    * plot files.
    * This default "HyperCLaw-V1.1" is for VisIt software and some of our
    * internal postprocessing routines
    */
    virtual std::string thePlotFileType () const
    {
        static const std::string the_plot_file_type("HyperCLaw-V1.1");
        return the_plot_file_type;
    }
    /**
    * \brief Write plot file stuff to specified directory.
    */
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                VisMF::How         how = VisMF::NFiles);

    //! Do pre-plotfile work to avoid synchronizations while writing the amr hierarchy
    virtual void writePlotFilePre (const std::string& dir,
                                   std::ostream&      os);

    //! Do post-plotfile work to avoid synchronizations while writing the amr hierarchy
    virtual void writePlotFilePost (const std::string& dir,
                                    std::ostream&      os);

    /**
    * \brief Write small plot file stuff to specified directory.
    */
    virtual void writeSmallPlotFile (const std::string& /*dir*/,
                                     std::ostream&      /*os*/,
                                     VisMF::How         /*how*/ = VisMF::NFiles) {}
    //! Write current state to checkpoint file.
    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os,
                             VisMF::How         how = VisMF::NFiles,
                             bool               dump_old = true);
    //! Do pre-checkpoint work to avoid synchronizations while writing the amr hierarchy
    virtual void checkPointPre (const std::string& dir,
                                std::ostream&      os);
    //! Do post-checkpoint work to avoid synchronizations while writing the amr hierarchy
    virtual void checkPointPost (const std::string& dir,
                                 std::ostream&      os);
    //! Restart from a checkpoint file.
    virtual void restart (Amr&          papa,
                          std::istream& is,
                          bool          bReadSpecial = false);

    //! Old checkpoint may have different number of states than the new source code.
    virtual void set_state_in_checkpoint (Vector<int>& state_in_checkpoint);

    //! Is name a state variable?
    static bool isStateVariable (const std::string& name,
                                int&               state_indx,
                                int&               ncomp);

    static void FlushFPICache ();
    /**
    * \brief Compute the initial time step.  This is a pure virtual function
    * and hence MUST be implemented by derived classes.
    */
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Vector<int>&           n_cycle,
                                   const Vector<IntVect>& ref_ratio,
                                   Vector<Real>&          dt_level,
                                   Real                  stop_time) = 0;
    /**
    * \brief Compute the next time step.  This is a pure virtual function
    * and hence MUST be implemented by derived classes.
    */
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_min,
                               Vector<Real>&          dt_level,
                               Real                  stop_time,
                               int                   post_regrid_flag) = 0;
    /**
    * \brief Do an integration step on this level.  Returns maximum safe
    * time step.  This is a pure virtual function and hence MUST
    * be implemented by derived classes.
    */
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) = 0;

    /**
    * \brief Contains operations to be done after a timestep.  If this
    * function is overridden, don't forget to reset FillPatcher.
    */
    virtual void post_timestep (int iteration);
    /**
    * \brief Contains operations to be done only after a full coarse
    * timestep.  The default implementation does nothing.
    */
    virtual void postCoarseTimeStep (Real time);
    /**
    * \brief Operations to be done after restart.
    */
    virtual void post_restart () {}
    /**
    * \brief Operations to be done after regridding
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual  void post_regrid (int lbase,
                               int new_finest) = 0;
    /**
    * \brief Operations to be done after initialization.
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual  void post_init (Real stop_time) = 0;
    /**
    * \brief Is it ok to continue the calculation?
    */
    virtual  int okToContinue () { return 1; }
    /**
    * \brief Should I regrid with this level as base level?
    * This test is only evaluated if regrid_int > 0 and
    * level_count >= regrid_int as well. Defaults to true.
    */
    virtual  int okToRegrid ();
    /**
    * \brief Init grid data at problem start-up.
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual void initData () = 0;
    //! Set the time levels of state data.
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new);
    //! Alloc space for old time data.
    virtual void allocOldData ();
    //! Delete old-time data.
    virtual void removeOldData ();
    /**
    * \brief Init data on this level from another AmrLevel (during regrid).
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual void init (AmrLevel &old) = 0;
    /**
    * Init data on this level after regridding if old AmrLevel
    * did not previously exist. This is a pure virtual function
    * and hence MUST be implemented by derived classes.
    */
    virtual void init () = 0;
    //! Reset data to initial time by swapping new and old time data.
    void reset ();
    //! Returns this AmrLevel.
    int Level () const noexcept { return level; }
    //! List of grids at this level.
    const BoxArray& boxArray () const noexcept { return grids; }
    const BoxArray& getEdgeBoxArray (int dir) const noexcept;
    const BoxArray& getNodalBoxArray () const noexcept;
    //
    const DistributionMapping& DistributionMap () const noexcept { return dmap; }
    //
    const FabFactory<FArrayBox>& Factory () const noexcept { return *m_factory; }
    //! Number of grids at this level.
    int numGrids () const noexcept { return static_cast<int>(grids.size()); }
    //! Number of states at this level.
    int numStates () const noexcept { return static_cast<int>(state.size()); }
    //! Returns the indices defining physical domain.
    const Box& Domain () const noexcept { return geom.Domain(); }
    //! Timestep n at this level.
    int nStep () const noexcept { return parent->levelSteps(level); }
    //! Returns the geometry object.
    const Geometry& Geom () const noexcept { return geom; }
    //
    const IntVect& fineRatio () const noexcept { return fine_ratio; }
    //! Returns number of cells on level.
    Long countCells () const noexcept;

    //! Get the area not to tag.
    const BoxArray& getAreaNotToTag () noexcept;
    const Box& getAreaToTag () noexcept;
    //! Construct the area not to tag.
    void constructAreaNotToTag ();
    //! Set the area not to tag.
    void setAreaNotToTag (BoxArray& ba) noexcept;

    void resetFillPatcher ();

    /**
    * \brief Error estimation for regridding. This is a pure virtual
    * function and hence MUST be implemented by derived classes.
    */
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
                           int          n_error_buf = 0,
                           int          ngrow = 0) = 0;
    //! Interpolate from coarse level to the valid area in mf.
    void FillCoarsePatch (MultiFab& mf,
                          int       dcomp,
                          Real      time,
                          int       state_idx,
                          int       scomp,
                          int       ncomp,
                          int       nghost = 0);
    //! Function to set physical boundary conditions.
    virtual void setPhysBoundaryValues (FArrayBox& dest,
                                        int        state_indx,
                                        Real       time,
                                        int        dest_comp,
                                        int        src_comp,
                                        int        num_comp);
    /**
    * \brief Returns a MultiFab containing the derived data for this level.
    * The user is responsible for deleting this pointer when done
    * with it.  If ngrow>0 the MultiFab is built on the appropriately
    * grown BoxArray.
    */
    virtual std::unique_ptr<MultiFab> derive (const std::string& name,
                                              Real               time,
                                              int                ngrow);
    /**
    * \brief This version of derive() fills the dcomp'th component of mf
    * with the derived quantity.
    */
    virtual void derive (const std::string& name,
                         Real               time,
                         MultiFab&          mf,
                         int                dcomp);
    //! State data object.
    StateData& get_state_data (int state_indx) noexcept { return state[state_indx]; }
    //! State data at old time.
    MultiFab& get_old_data (int state_indx) noexcept { return state[state_indx].oldData(); }
    //! State data at old time.
    const MultiFab& get_old_data (int state_indx) const noexcept { return state[state_indx].oldData(); }
    //! State data at new time.
    MultiFab& get_new_data (int state_indx) noexcept { return state[state_indx].newData(); }
    //! State data at new time.
    const MultiFab& get_new_data (int state_indx) const noexcept { return state[state_indx].newData(); }
    //! Returns list of Descriptors.
    static const DescriptorList& get_desc_lst () noexcept { return desc_lst; }
    //! Returns list of derived variables.
    static DeriveList& get_derive_lst () noexcept;
    //! Returns whether or not we want a post-timestep regrid.
    int postStepRegrid () const noexcept { return post_step_regrid; }
    //! Sets a new value for the post-timestep regrid trigger.
    void setPostStepRegrid (int new_val) noexcept { post_step_regrid = new_val; }

    //! Update the distribution maps in StateData based on the size of the map
    void UpdateDistributionMaps ( DistributionMapping& dmap );

    //! Boundary condition access function.
    Vector<int> getBCArray (int State_Type,
                            int gridno,
                            int strt_comp,
                            int ncomp);
    //! Get state data at specified index and time.
    MultiFab& get_data (int  state_indx, Real time) noexcept;
    //! Hack to allow override of (non-fine-fine) fillpatched boundary data
    virtual void set_preferred_boundary_values (MultiFab& S,
                                                int       state_index,
                                                int       scomp,
                                                int       dcomp,
                                                int       ncomp,
                                                Real      time) const;
    /**
    * \brief Called in grid_places after other tagging routines to modify
    * the list of tagged points.  Default implementation does nothing.
    */
    virtual void manual_tags_placement (TagBoxArray&    tags,
                                        const Vector<IntVect>& bf_lev);
    //! Modify list of variables to be plotted
    virtual void setPlotVariables ();
    //! Modify list of variables to be plotted
    virtual void setSmallPlotVariables ();
    /**
    * \brief Estimate the amount of work required to advance Just this level
    * based on the number of cells.
    * This estimate can be overwritten with different methods
    */
    virtual Real estimateWork();

    //! Which state data type is for work estimates? -1 means none
    virtual int WorkEstType () { return -1; }

    /**
    * \brief Returns one the TimeLevel enums.
    * Asserts that time is between AmrOldTime and AmrNewTime.
    */
    TimeLevel which_time (int  state_indx, Real time) const noexcept;

    //! Does the AmrLevel want Amr to write a plotfile now?
    virtual bool writePlotNow ();

    //! Does the AmrLevel want Amr to write a small plotfile now?
    virtual bool writeSmallPlotNow ();

#ifdef AMREX_PARTICLES
    //! This function can be called from the parent
    virtual void particle_redistribute (int /*lbase*/ = 0, bool /*a_init*/ = false) {;}
#endif

    /**
     * \brief Fill with FillPatcher on level > 0 and AmrLevel::FillPatch on level 0.
     *
     * \param mf     destination MultiFab
     * \param dcomp  starting component for the destination
     * \param ncomp  number of component to fill
     * \param nghost number of ghost cells to fill
     * \param time   time
     * \param state_index StateData index
     * \param scomp  starting component in the StateData
     */
    void FillPatcherFill (amrex::MultiFab& mf, int dcomp, int ncomp, int nghost,
                          amrex::Real time, int state_index, int scomp);

    static void FillPatch (AmrLevel& amrlevel,
                           MultiFab& leveldata,
                           int       boxGrow,
                           Real      time,
                           int       index,
                           int       scomp,
                           int       ncomp,
                           int       dcomp=0);

    static void FillPatchAdd (AmrLevel& amrlevel,
                              MultiFab& leveldata,
                              int       boxGrow,
                              Real      time,
                              int       index,
                              int       scomp,
                              int       ncomp,
                              int       dcomp=0);

    /**
     * \brief Evolve one step with Runge-Kutta (2, 3, or 4)
     *
     * To use RK, the StateData must have all the ghost cells needed.  See
     * namespace RungeKutta for expected function signatures of the callable
     * parameters.
     *
     * \param order      order of RK
     * \param state_type index of StateData
     * \param time       time at the beginning of the step.
     * \param dt         time step
     * \param iteration  iteration number on fine level during a coarse time
     *                   step.  For an AMR simulation with subcycling and a
     *                   refinement ratio of 2, the number is either 1 or 2,
     *                   denoting the first and second substep, respectively.
     * \param ncycle     number of subcyling steps.  It's usually 2 or 4.
     *                   Without subcycling, this will be 1.
     * \param f          computing right-hand side for evolving the StateData.
     *                   One can also register data for flux registers in this.
     * \param p          optionally post-processing RK stage results
     */
    template <typename F, typename P = RungeKutta::PostStageNoOp>
    void RK (int order, int state_type, Real time, Real dt, int iteration,
             int ncycle, F&& f, P&& p = RungeKutta::PostStageNoOp());

#ifdef AMREX_USE_EB
    static void SetEBMaxGrowCells (int nbasic, int nvolume, int nfull) noexcept {
        m_eb_basic_grow_cells = nbasic;
        m_eb_volume_grow_cells = nvolume;
        m_eb_full_grow_cells = nfull;
    }
    static int            m_eb_basic_grow_cells;
    static int            m_eb_volume_grow_cells;
    static int            m_eb_full_grow_cells;
    static void SetEBSupportLevel (EBSupport ebs) { m_eb_support_level = ebs; }
    static EBSupport      m_eb_support_level;
#endif

    //! Recommendation of a proper blocking factor
    static IntVect ProperBlockingFactor (AmrLevel const& amr_level, int boxGrow,
                                         IndexType const& boxType,
                                         StateDescriptor const& desc, int SComp);

protected:
    //! The constructors -- for derived classes.
    AmrLevel () noexcept {} // NOLINT

    AmrLevel (Amr&            papa,
              int             lev,
              const Geometry& level_geom,
              const BoxArray& ba,
              const DistributionMapping& dm,
              Real            time);

    //! Common code used by all constructors.
    void finishConstructor ();

    //
    // The Data.
    //
    int                   level{-1};    // AMR level (0 is coarsest).
    Geometry              geom;         // Geom at this level.
    BoxArray              grids;        // Cell-centered locations of grids.
    DistributionMapping   dmap;         // Distribution of grids among processes
    Amr*                  parent{nullptr};// Pointer to parent AMR structure.
    IntVect               crse_ratio;   // Refinement ratio to coarser level.
    IntVect               fine_ratio;   // Refinement ratio to finer level.
    static DeriveList     derive_lst;   // List of derived quantities.
    static DescriptorList desc_lst;     // List of state variables.
    Vector<StateData>     state;        // Array of state data.

    BoxArray m_AreaNotToTag; //Area which shouldn't be tagged on this level.
    Box      m_AreaToTag;    //Area which is allowed to be tagged on this level.

    int post_step_regrid{0}; // Whether or not to do a regrid after the timestep.

    bool levelDirectoryCreated{false}; // for checkpoints and plotfiles

    std::unique_ptr<FabFactory<FArrayBox> > m_factory;

    Vector<std::unique_ptr<FillPatcher<MultiFab>>> m_fillpatcher;

private:

    template <std::size_t order>
    void storeRKCoarseData (int state_type, Real time, Real dt,
                            MultiFab const& S_old,
                            Array<MultiFab,order> const& rkk);

    void FillRKPatch (int state_index, MultiFab& S, Real time,
                      int stage, int iteration, int ncycle);

    mutable BoxArray      edge_grids[AMREX_SPACEDIM];  // face-centered grids
    mutable BoxArray      nodal_grids;              // all nodal grids
};

//
// Forward declaration.
//
class FillPatchIteratorHelper;

class FillPatchIterator
    :
    public MFIter
{
  public:

    friend class AmrLevel;

    FillPatchIterator (AmrLevel& amrlevel,
                       MultiFab& leveldata);

    FillPatchIterator (AmrLevel& amrlevel,
                       MultiFab& leveldata,
                       int       boxGrow,
                       Real      time,
                       int       idx,
                       int       scomp,
                       int       ncomp);

    void Initialize (int  boxGrow,
                     Real time,
                     int  idx,
                     int  scomp,
                     int  ncomp);

    FArrayBox& operator() () noexcept { return m_fabs[MFIter::index()]; }

    Box UngrownBox () const noexcept { return MFIter::validbox(); }

    MultiFab& get_mf() noexcept { return m_fabs; }

private:

    void FillFromLevel0 (Real time, int idx, int scomp, int dcomp, int ncomp);
    void FillFromTwoLevels (Real time, int idx, int scomp, int dcomp, int ncomp);

    //
    // The data.
    //
    AmrLevel&                         m_amrlevel;
    MultiFab&                         m_leveldata;
    std::vector< std::pair<int,int> > m_range;
    MultiFab                          m_fabs;
    int                               m_ncomp;
};

class FillPatchIteratorHelper
{
public:

    friend class FillPatchIterator;

    FillPatchIteratorHelper (AmrLevel& amrlevel,
                             MultiFab& leveldata);

    FillPatchIteratorHelper (AmrLevel&     amrlevel,
                             MultiFab&     leveldata,
                             int           boxGrow,
                             Real          time,
                             int           index,
                             int           scomp,
                             int           ncomp,
                             InterpBase*   mapper);

    void Initialize (int           boxGrow,
                     Real          time,
                     int           idx,
                     int           scomp,
                     int           ncomp,
                     InterpBase*   mapper);

    void fill (FArrayBox& fab, int dcomp, int idx);

private:
    //
    // The data.
    //
    AmrLevel&                  m_amrlevel;
    MultiFab&                  m_leveldata;
    MultiFabCopyDescriptor     m_mfcd;
    Vector< Vector<MultiFabId> > m_mfid;     // [level][oldnew]
    Interpolater*              m_map = nullptr;
    std::map<int,Box>          m_ba;
    Real                       m_time = std::numeric_limits<Real>::lowest();
    int                        m_growsize = std::numeric_limits<int>::lowest();
    int                        m_index = -1;
    int                        m_scomp = -1;
    int                        m_ncomp = -1;
    bool                       m_FixUpCorners = false;

    std::map< int,Vector< Vector<Box> > >                m_fbox; // [grid][level][validregion]
    std::map< int,Vector< Vector<Box> > >                m_cbox; // [grid][level][fillablesubbox]
    std::map< int,Vector< Vector< Vector<FillBoxId> > > > m_fbid; // [grid][level][fillablesubbox][oldnew]
};

template <typename F, typename P>
void AmrLevel::RK (int order, int state_type, Real time, Real dt, int iteration,
                   int ncycle, F&& f, P&& p)
{
    BL_PROFILE("AmrLevel::RK()");

    AMREX_ASSERT(AmrLevel::desc_lst[state_type].nExtra() > 0); // Need ghost cells in StateData

    MultiFab& S_old = get_old_data(state_type);
    MultiFab& S_new = get_new_data(state_type);
    const Real t_old = state[state_type].prevTime();
    const Real t_new = state[state_type].curTime();
    AMREX_ALWAYS_ASSERT(amrex::almostEqual(time,t_old,10) && amrex::almostEqual(time+dt,t_new,10));

    if (order == 2) {
        RungeKutta::RK2(S_old, S_new, time, dt, std::forward<F>(f),
                        [&] (int /*stage*/, MultiFab& mf, Real t) {
                            FillPatcherFill(mf, 0, mf.nComp(), mf.nGrow(), t,
                                            state_type, 0); },
                        std::forward<P>(p));
    } else if (order == 3) {
        RungeKutta::RK3(S_old, S_new, time, dt, std::forward<F>(f),
                        [&] (int stage, MultiFab& mf, Real t) {
                            FillRKPatch(state_type, mf, t, stage, iteration, ncycle);
                        },
                        [&] (Array<MultiFab,3> const& rkk) {
                            if (level < parent->finestLevel()) {
                                storeRKCoarseData(state_type, time, dt, S_old, rkk);
                            }
                        },
                        std::forward<P>(p));
    } else if (order == 4) {
        RungeKutta::RK4(S_old, S_new, time, dt, std::forward<F>(f),
                        [&] (int stage, MultiFab& mf, Real t) {
                            FillRKPatch(state_type, mf, t, stage, iteration, ncycle);
                        },
                        [&] (Array<MultiFab,4> const& rkk) {
                            if (level < parent->finestLevel()) {
                                storeRKCoarseData(state_type, time, dt, S_old, rkk);
                            }
                        },
                        std::forward<P>(p));
    } else {
        amrex::Abort("AmrLevel::RK: order = "+std::to_string(order)+" is not supported");
    }
}

template <std::size_t order>
void AmrLevel::storeRKCoarseData (int state_type, Real time, Real dt,
                                  MultiFab const& S_old,
                                  Array<MultiFab,order> const& rkk)
{
    if (level == parent->finestLevel()) { return; }

    const StateDescriptor& desc = AmrLevel::desc_lst[state_type];

    auto& fillpatcher = parent->getLevel(level+1).m_fillpatcher[state_type];
    fillpatcher = std::make_unique<FillPatcher<MultiFab>>
        (parent->boxArray(level+1), parent->DistributionMap(level+1),
         parent->Geom(level+1),
         parent->boxArray(level), parent->DistributionMap(level),
         parent->Geom(level),
         IntVect(desc.nExtra()), desc.nComp(), desc.interp(0));

    fillpatcher->storeRKCoarseData(time, dt, S_old, rkk);
}


}

#endif /*_AmrLevel_H_*/
