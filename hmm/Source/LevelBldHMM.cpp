
#include <AMReX_LevelBld.H>
#include <HMM.H>

using namespace amrex;

class LevelBldAdv
    :
    public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
				  const DistributionMapping& dm,
                                  Real            time) override;
};

LevelBldAdv Adv_bld;

LevelBld*
getLevelBld ()
{
    return &Adv_bld;
}

void
LevelBldAdv::variableSetUp ()
{
    HMM::variableSetUp();
}

void
LevelBldAdv::variableCleanUp ()
{
    HMM::variableCleanUp();
}

AmrLevel*
LevelBldAdv::operator() ()
{
    return new HMM;
}

AmrLevel*
LevelBldAdv::operator() (Amr&            papa,
        	   	         int             lev,
                         const Geometry& level_geom,
                         const BoxArray& ba,
                         const DistributionMapping& dm,
                         Real            time)
{
    return new HMM(papa, lev, level_geom, ba, dm, time);
}
