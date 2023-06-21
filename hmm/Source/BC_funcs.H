#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>


// Functions that implement the boundary conditions. 
// These functions are used to construct a "BndryFunc" object.
// The "BndryFunc" object is then supplied to the setComponent function when defining a state variable. 


// boundary function.
// - First fills all ghost cells, ignoring external BC 
// - then calls custom function to overwrite external BC ghost cells 
void hmm_statefill(amrex::Box const& bx, amrex::FArrayBox& data,
                  const int dcomp, const int numcomp,
                  amrex::Geometry const& geom, const amrex::Real time,
                  const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                  const int scomp);

// boundary function v2
// - fills external BC ghost cells in one sweep 
void hmm_statefill_v2(amrex::Box const& bx, amrex::FArrayBox& data,
                  const int dcomp, const int numcomp,
                  amrex::Geometry const& geom, const amrex::Real time,
                  const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                  const int scomp);

