#include "libint2/config.h"

#include <stddef.h>

typedef struct SystemData SystemData;

/// Initializes the Libint molecular integral environment, atomic coordinates, and basis set.
SystemData* libint_init(const char *system, const char *basis, double bohr_to_angstrom);

/// Releases Libint molecular system resources and frees allocated memory.
void libint_deinit(SystemData *sys);

/// Creates a deep copy of the molecular system and orbital basis structures.
SystemData* libint_clone(SystemData *sys);

/// Returns the total number of contracted basis functions in the primary basis.
size_t libint_nbf(SystemData *sys);

/// Returns the number of atoms in the molecular system.
size_t libint_nat(SystemData *sys);

/// Returns the number of contracted basis functions in the auxiliary basis set.
size_t libint_aux_nbf(SystemData *sys, const char *aux_basis);

/// Computes the generalized Hartree-Fock two-electron Fock matrix from the density matrix.
void libint_fock_ghf(double *F, const double *P, double exch_factor, SystemData *sys, size_t nthreads);

/// Computes the restricted Hartree-Fock two-electron Fock matrix from the density matrix.
void libint_fock_rhf(double *F, const double *P, double exch_factor, SystemData *sys, size_t nthreads);

/// Computes four-center two-electron electron repulsion integrals over contracted basis functions.
void libint_coulomb(double *I, SystemData *sys, size_t nthreads);

/// Computes one-electron kinetic energy integrals over contracted basis functions.
void libint_kinetic(double *I, SystemData *sys, size_t nthreads);

/// Computes one-electron nuclear attraction integrals over contracted basis functions.
void libint_nuclear(double *I, SystemData *sys, size_t nthreads);

/// Computes one-electron overlap integrals between contracted basis functions.
void libint_overlap(double *I, SystemData *sys, size_t nthreads);

/// Computes two-center two-electron Coulomb metric integrals over auxiliary basis functions.
void libint_coulomb_2c(double *J, SystemData *sys, const char *aux_basis, size_t nthreads);

/// Computes three-center two-electron Coulomb integrals over orbital and auxiliary bases.
void libint_coulomb_3c(double *I, SystemData *sys, const char *aux_basis, size_t nthreads);

/// Computes the overlap integral matrix between basis sets of two distinct molecular systems.
void libint_overlap_2sys(double *I, SystemData *sys1, SystemData *sys2, size_t nthreads);

/// Computes the generalized Hartree-Fock two-electron energy gradient with respect to nuclear coordinates.
void libint_coulomb_gradient_ghf(double *G, const double *P, double exch_factor, SystemData *sys, size_t nthreads);

/// Computes the restricted Hartree-Fock two-electron energy gradient with respect to nuclear coordinates.
void libint_coulomb_gradient_rhf(double *G, const double *P, double exch_factor, SystemData *sys, size_t nthreads);

/// Computes first geometric derivatives of four-center two-electron Coulomb repulsion integrals.
void libint_coulomb_deriv(double *I, SystemData *sys, size_t nthreads);

/// Computes first geometric derivatives of one-electron kinetic energy integrals.
void libint_kinetic_deriv(double *I, SystemData *sys, size_t nthreads);

/// Computes first geometric derivatives of one-electron nuclear attraction integrals.
void libint_nuclear_deriv(double *I, SystemData *sys, size_t nthreads);

/// Computes first geometric derivatives of one-electron overlap integrals.
void libint_overlap_deriv(double *I, SystemData *sys, size_t nthreads);

/// Computes first geometric derivatives of two-center two-electron Coulomb metric integrals.
void libint_coulomb_2c_deriv(double *I, SystemData *sys, const char *aux_basis, size_t nthreads);

/// Computes first geometric derivatives of three-center two-electron Coulomb integrals.
void libint_coulomb_3c_deriv(double *I, SystemData *sys, const char *aux_basis, size_t nthreads);

/// Retrieves atomic numbers for all atoms in the molecular system.
void libint_atoms(int *atoms, SystemData *sys);

/// Retrieves Cartesian coordinates in Bohr for all atoms in the molecular system.
void libint_coors(double *coors, SystemData *sys);

/// Maps each basis function index to its parent atom index in the molecular system.
void libint_bf2at(int *bf2at, SystemData *sys);

/// Updates nuclear Cartesian coordinates and recalculates basis shell center origins.
void libint_update_coords(SystemData *sys, const double *coords);

/// Evaluates molecular orbital basis function amplitudes at a Cartesian spatial coordinate.
void libint_evaluate_basis_v(double *phi, double x, double y, double z, SystemData *sys);

/// Evaluates basis function amplitudes and spatial gradient components at a Cartesian coordinate.
void libint_evaluate_basis_vg(double *phi, double *dphi_dx, double *dphi_dy, double *dphi_dz, double x, double y, double z, SystemData *sys);

/// Evaluates basis function amplitudes, spatial gradients, and Laplacian values at a point.
void libint_evaluate_basis_vgl(double *phi, double *dphi_dx, double *dphi_dy, double *dphi_dz, double *lapl_phi, double x, double y, double z, SystemData *sys);

/// Computes the number of atomic orbitals for an element in the STO-3G minimal basis.
size_t libint_num_ao(int Z);

/// Constructs the initial Superposition of Atomic Densities (SAD) guess density matrix.
void libint_soad_density(double *D, SystemData *sys);
