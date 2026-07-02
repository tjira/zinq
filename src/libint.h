#include <stddef.h>

typedef struct SystemData SystemData;

SystemData* libint_init(const char *system, const char *basis); void libint_deinit(SystemData *sys);

size_t libint_nbf(SystemData *sys);
size_t libint_nat(SystemData *sys);

void libint_coulomb(double *I, SystemData *sys);
void libint_kinetic(double *I, SystemData *sys);
void libint_nuclear(double *I, SystemData *sys);
void libint_overlap(double *I, SystemData *sys);

void libint_coulomb_deriv(double *I, SystemData *sys);
void libint_kinetic_deriv(double *I, SystemData *sys);
void libint_nuclear_deriv(double *I, SystemData *sys);
void libint_overlap_deriv(double *I, SystemData *sys);

void libint_atoms(int    *atoms, SystemData *sys);
void libint_coors(double *coors, SystemData *sys);
void libint_bf2at(int    *bf2at, SystemData *sys);

void libint_evaluate_basis_v  (double *phi,                                                                      double x, double y, double z, SystemData *sys);
void libint_evaluate_basis_vg (double *phi, double *dphi_dx, double *dphi_dy, double *dphi_dz,                   double x, double y, double z, SystemData *sys);
void libint_evaluate_basis_vgl(double *phi, double *dphi_dx, double *dphi_dy, double *dphi_dz, double *lapl_phi, double x, double y, double z, SystemData *sys);
