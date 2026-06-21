typedef unsigned long ulong;

typedef struct SystemData SystemData;

SystemData* init(const char *system, const char *basis); void deinit(SystemData *sys);

ulong nbf(SystemData *sys);
ulong nat(SystemData *sys);

void coulomb(double *I, SystemData *sys);
void kinetic(double *I, SystemData *sys);
void nuclear(double *I, SystemData *sys);
void overlap(double *I, SystemData *sys);

void coulomb_deriv(double *I, SystemData *sys);
void kinetic_deriv(double *I, SystemData *sys);
void nuclear_deriv(double *I, SystemData *sys);
void overlap_deriv(double *I, SystemData *sys);

void atoms(int    *atoms, SystemData *sys);
void coors(double *coors, SystemData *sys);

void evaluate_basis           (double *phi,                                                    double x, double y, double z, SystemData *sys);
void evaluate_basis_derivative(double *phi, double *dphi_dx, double *dphi_dy, double *dphi_dz, double x, double y, double z, SystemData *sys);

