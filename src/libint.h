#include <stdbool.h>

typedef unsigned long ulong;

typedef struct SystemData SystemData;

SystemData* init(const char *system, const char *basis); void deinit(SystemData *sys);

ulong nbf(SystemData *sys);

void coulomb(double *I, SystemData *sys);
void kinetic(double *I, SystemData *sys);
void nuclear(double *I, SystemData *sys);
void overlap(double *I, SystemData *sys);
