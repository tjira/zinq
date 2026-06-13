#include <libint2.hpp>

#include <fstream>

struct SystemData {
    std::vector<libint2::Atom> atoms; libint2::BasisSet obs;
};

extern "C" {
    using namespace libint2; typedef unsigned long ulong;

    SystemData* init(const char *system, const char *basis) {
        std::ifstream file(system);

        if (!file.is_open()) return nullptr;

        std::vector<Atom> atoms = read_dotxyz(file);

        return new SystemData{atoms, BasisSet(atoms, BasisSet::read_g94_basis_library(basis))};
    }

    void deinit(SystemData *sys) {
        if (!sys) return; delete sys;
    }

    ulong nat(SystemData *sys) {
        if (!sys) return 0; return sys->atoms.size();
    }

    ulong nbf(SystemData *sys) {
        if (!sys) return 0; return sys->obs.nbf();
    }

    void atoms(int *atoms, SystemData *sys) {
        if (!sys) return;

        for (size_t i = 0; i < sys->atoms.size(); ++i) {
            atoms[i] = sys->atoms[i].atomic_number;
        }
    }

    void coors(double *coords, SystemData *sys) {
        if (!sys) return;

        for (size_t i = 0; i < sys->atoms.size(); i++) {
            coords[3 * i + 0] = sys->atoms[i].x;
            coords[3 * i + 1] = sys->atoms[i].y;
            coords[3 * i + 2] = sys->atoms[i].z;
        }
    }
}

extern "C" {
    using namespace libint2; typedef unsigned long ulong;

    void oneelec(double *I, libint2::Engine &engine, const BasisSet &obs) {
        std::vector<Engine> engines(1, engine); size_t nbf = obs.nbf(); auto sh2bf = obs.shell2bf();

        for (size_t i = 0; i < obs.size(); i++) {
            for (size_t j = i; j < obs.size(); j++) {

                int id = 0; int idx = 0;

                engines.at(id).compute(obs.at(i), obs.at(j));

                if (engines.at(id).results().at(0) == nullptr) continue;

                for (size_t k = 0; k < obs.at(i).size(); k++) {
                    for (size_t l = 0; l < obs.at(j).size(); l++) {
                        double val = engines.at(id).results().at(0)[idx++];

                        I[(k + sh2bf.at(i)) * nbf + (l + sh2bf.at(j))] = val;
                        I[(l + sh2bf.at(j)) * nbf + (k + sh2bf.at(i))] = val;
                    }
                }
            }
        }
    }

    void twoelec(double *I, libint2::Engine &engine, const BasisSet &obs) {
        std::vector<Engine> engines(1, engine); size_t nbf = obs.nbf(); auto sh2bf = obs.shell2bf();

        for (size_t i = 0; i < obs.size(); i++) {
            for (size_t j = i; j < obs.size(); j++) {
                for (size_t k = i; k < obs.size(); k++) {

                    for (size_t l = (i == k ? j : k); l < obs.size(); l++) {

                        int id = 0; int idx = 0;

                        engines.at(id).compute(obs.at(i), obs.at(j), obs.at(k), obs.at(l));

                        if (engines.at(id).results().at(0) == nullptr) continue;

                        for (size_t m = 0; m < obs.at(i).size(); m++) {
                            size_t bf1 = m + sh2bf.at(i);

                            for (size_t n = 0; n < obs.at(j).size(); n++) {
                                size_t bf2 = n + sh2bf.at(j);

                                for (size_t o = 0; o < obs.at(k).size(); o++) {
                                    size_t bf3 = o + sh2bf.at(k);

                                    for (size_t p = 0; p < obs.at(l).size(); p++) {
                                        size_t bf4 = p + sh2bf.at(l);

                                        double val = engines.at(id).results().at(0)[idx++];

                                        I[bf1 * nbf * nbf * nbf + bf3 * nbf * nbf + bf2 * nbf + bf4] = val;
                                        I[bf1 * nbf * nbf * nbf + bf4 * nbf * nbf + bf2 * nbf + bf3] = val;
                                        I[bf2 * nbf * nbf * nbf + bf3 * nbf * nbf + bf1 * nbf + bf4] = val;
                                        I[bf2 * nbf * nbf * nbf + bf4 * nbf * nbf + bf1 * nbf + bf3] = val;
                                        I[bf3 * nbf * nbf * nbf + bf1 * nbf * nbf + bf4 * nbf + bf2] = val;
                                        I[bf3 * nbf * nbf * nbf + bf2 * nbf * nbf + bf4 * nbf + bf1] = val;
                                        I[bf4 * nbf * nbf * nbf + bf1 * nbf * nbf + bf3 * nbf + bf2] = val;
                                        I[bf4 * nbf * nbf * nbf + bf2 * nbf * nbf + bf3 * nbf + bf1] = val;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void coulomb(double *I, SystemData *sys) {
        if (!sys) return; libint2::initialize();

        Engine engine(Operator::coulomb, sys->obs.max_nprim(), sys->obs.max_l(), 0, 1e-14);

        twoelec(I, engine, sys->obs); libint2::finalize();
    }

    void kinetic(double *I, SystemData *sys) {
        if (!sys) return; libint2::initialize();

        Engine engine(Operator::kinetic, sys->obs.max_nprim(), sys->obs.max_l(), 0, 1e-14);

        oneelec(I, engine, sys->obs); libint2::finalize();
    }

    void nuclear(double *I, SystemData *sys) {
        if (!sys) return; libint2::initialize();

        Engine engine(Operator::nuclear, sys->obs.max_nprim(), sys->obs.max_l(), 0, 1e-14);

        engine.set_params(make_point_charges(sys->atoms));

        oneelec(I, engine, sys->obs); libint2::finalize();
    }

    void overlap(double *I, SystemData *sys) {
        if (!sys) return; libint2::initialize();

        Engine engine(Operator::overlap, sys->obs.max_nprim(), sys->obs.max_l(), 0, 1e-14);

        oneelec(I, engine, sys->obs); libint2::finalize();
    }
}
