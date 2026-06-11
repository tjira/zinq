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

        std::vector<libint2::Atom> atoms = libint2::read_dotxyz(file);

        return new SystemData{atoms, libint2::BasisSet(basis, atoms)};
    }

    void deinit(SystemData *sys) {
        delete sys;
    }

    ulong nbf(SystemData *sys) {
        return sys->obs.nbf();
    }

    void oneelec(double *I, libint2::Engine &engine, const BasisSet &obs) {
        std::vector<Engine> engines(1, engine); size_t nbf = obs.nbf(); auto sh2bf = obs.shell2bf();

        for (size_t i = 0; i < obs.size(); i++) {
            for (size_t j = i; j < obs.size(); j++) {

                int id = 0; int idx = 0;

                engines.at(id).compute(obs.at(i), obs.at(j));

                if (engines.at(id).results().at(0) == nullptr) continue;

                for (size_t k = 0; k < obs.at(i).size(); k++) {
                    for (size_t l = 0; l < obs.at(j).size(); l++) {
                        I[(k + sh2bf.at(i)) * nbf + (l + sh2bf.at(j))] = engines.at(id).results().at(0)[idx  ];
                        I[(l + sh2bf.at(j)) * nbf + (k + sh2bf.at(i))] = engines.at(id).results().at(0)[idx++];
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
                            for (size_t n = 0; n < obs.at(j).size(); n++) {
                                for (size_t o = 0; o < obs.at(k).size(); o++) {
                                    for (size_t p = 0; p < obs.at(l).size(); p++) {

                                        size_t bf1 = m + sh2bf.at(i);
                                        size_t bf2 = n + sh2bf.at(j);
                                        size_t bf3 = o + sh2bf.at(k);
                                        size_t bf4 = p + sh2bf.at(l);

                                        size_t idx1 = bf1 * nbf * nbf * nbf + bf2 * nbf * nbf + bf3 * nbf + bf4;
                                        size_t idx2 = bf1 * nbf * nbf * nbf + bf2 * nbf * nbf + bf4 * nbf + bf3;
                                        size_t idx3 = bf2 * nbf * nbf * nbf + bf1 * nbf * nbf + bf3 * nbf + bf4;
                                        size_t idx4 = bf2 * nbf * nbf * nbf + bf1 * nbf * nbf + bf4 * nbf + bf3;
                                        size_t idx5 = bf3 * nbf * nbf * nbf + bf4 * nbf * nbf + bf1 * nbf + bf2;
                                        size_t idx6 = bf3 * nbf * nbf * nbf + bf4 * nbf * nbf + bf2 * nbf + bf1;
                                        size_t idx7 = bf4 * nbf * nbf * nbf + bf3 * nbf * nbf + bf1 * nbf + bf2;
                                        size_t idx8 = bf4 * nbf * nbf * nbf + bf3 * nbf * nbf + bf2 * nbf + bf1;

                                        I[idx1] = engines.at(id).results().at(0)[idx  ];
                                        I[idx2] = engines.at(id).results().at(0)[idx  ];
                                        I[idx3] = engines.at(id).results().at(0)[idx  ];
                                        I[idx4] = engines.at(id).results().at(0)[idx  ];
                                        I[idx5] = engines.at(id).results().at(0)[idx  ];
                                        I[idx6] = engines.at(id).results().at(0)[idx  ];
                                        I[idx7] = engines.at(id).results().at(0)[idx  ];
                                        I[idx8] = engines.at(id).results().at(0)[idx++];
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
