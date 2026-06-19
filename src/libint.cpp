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

                const auto& res = engines.at(id).results();

                if (res.at(0) == nullptr) continue;

                for (size_t k = 0; k < obs.at(i).size(); k++) {
                    for (size_t l = 0; l < obs.at(j).size(); l++) {
                        double val = res.at(0)[idx++];

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

                        const auto& res = engines.at(id).results();

                        if (res.at(0) == nullptr) continue;

                        for (size_t m = 0; m < obs.at(i).size(); m++) {
                            size_t bf1 = m + sh2bf.at(i);

                            for (size_t n = 0; n < obs.at(j).size(); n++) {
                                size_t bf2 = n + sh2bf.at(j);

                                for (size_t o = 0; o < obs.at(k).size(); o++) {
                                    size_t bf3 = o + sh2bf.at(k);

                                    for (size_t p = 0; p < obs.at(l).size(); p++) {
                                        size_t bf4 = p + sh2bf.at(l);

                                        double val = res.at(0)[idx++];

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

extern "C" {
    using namespace libint2; typedef unsigned long ulong;

    void oneelec_deriv(double *I, libint2::Engine &engine, const BasisSet &obs, const std::vector<Atom> &atoms) {
        std::vector<Engine> engines(1, engine); bool is_nuclear = engine.oper() == Operator::nuclear;

        size_t nbf = obs   .nbf(); auto sh2bf = obs  .shell2bf(     );
        size_t nat = atoms.size(); auto sh2at = obs.shell2atom(atoms);

        for (size_t i = 0; i < obs.size(); i++) {
            size_t atom_i = sh2at[i];

            for (size_t j = i; j < obs.size(); j++) {
                size_t atom_j = sh2at[j];

                int id = 0; int idx = 0;

                engines.at(id).compute(obs.at(i), obs.at(j));

                const auto& res = engines.at(id).results();

                if (res.at(0) == nullptr) continue;

                for (size_t k = 0; k < obs.at(i).size(); k++) {
                    size_t bf1 = k + sh2bf.at(i);

                    for (size_t l = 0; l < obs.at(j).size(); l++) {
                        size_t bf2 = l + sh2bf.at(j);

                        for (size_t c = 0; c < nat; c++) for (size_t d = 0; d < 3; d++) {
                            double val = 0;

                            if (is_nuclear) {
                                val += res.at(6 + 3 * c + d)[idx];
                            }

                            if (atom_i == c) val += res.at(d + 0)[idx];
                            if (atom_j == c) val += res.at(d + 3)[idx];

                            I[(3 * c + d) * nbf * nbf + bf1 * nbf + bf2] = val;

                            if (bf1 != bf2) {
                                I[(3 * c + d) * nbf * nbf + bf2 * nbf + bf1] = val;
                            }
                        }

                        idx++;
                    }
                }
            }
        }
    }

    void twoelec_deriv(double *I, libint2::Engine &engine, const BasisSet &obs, const std::vector<Atom> &atoms) {
        std::vector<libint2::Engine> engines(1, engine);

        size_t nbf = obs   .nbf(); auto sh2bf = obs  .shell2bf(     );
        size_t nat = atoms.size(); auto sh2at = obs.shell2atom(atoms);

        for (size_t i = 0; i < obs.size(); i++) {
            size_t atom_i = sh2at[i];

            for (size_t j = i; j < obs.size(); j++) {
                size_t atom_j = sh2at[j];

                for (size_t k = i; k < obs.size(); k++) {
                    size_t atom_k = sh2at[k];

                    for (size_t l = (i == k ? j : k); l < obs.size(); l++) {
                        size_t atom_l = sh2at[l];

                        int id = 0;  int idx = 0;

                        engines.at(id).compute(obs.at(i), obs.at(j), obs.at(k), obs.at(l));

                        const auto& res = engines.at(id).results();

                        if (res.at(0) == nullptr) continue;

                        auto apply = [&](size_t p, size_t q, size_t r, size_t s) {
                            size_t base_idx = p * nbf * nbf * nbf + q * nbf * nbf + r * nbf + s;

                            for (size_t c = 0; c < nat; c++) for (size_t d = 0; d < 3; d++) {
                                double val = 0;

                                if (atom_i == c) val += res.at(d + 0)[idx];
                                if (atom_j == c) val += res.at(d + 3)[idx];
                                if (atom_k == c) val += res.at(d + 6)[idx];
                                if (atom_l == c) val += res.at(d + 9)[idx];

                                I[(3 * c + d) * nbf * nbf * nbf * nbf + base_idx] = val;
                            }
                        };

                        for (size_t m = 0; m < obs.at(i).size(); m++) {
                            size_t bf1 = m + sh2bf.at(i);

                            for (size_t n = 0; n < obs.at(j).size(); n++) {
                                size_t bf2 = n + sh2bf.at(j);

                                for (size_t o = 0; o < obs.at(k).size(); o++) {
                                    size_t bf3 = o + sh2bf.at(k);

                                    for (size_t p_val = 0; p_val < obs.at(l).size(); p_val++) {
                                        size_t bf4 = p_val + sh2bf.at(l);

                                        apply(bf1, bf3, bf2, bf4); 
                                        
                                        if (bf3 != bf4) {
                                            apply(bf1, bf4, bf2, bf3); 
                                        }
                                        
                                        if (bf1 != bf2) {
                                            apply(bf2, bf3, bf1, bf4); 
                                        }
                                        
                                        if (bf1 != bf2 && bf3 != bf4) {
                                            apply(bf2, bf4, bf1, bf3); 
                                        }

                                        if (!(bf1 == bf3 && bf2 == bf4) && !(bf1 == bf4 && bf2 == bf3)) {
                                            
                                            apply(bf3, bf1, bf4, bf2); 
                                            
                                            if (bf3 != bf4) {
                                                apply(bf4, bf1, bf3, bf2); 
                                            }
                                            
                                            if (bf1 != bf2) {
                                                apply(bf3, bf2, bf4, bf1); 
                                            }
                                            
                                            if (bf1 != bf2 && bf3 != bf4) {
                                                apply(bf4, bf2, bf3, bf1); 
                                            }
                                        }

                                        idx++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void coulomb_deriv(double *I, SystemData *sys) {
        if (!sys) return; libint2::initialize();

        Engine engine(Operator::coulomb, sys->obs.max_nprim(), sys->obs.max_l(), 1, 1e-14);

        twoelec_deriv(I, engine, sys->obs, sys->atoms); libint2::finalize();
    }

    void kinetic_deriv(double *I, SystemData *sys) {
        if (!sys) return; libint2::initialize();

        Engine engine(Operator::kinetic, sys->obs.max_nprim(), sys->obs.max_l(), 1, 1e-14);

        oneelec_deriv(I, engine, sys->obs, sys->atoms); libint2::finalize();
    }

    void overlap_deriv(double *I, SystemData *sys) {
        if (!sys) return; libint2::initialize();

        Engine engine(Operator::overlap, sys->obs.max_nprim(), sys->obs.max_l(), 1, 1e-14);

        oneelec_deriv(I, engine, sys->obs, sys->atoms); libint2::finalize();
    }

    void nuclear_deriv(double *I, SystemData *sys) {
        if (!sys) return; libint2::initialize();

        Engine engine(Operator::nuclear, sys->obs.max_nprim(), sys->obs.max_l(), 1, 1e-14);

        engine.set_params(make_point_charges(sys->atoms));

        oneelec_deriv(I, engine, sys->obs, sys->atoms); libint2::finalize();
    }
}

extern "C" {
    using namespace libint2; using namespace libint2::solidharmonics;

    void evaluate_basis(double *phi, double x, double y, double z, SystemData *sys) {
        if (!sys) return;

        size_t bf_idx = 0; int max_nprim = sys->obs.max_nprim();

        std::vector<double> px(sys->obs.max_l() + 1);
        std::vector<double> py(sys->obs.max_l() + 1);
        std::vector<double> pz(sys->obs.max_l() + 1);

        int max_ncart = (sys->obs.max_l() + 1) * (sys->obs.max_l() + 2) / 2;

        std::vector<double> exp_vals(max_nprim);
        std::vector<double> crt_vals(max_ncart);

        for (const auto& shell : sys->obs) {
            double dx = x - shell.O[0];
            double dy = y - shell.O[1];
            double dz = z - shell.O[2];

            double r2 = dx * dx + dy * dy + dz * dz;

            for (size_t i = 0; i < shell.alpha.size(); i++) {
                exp_vals[i] = std::exp(-shell.alpha[i] * r2);
            }

            for (size_t i = 0; i < shell.contr.size(); i++) {
                int L = shell.contr[i].l; double sum = 0;

                for (size_t j = 0; j < shell.alpha.size(); j++) {
                    sum += shell.contr[i].coeff[j] * exp_vals[j];
                }

                px[0] = 1; py[0] = 1; pz[0] = 1;

                for (int j = 1; j <= L; j++) {
                    px[j] = px[j - 1] * dx;
                    py[j] = py[j - 1] * dy;
                    pz[j] = pz[j - 1] * dz;
                }

                for (int j = L, c = 0; j >= 0; j--) for (int k = L - j; k >= 0; k--) {
                    crt_vals[c++] = sum * px[j] * py[k] * pz[L - j - k];
                }

                if (shell.contr[i].pure) {
                    const auto& coefs = SolidHarmonicsCoefficients<double>::instance(L);

                    for (int j = 0; j < 2 * L + 1; j++) {
                        double val = 0.0;

                        for (unsigned int k = 0; k < coefs.nnz(j); k++) {
                            val += coefs.row_values(j)[k] * crt_vals[coefs.row_idx(j)[k]];
                        }

                        phi[bf_idx++] = val;
                    }
                }

                else for (int j = 0; j < (L + 1) * (L + 2) / 2; j++) {
                    phi[bf_idx++] = crt_vals[j];
                }
            }
        }
    }
}
