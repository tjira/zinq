# Population Analysis

Population analysis is a mathematical procedure used to partition the total electronic density of a molecule among its constituent atoms or basis functions. This enables the calculation of effective atomic charges, bond orders, and orbital populations, providing physical insights into chemical bonding and polarization.

---

## I. Mulliken Population Analysis

### 1. Basis Function and Atomic Populations

Mulliken population analysis partitions the total electronic charge by utilizing the non-orthogonal atomic orbital basis set. The net population $N_\mu$ of basis function $\chi_\mu$ is computed from the density matrix $\mathbf{P}$ and overlap matrix $\mathbf{S}$ as

$$
N_\mu=\sum_\nu P_{\mu\nu}S_{\mu\nu}
$$

which corresponds to the diagonal element of the matrix product $\mathbf{P}\mathbf{S}$. To avoid the expensive computational cost of performing a full matrix multiplication, the codebase implements a trick where it directly evaluates the dot product of the corresponding rows of the density and overlap matrices, extracting the diagonal elements in linear time. The total electron population $N_A$ associated with atom $A$ is then obtained by summing the contributions of all basis functions belonging to that atom as

$$
N_A=\sum_{\mu\in A}N_\mu
$$

where the shared overlap populations are divided equally between the participating atoms. The codebase maps basis functions to their corresponding atomic centers using a lookup array `sys.bf2at` representing the basis function to atom association.

### 2. Net Atomic Charges

The net atomic charge $q_A$ of atom $A$ is calculated by subtracting the total electron population $N_A$ from the nuclear charge $Z_A$ as

$$
q_A=Z_A-N_A
$$

which provides a quantitative measure of charge transfer and molecular dipole moments.
