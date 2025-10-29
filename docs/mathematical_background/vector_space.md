---
title: Vector Spaces and Linear Operators
parent: Mathematical Background
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Vector Spaces and Linear Operators<!--\label{sec:vector_spaces}-->
Many concepts in quantum mechanics rest on the algebra of vector spaces and the action of linear operators. We review what is needed here. Fuller accounts can be found in standard linear algebra and functional analysis texts. Throughout, the scalar field is a general field $$F$$, but in quantum mechanics we take $$F=\mathbb{C}$$ unless stated otherwise.

## Vector Spaces and Basic Properties

The object of interest is a vector space: a set with addition and scalar multiplication satisfying compatibility rules. The following definition states these axioms explicitly.

{:.definition}
> A vector space $$V$$ over a field $$F$$ is a non-empty set equipped with vector addition and scalar multiplication, which must satisfy the following axioms for all $$\mathbf{u},\mathbf{v},\mathbf{w}\in V$$ and all scalars $$a,b\in F$$:
>
> 1. Associativity of vector addition: $$(\mathbf{u}+\mathbf{v})+\mathbf{w}=\mathbf{u}+(\mathbf{v}+\mathbf{w})$$
> 2. Commutativity of vector addition: $$\mathbf{u}+\mathbf{v}=\mathbf{v}+\mathbf{u}$$
> 3. Identity element of vector addition: $$\exists\mathbf{0}\in V$$ such that $$\mathbf{u}+\mathbf{0}=\mathbf{u}$$
> 4. Inverse element of vector addition: $$\exists(-\mathbf{u})\in V$$ such that $$\mathbf{u}+(-\mathbf{u})=\mathbf{0}$$
> 5. Compatibility of scalar multiplication with field multiplication: $$a(b\mathbf{u})=(ab)\mathbf{u}$$
> 6. Identity element of scalar multiplication: $$1\mathbf{u}=\mathbf{u}$$, where $$1$$ is the multiplicative identity in $$F$$
> 7. Distributivity of scalar multiplication with respect to vector addition: $$a(\mathbf{u}+\mathbf{v})=a\mathbf{u}+a\mathbf{v}$$
> 8. Distributivity of scalar multiplication with respect to field addition: $$(a+b)\mathbf{u}=a\mathbf{u}+b\mathbf{u}$$

These axioms are necessary for a set to be considered a vector space. Examples of vector spaces include the Euclidean space $$\mathbb{R}^n$$, the set of all functions with continuous $$m$$-th derivatives $$C^m(\mathbb{R}^n)$$ or $$C^\infty(\mathbb{R}^n)$$, which is the set of all infinitely differentiable functions. In quantum mechanics we usually work in Hilbert spaces, which will be described in more detail in the next section.

Before we can use vector spaces effectively, we must identify and characterize smaller structures that inherit their properties. These are known as subspaces.

{:.definition}
> A subset $$W\subset V$$ is a subspace if it is non-empty and closed under addition and scalar multiplication. For all $$\mathbf{u},\mathbf{v}\in W$$ and $$a\in F$$, one has $$\mathbf{u}+\mathbf{v}\in W$$ and $$a\mathbf{u}\in W$$. Equivalently, $$W$$ is a subspace iff $$\mathbf{0}\in W$$, $$\mathbf{u}+\mathbf{v}\in W$$ and $$a\mathbf{u}\in W$$.

Every subspace necessarily contains the zero vector and is closed under both addition and scalar multiplication. Subspaces are important because they capture smaller, self-contained portions of a larger space where similar operations remain valid.

We now proceed to the concept of linear combinations, which provides a mechanism for generating new vectors from existing ones, and the related concept of the span, which describes the set of all such combinations.

{:.definition}
> Given a vector space $$V$$ over a field $$F$$, a linear combination of vectors $$\mathbf{v}_1,\mathbf{v}_2,\ldots,\mathbf{v}_n\in V$$ is a vector
>
> $$
> \begin{equation}
> \mathbf{u}=a_1\mathbf{v}_1+a_2\mathbf{v}_2+\ldots+a_n\mathbf{v}_n,
> \end{equation}
> $$
>
> where $$a_1,a_2,\ldots,a_n\in F$$. The set of all linear combinations of $${\mathbf{v}_1,\mathbf{v}_2,\ldots,\mathbf{v}_n}$$ is called the span of these vectors
>
> $$
> \begin{equation}
> \text{span}\lbrace\mathbf{v}_1,\mathbf{v}_2,\ldots,\mathbf{v}_n\rbrace=\left\lbrace\sum_{i=1}^n a_i\mathbf{v}_i\middle| a_i\in F\right\rbrace.
> \end{equation}
> $$
>
> The span of a set of vectors is the smallest subspace of $$V$$ that contains all the vectors in the set.

The span formalizes the idea of building new vectors from known ones. In finite-dimensional spaces, the span of a finite set can fill the entire space (as in the case of a basis) or a lower-dimensional subspace.

To understand how vectors relate to one another within a span, we must determine whether some vectors can be expressed as combinations of others. This leads us to the notion of linear independence.

{:.definition}
> A set of vectors $$\lbrace\mathbf{v}_1,\mathbf{v}_2,\ldots,\mathbf{v}_n\rbrace$$ in a vector space $$V$$ over a field $$F$$ is said to be linearly independent if the only solution to the equation
>
> $$
> \begin{equation}
> a_1\mathbf{v}_1+a_2\mathbf{v}_2+\ldots+a_n\mathbf{v}_n=\mathbf{0}
> \end{equation}
> $$
>
> is $$a_1=a_2=\ldots=a_n=0$$, where $$a_1,a_2,\ldots,a_n\in F$$. If there exists a non-trivial solution (i.e., not all $$a_i$$ are zero), then the set is said to be linearly dependent.

Linear independence ensures that no vector in the set can be constructed from others. This property is essential when defining a minimal generating set for the entire space. For linear operators (defined later), eigenvectors associated with distinct eigenvalues are linearly independent. No orthogonality is assumed here.

Finally, the concepts of basis and dimension summarize the structure of a vector space in a concise form. A basis provides coordinates for states and a matrix representation for operators.

{:.definition}
> A basis of a vector space $$V$$ over a field $$F$$ is a set of vectors $$\lbrace\mathbf{v}_1,\mathbf{v}_2,\ldots,\mathbf{v}_n\rbrace$$ in $$V$$ that is linearly independent and spans $$V$$. The dimension of $$V$$, denoted as $$\text{dim}(V)$$, is the number of vectors in any basis of $$V$$.

Every vector in a vector space can be written uniquely as a linear combination of basis vectors. The dimension quantifies the minimal number of coordinates required to specify any element of the space, which is a fundamental concept underlying state representation in quantum mechanics.

## Normed, Inner Product and Complete Vector Space

To introduce geometric concepts into vector spaces, we first define the normed vector space.

{:.definition}
> A normed vector space $$(V,\lvert\cdot\rvert)$$ is a vector space $$V$$ over a field $$F$$ (here, $$F$$ is $$\mathbb{R}$$ or $$\mathbb{C}$$) equipped with a norm, which is a function $$\lvert\cdot\rvert:V\to\mathbb{R}$$ that satisfies the following properties for all $$\mathbf{u},\mathbf{v}\in V$$ and all scalars $$a\in F$$:
>
> 1. Non-negativity: $$\lvert\mathbf{u}\rvert\geq 0$$, with equality if and only if $$\mathbf{u}=\mathbf{0}$$
> 2. Absolute scalability: $$\lvert a\mathbf{u}\rvert=\lvert a\rvert\lvert\mathbf{u}\rvert$$
> 3. Triangle inequality: $$\lvert\mathbf{u}+\mathbf{v}\rvert\leq\lvert\mathbf{u}\rvert+\lvert\mathbf{v}\rvert$$
> 4. Positive-definiteness: $$\lvert\mathbf{u}\rvert=0$$ if and only if $$\mathbf{u}=\mathbf{0}$$

A norm provides a measure of the "length" or "magnitude" of vectors in the space, allowing us to discuss concepts such as convergence and continuity. The norm induces a metric (or distance function) on the vector space defined briefly as

$$
\begin{equation}
d(\mathbf{u},\mathbf{v})=\lvert\mathbf{u}-\mathbf{v}\rvert,
\end{equation}
$$

which satisfies all the properties of a metric, thereby enabling the study of geometric properties within the vector space. More specialized normed vector space emerges after we define an inner product.

The inner product generalizes the familiar Euclidean dot product to abstract vector spaces.

{:.definition}
> An inner product on a vector space $$V$$ over a field $$F$$ is a function $$\langle\cdot,\cdot\rangle:V\times V\to F$$ that satisfies the following properties for all $$\mathbf{u},\mathbf{v},\mathbf{w}\in V$$ and all scalars $$a\in F$$:
>
> 1. Conjugate symmetry: $$\langle\mathbf{u},\mathbf{v}\rangle=\overline{\langle\mathbf{v},\mathbf{u}\rangle}$$
> 2. Linearity in the first argument: $$\langle a\mathbf{u}+\mathbf{v},\mathbf{w}\rangle=a\langle\mathbf{u},\mathbf{w}\rangle+\langle\mathbf{v},\mathbf{w}\rangle$$
> 3. Positive-definiteness: $$\langle\mathbf{u},\mathbf{u}\rangle\geq 0$$, with equality if and only if $$\mathbf{u}=\mathbf{0}$$

For real vector spaces, the conjugate symmetry reduces to ordinary symmetry. The inner product induces a norm (or length) of a vector $$\mathbf{u}$$ defined by

$$
\begin{equation}
\lvert\mathbf{u}\rvert=\sqrt{\langle\mathbf{u},\mathbf{u}\rangle}.
\end{equation}
$$

We can also generalize the concept of perpendicularity to abstract vector spaces using the inner product. Two vectors $$\mathbf{u}$$ and $$\mathbf{v}$$ are said to be orthogonal if their inner product is zero, i.e., $$\langle\mathbf{u},\mathbf{v}\rangle=0$$. If their lengths are both one, they are called orthonormal.

The inner product gives rise to two key inequalities that determine much of the geometric structure of the space. The Cauchy-Schwarz inequality states that for any vectors $$\mathbf{u},\mathbf{v}\in V$$,

$$
\begin{equation}
\lvert\langle\mathbf{u},\mathbf{v}\rangle\rvert\leq\lvert\mathbf{u}\rvert\lvert\mathbf{v}\rvert.
\end{equation}
$$

This inequality ensures that the angle defined in Equation \eqref{eq:vector_angle} is well-defined and implies the triangle inequality for the induced norm. A vector space with an inner product is called an inner product space. Together, these results ensure that the norm induced by an inner product satisfies all metric axioms, endowing the space with a well-defined notion of distance. Thus, every inner product space is a normed vector space, although the converse is not generally true.

Before we get to linear operators, we add one final definition that will become important when we discuss infinite-dimensional spaces.

{:.definition}
> Let $$(V,\lvert\cdot\rvert)$$ be a normed vector space. A sequence $$\lbrace\mathbf{v}_n\rbrace$$ in $$V$$ is called a Cauchy sequence if for every $$\epsilon>0$$, there exists an integer $$N$$ such that $$\lvert\mathbf{v}_m-\mathbf{v}_n\rvert<\epsilon$$ for all $$m,n>N$$. The space $$V$$ is said to be complete if every Cauchy sequence in $$V$$ converges to a limit that is also in $$V$$. A complete normed vector space is called a Banach space.

## Linear Operators

Linear operators lie at the core of quantum mechanics, representing physical observables and the evolution of quantum states. We now define linear operators and explore their properties.

{:.definition}
> A linear operator on avector space $$V$$ over a field $$F$$ is a map $$A:V\to V$$ with $$A(\mathbf{u}+\mathbf{v})=A\mathbf{u}+A\mathbf{v}$$ and $$A(a\mathbf{u})=aA\mathbf{u}$$ for all $$\mathbf{u},\mathbf{v}\in V$$, $$a\in F$$.

Using a basis of a vector space, we can represent linear operators as matrices. If $$\lbrace\mathbf{e}_i\rbrace_{i=1}^n$$ is a basis of an $$n$$-dimensional vector space $$V$$, an operator $$A$$ can be represented by a matrix $$[A]_{ij}=a_{ij}$$ with the coefficients defined by

$$
\begin{equation}
A\mathbf{e}_j=\sum_{i=1}^n a_{ij}\mathbf{e}_i.
\end{equation}
$$

The action of $$A$$ on any vector $$\mathbf{v}=\sum_{j=1}^n v_j\mathbf{e}_j$$ can then be computed as

$$
\begin{equation}
A\mathbf{v}=\sum_{i=1}^n\left(\sum_{j=1}^n a_{ij}v_j\right)\mathbf{e}_i.
\end{equation}
$$

If the vector space $$V$$ is equipped with an inner product and the basis $$\lbrace\mathbf{e}_i\rbrace_{i=1}^n$$ is orthonormal, the matrix elements can also be expressed as

$$
\begin{equation}
a_{ij}=\langle\mathbf{e}_i,A\mathbf{e}_j\rangle.
\end{equation}
$$

Note that the matrix representation of an operator depends on the choice of basis. If $$U$$ is a change of basis matrix from a new basis $$\lbrace\mathbf{f}_i\rbrace_{i=1}^n$$ to the old basis $$\lbrace\mathbf{e}_i\rbrace_{i=1}^n$$ (i.e., $$\mathbf{e}_i=\sum_{j=1}^n U_{ji}\mathbf{f}_j$$), the matrix representation of $$A$$ in the new basis is given by

$$
\begin{equation}
[A]_{\mathbf{e}}=U^{-1}[A]_{\mathbf{f}}U.
\end{equation}
$$

Let's now define some additional properties of linear operators.

{:.definition}
> For $$A:V\to V$$, the kernel (or null space) of $$A$$ is the set of vectors $$\mathbf{v}\in V$$ such that $$A\mathbf{v}=\mathbf{0}$$. The image (or range) of $$A$$ is the set of vectors $$\mathbf{w}\in V$$ such that $$\mathbf{w}=A\mathbf{v}$$ for some $$\mathbf{v}\in V$$. If $$\text{dim}(V)$$ is finite, the rank-nullity theorem states that $$\text{dim}(\text{ker}(A))+\text{dim}(\text{im}(A))=\text{dim}(V)$$.

Another important concept is the concept of eigenvalues and eigenvectors.

{:.definition}
> A scalar $$\lambda\in F$$ is called an eigenvalue of a linear operator $$A:V\to V$$ if there exists a non-zero vector $$\mathbf{v}\in V$$ such that $$A\mathbf{v}=\lambda\mathbf{v}$$. The vector $$\mathbf{v}$$ is called an eigenvector associated with the eigenvalue $$\lambda$$. For a finite-dimensional vector space, the set of all eigenvalues of $$A$$ is called the spectrum of $$A$$. The subspace $$E_\lambda=\text{ker}(A-\lambda I)$$ is called the eigenspace associated with the eigenvalue $$\lambda$$, where $$I$$ is the identity operator on $$V$$.

If $$V$$ is finite-dimensional, the eigenvalues of $$A$$ can be found by solving the characteristic polynomial equation $$\text{det}(A-\lambda I)=0$$. The eigenvalues may be real or complex, depending on the operator and the field $$F$$.
