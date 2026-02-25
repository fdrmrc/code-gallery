# A Discontinuous Galerkin solver for the Poisson problem on general polytopal meshes generated through mesh agglomeration


This program solves a Poisson problem on an agglomerated polytopal mesh
using a symmetric interior penalty discontinuous Galerkin (SIPG) method.
Agglomerates are constructed by an R-tree–based spatial indexing
strategy, following the approach proposed in [2]. 
In addition, a graph-based METIS partitioner is also provided in this program for comparison.



## Running the code:

As in the tutorial programs, type 

`cmake -DDEAL_II_DIR=/path/to/deal.II .` 

on the command line to configure the program. After that you can compile with `make` and run with either `make run` or using 

`./agglomeration_poisson`

on the command line. 


## Program output

Running the program produces two kinds of output. Terminal output (text summary)

The program prints a short summary including:

- the finite element degree (`FE degree`);
- the triangulation size (`Size of tria`);
- the agglomeration construction time (R-tree build time);
- the number of agglomerated subdomains (`N subdomains`);
- the number of DoFs per cell;
- the assembly time;
- a convergence table with `#DoFs`, `L2 error`, and `H1 error`.

And also the Visualization files (.vtu)



## Problem description:

We consider the Poisson problem in a bounded, simply connected domain
@f$\Omega \subset \mathbb{R}^d @f$, @f$d = 2,3@f$.
The strong formulation reads
@f{align*}
  -\Delta u &= f  && \text{in } \Omega, \\
           u &= u_D && \text{on } \partial\Omega,
@f}
where the right-hand side satisfies @f$f \in L^2(\Omega) @f$ and the prescribed
Dirichlet data satisfy @f$u_D \in H^{1/2}(\partial\Omega) @f$.

The corresponding weak formulation is: find @f$u \in H^1(\Omega)@f$ with
@f$u = u_D @f$ on @f$\partial\Omega @f$ such that
@f{equation}
  \int_{\Omega} \nabla u \cdot \nabla v \,\mathrm d\mathbf{x}
  =\int_{\Omega} f\, v \,\mathrm d\mathbf{x}
  \qquad \text{for all } v \in H_0^1(\Omega).
@f}


## SIPG discretization on agglomerated polytopal meshes:

We discretize the weak formulation by a symmetric interior penalty
discontinuous Galerkin (SIPG) method on the agglomerated polytopal mesh
@f$\mathcal T_h @f$, whose elements @f$K \in \mathcal{T}_h @f$ are mutually disjoint
open polygons (for @f$d=2 @f$) or polyhedra (for @f$d=3 @f$).
For each element we denote its diameter by
@f[
  h_K := \operatorname{diam}(K).
@f]
The mesh skeleton is given by
@f[
  \Gamma := \bigcup_{K \in \mathcal T_h} \partial K,
@f]
and we denote by @f$\Gamma_{\mathrm{int}} @f$ the union of interior faces,
while @f$\Gamma_{\mathrm D} := \Gamma \cap \partial\Omega @f$ collects the
Dirichlet boundary faces.

The discrete space @f$V_h @f$ consists of element-wise polynomials of degree
at most $p$ on each @f$K \in \mathcal T_h @f$. For @f$u_h, v_h \in V_h @f$ we use
the broken gradient @f$\nabla_h @f$ and the standard jump and average
operators @f$[\![\cdot]\!]@f$ and @f$\{\!\!\{\cdot\}\!\!\}@f$ on faces.

The DG formulation reads: find @f$u_h \in V_h @f$ such that
@f{equation}
  B(u_h,v_h) = l(v_h)
  \qquad \forall\, v_h \in V_h,
@f}
with
@f{align*}
  B(u_h,v_h)
  &=
  \int_{\Omega} \nabla_h u_h \cdot \nabla_h v_h \,\mathrm d\mathbf{x}
  \\
  &\quad
  - \int_{\Gamma}
    \Bigl(
      \{\!\!\{\nabla u_h\}\!\!\} \cdot [\![v_h]\!]
      +
      \{\!\!\{\nabla v_h\}\!\!\} \cdot [\![u_h]\!]
    \Bigr)\,\mathrm d s
  \\
  &\quad
  + \int_{\Gamma} \sigma \,[\![u_h]\!] \cdot [\![v_h]\!] \,\mathrm d s,
@f}
and
@f{equation}
  l(v_h)
  =
  \int_\Omega f\, v_h \,\mathrm d\mathbf{x}
  +
  \int_{\Gamma_{\mathrm D}}
    u_D \bigl(\sigma v_h - \nabla v_h \cdot \mathbf n\bigr)\,\mathrm d s.
@f}

The penalty parameter is chosen as
@f{equation}
  \sigma(\mathbf x) = C_\sigma
  \begin{cases}
    \dfrac{p^2}{h_K}, &
      \text{if } \mathbf x \in \partial K \cap \partial\Omega, \\[0.5em]
    \dfrac{p^2}{\min\{h_K^+,h_K^-\}}, &
      \text{if } \mathbf x \in \Gamma_{\mathrm{int}},
  \end{cases}
@f}
where @f$h_K^\pm @f$ are the diameters of the two elements sharing the
interior face, $p$ is the polynomial degree, and we fix @f$C_\sigma = 10 @f$ in this program.

This scheme is well posed and admits optimal-order a priori error
estimates. More precisely, assuming that @f$u|_K \in H^{s+1}(K) @f$ for all
@f$K \in \mathcal T_h @f$ and some @f$1 \le s \le p @f$, there exists a constant
@f$C > 0 @f$, independent of @f$h @f$, such that
@f[
  \|u - u_h\|_{L^2(\Omega)}
  \le C\, h^{s+1} \, |u|_{H^{s+1}(\Omega)},
@f]
and
@f[
  \|\nabla(u - u_h)\|_{L^2(\Omega)}
  \le C\, h^{s} \, |u|_{H^{s+1}(\Omega)}.
@f]
We refer to [2] for details of the analysis.

## Agglomeration strategies
Agglomeration is a key ingredient for constructing polytopic meshes.
This program supports two strategies for generating agglomerates,
corresponding to the choices ‘metis’ and ‘rtree’.

### R-tree-based agglomeration

This example implements an R-tree based agglomeration strategy for the construction of nested polytopic meshes. The main goal is to generate a multilevel hierarchy that preserves clear parent-child relationships, so that intergrid transfers remain simple and efficient.

In the ‘rtree’ option, axis-aligned bounding boxes of all fine cells are
inserted into a spatial R-tree. Agglomerates are obtained by grouping the
cells whose bounding boxes belong to the same node at a user-selected
extraction level of the tree. This purely geometric strategy does not
require external graph partitioners and is typically fast and scalable.
The number and shape of the agglomerates are determined by the R-tree
structure and the chosen level; see~[2] for details.


#### Basic idea and data structure

The approach is based on bounding-box packing and recursive traversal of an R-tree built from the active cells (or cut-cell bounding boxes). It is particularly suitable for multilevel methods, where a nested agglomerated structure is desirable.

Construct leaf nodes and internal nodes.

- **Leaf nodes** store the geometric objects (here, mesh cells or their bounding boxes).
- **Internal nodes** store:
  - a pointer (or reference) to a child node,
  - a bounding box that encloses all entries contained in that child subtree.
  
As a result, each internal node represents a spatial grouping of the objects below it. 
The R-tree is used here as a geometry-aware structure for grouping cell bounding boxes.  
The main design targets are:

- **Minimize box area**: reduce the area covered by each bounding box,
- **Minimize overlap**: reduce overlap between neighboring boxes,
- **Improve compactness**: reduce box perimeters (equivalently, favor more compact boxes).

These criteria improve the spatial quality of the hierarchy and typically lead to better grouping and query behavior.

Given a collection of cells (or cut-cell bounding boxes), we:

1. build an R-tree from their bounding boxes,
2. choose a target tree level,
3. recursively descend each node at that level to its leaf descendants,
4. agglomerate leaf cells that share the same ancestor.

This produces a sequence of nested agglomerated meshes, which can be used in multilevel solvers and preconditioners.


<div align="center">
  <img src="./doc/images/rtree_example.png" width="250">
  <img src="./doc/images/Tree_A_F.png" width="280">
  <br>
  <span style="display:inline-block; width:250px;"><em>(a) Examples of MBRs holding geometric data and their MBRs</em></span>
  <span style="display:inline-block; width:200px;"><em>(b) Corresponding R-tree data structure elements</em></span>
  <br>
</div>

#### Design targets and practical properties



##### Practical properties in this example

For the agglomeration workflow considered here, the R-tree-based extraction has the following practical features:

- **Level-independent extraction cost (observed)**: the wall-clock time is approximately constant with respect to the chosen extraction level.
- **Boost.Geometry backend**: the implementation relies on the `Boost.Geometry` R-tree.
- **Custom traversal logic**: hierarchy traversal in the form needed for agglomeration is not directly exposed, so a custom node visitor is implemented.

These properties make the R-tree approach convenient for constructing nested agglomerated meshes in multilevel finite element and DG settings.

The following images illustrate the R-tree-based agglomeration on a
structured fine mesh:

<div align="center">
  <img src="./doc/images/grid_raw.png" width="250">
  <img src="./doc/images/grid_raw_rtree.png" width="280">
  <br>
  <span style="display:inline-block; width:250px;"><em>(a) Original 8x8 mesh</em></span>
  <span style="display:inline-block; width:300px;"><em>(b) Bounding boxes and mesh elements</em></span>
  <br>
  <img src="./doc/images/tree_structure.png" width="550">
  <br>
  <em>(c) Tree hierarchy: root node, internal nodes, and leaf nodes from top to bottom (only two child subtrees of the root are shown for clarity)</em>
  <br>
</div>


From left to right, these plots show the original fine mesh, the blocks
induced by the R-tree on the cell bounding boxes, and the corresponding
tree structure.

#### Algorithmic pipeline

Let $\mathcal{T}_k$ denote the agglomerated mesh sequence.

- **Step 1: Build the R-tree**  
  Construct an R-tree from the set of bounding boxes associated with the fine-level cells.

- **Step 2: Select a target level**  
  Choose a tree level $l\in {1,\dots,L}$ to control the agglomeration granularity.

- **Step 3: Collect leaf descendants**  
  For each node on level $l$, recursively traverse its children until leaf nodes are reached.

- **Step 4: Agglomerate by common ancestor**  
  Merge leaf cells that belong to the same subtree (same ancestor at level \(l\)).

This yields a nested hierarchy with natural parent-child relations across levels.

Using an R-tree in the agglomeration pipeline has two practical advantages:

- it provides a natural hierarchy over the fine-level mesh entities,
- it enables geometry-aware grouping through recursive subtree extraction.

This is the key mechanism used later to build nested agglomerated meshes for multilevel methods and and DG frameworks.


### METIS-based partitioning 
In the ‘metis’ option, the adjacency graph of the fine mesh is
constructed with one vertex per cell and edges between face-neighbouring
cells. This graph is then partitioned by the multilevel graph partitioner
METIS into a prescribed number of parts, and each part defines one
agglomerate; see~[3] for details.


## Test case:

We consider the Poisson problem on the unit square @f$\Omega = (0,1)^2 @f$
with the manufactured exact solution
@f[
  u(x,y) = \sin(\pi x)\sin(\pi y).
@f]
The corresponding right-hand side is
@f[
  f(x,y) = 2\pi^2 \sin(\pi x)\sin(\pi y).
@f]
This manufactured solution allows us to compute the global
@f$L^2 @f$- and @f$H^1 @f$-seminorm errors of the discrete solution in order to
assess the quality of the numerical approximation.

In this example, an unstructured fine mesh (e.g., a triangular mesh) is
used as the starting point. Agglomerates are then constructed by METIS
and by the R-tree strategy, leading to different polytopal meshes. The
following images compare the resulting agglomerates for two different
numbers of agglomerates:

<h4>Comparison of agglomeration strategies</h4>

<div align="center" style="white-space:nowrap;">
  <img src="./doc/images/polygonmetis_91.png"
       style="width:240px; display:inline-block; vertical-align:top; margin-right:18px;">
  <img src="./doc/images/polygonmetis_364.png"
       style="width:240px; display:inline-block; vertical-align:top;">
  <br>
  <span style="display:inline-block; width:200px;"><em>(a) METIS, n_partitions = 91</em></span>
  <span style="display:inline-block; width:300px;"><em>(b) METIS, n_partitions = 364</em></span>
  <br>
</div>  

<div align="center" style="white-space:nowrap;">
  <img src="./doc/images/polygonrtree_91.png"
       style="width:240px; display:inline-block; vertical-align:top; margin-right:18px;">
  <img src="./doc/images/polygonrtree_364.png"
       style="width:240px; display:inline-block; vertical-align:top;">
  <br>
  <span style="display:inline-block; width:200px;"><em>(c) R-tree, extraction_level = 4</em></span>
  <span style="display:inline-block; width:300px;"><em>(d) R-tree, extraction_level = 5</em></span>
  <br>
</div>

These plots illustrate how the two strategies distribute and shape the
agglomerates on the same underlying unstructured mesh.



The corresponding error curves are shown below:
<div align="center">
  <img src="./doc/images/test_result.png" width="630">
  <br>
  <span style="display:inline-block; width:200px;"><em>(c) R-tree, extraction_level = 4</em></span>
</div>
The figure reports the @f$L^2@f$- and @f$H^1 @f$-seminorm errors with respect to
the manufactured solution @f$u @f$. Optimal convergence rates are observed
for all polynomial degrees and for both agglomeration strategies. In
addition, the curves associated with the R-tree approach are consistently
lower than or comparable to those obtained with METIS-based partitioning.

The construction of an R-tree spatial index on an arbitrary fine grid
provides a natural and efficient agglomeration strategy with the
following features:

- fully automated, robust, and dimension-independent;
- it produces a balanced and nested hierarchy of agglomerates;
- the shape of the agglomerates closely follows their axis-aligned
  bounding boxes.

These properties make the R-tree approach an attractive alternative to
graph-based agglomeration methods (see [3] for more details).



## References 
* [1] Di Pietro, Daniele Antonio and Ern, Alexandre (2012), Mathematical Aspects of Discontinuous Galerkin Methods. ISBN: [978-3-642-22980-0](https://www.springer.com/gp/book/9783642229794)
* [2] Marco Feder, Andrea Cangiani and Luca Heltai (2025), R3MG: R-tree based agglomeration of polytopal grids with applications to multilevel methods. DOI: [10.1016/j.jcp.2025.113773](https://doi.org/10.1016/j.jcp.2025.113773)
* [3] George Karypis and Vipin Kumar, A fast and high quality multilevel scheme for partitioning irregular graphs. DOI: [10.1137/S1064827595287997](https://doi.org/10.1137/S1064827595287997)

