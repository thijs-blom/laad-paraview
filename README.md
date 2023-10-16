# ParaView LAAD Filter
This repository contains an implementation of the locally-averaged angular deviation (LAAD) filter for use in [ParaView](https://www.paraview.org/).
The plug-in is written in C++, and may be built using the [ParaViewEasyPluginBuilder](https://gitlab.kitware.com/paraview/paraview-easy-plugin-builder).
Please make sure that the ParaView version used to build the plug-in (exactly) matches your version of ParaView.

The computation is parallelised using OpenMP.
To explicitly specify the number of threads used, set the [OMP_NUM_THREADS](https://www.openmp.org/spec-html/5.0/openmpse50.html) environment variable.

# About the LAAD Filter
The LAAD filter is a measurement for the average local deviation in flow direction at a fixed point in time.
At each point, the filter computes the value
$$f_s(\mathbf{x}) = \frac{1}{N\cdot \pi}\sum_{\mathbf{y} \in B_s(\mathbf{x})} \arccos{\frac{\langle \mathbf{v}, \mathbf{u}(\mathbf{y}) \rangle}{\|\mathbf{v}\|\|\mathbf{u}(\mathbf{y})\|}}$$
for each point $\mathbf{x}$ in the mesh, where

- $\mathbf{u}$ is the (three-dimensional) velocity field,
- $s$ is a user-supplied scale,
- $B_s(\mathbf{x})$ is a sphere around $\mathbf{x}$ with radius $s$,
- $N$ is the number of points in the mesh that lie in $B_s(\mathbf{x})$,
- $\mathbf{v}$ is the average vector in the sphere, given by $\mathbf{v} = \frac{1}{N} \sum_{y \in B_s(\mathbf{x})} \mathbf{u}(\mathbf{y})$,
- $\langle \cdot, \cdot \rangle$ is the standard inner product in $\mathbb{R}^3$.

# Limitations
The filter currently requires a [vtkUnstructuredGrid](https://vtk.org/doc/nightly/html/classvtkUnstructuredGrid.html) as input.
Support for more general input types, such as a [vtkDataSet](https://vtk.org/doc/nightly/html/classvtkDataSet.html) or [vtkMultiBlockDataSet](https://vtk.org/doc/nightly/html/classvtkMultiBlockDataSet.html) may be added in a future release.

The filter acts on a three-dimensional velocity field.
It is assumed that this is specified per component as one-dimensional datasets named "u", "v", and "w".
The naming of these datasets is important, and is not user-changeable at this time.
