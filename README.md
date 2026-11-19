<div align="center">
<h1>svSlicer</h1>
</div>

svSlicer is a tool for mapping 3D simulation results from SimVascular
software to a centerline. The quantities of interest, namely "pressure" and
"velocity" area extracted from the 3D results, sliced at each centerline point,
and integrated over the respective area of the slice to derive scaler quantities
"pressure" and "flow".

> Modified version of slice extractor from
https://github.com/lucapegolotti/cardiovascular/tree/master/slice-vtu-results-centerlines/c++

## Build

### Dependencies

VTK is required for svSlicer. You can either install it with a package
manager (e.g. `brew install vtk` on macOS, `apt install libvtk9-dev` on Debian) or build it using:

```
mkdir -p ./vtk/src
git clone --recursive https://gitlab.kitware.com/vtk/vtk.git ./vtk/src
mkdir -p ./vtk/build
cmake -DBUILD_SHARED_LIBS:BOOL=ON \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_COMPILER=/share/software/user/open/gcc/12.1.0/bin/g++ \
   -DCMAKE_C_COMPILER=/share/software/user/open/gcc/12.1.0/bin/gcc \
   -S ./vtk/src -B ./vtk/build/
cmake --build ./vtk/build/ -- -j4
```

### Build svSlicer

You can build svSlicer with the following commands:

```bash
mkdir Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cmake --install .
```

> For macOS, the OpenMP runtime is required. It can be installed using
> `brew install libomp` (the build picks it up automatically via Homebrew).

> The `libvtk9-dev` package of Ubuntu 24.04 cannot read VTK XML files with raw
> appended data since a libexpat security update
> ([Launchpad bug 2058415](https://bugs.launchpad.net/ubuntu/+source/paraview/+bug/2058415)).
> Use Debian's `libvtk9-dev` (bookworm or later), Homebrew, or a VTK built from
> source instead.

## Run

```bash
OMP_NUM_THREADS=<NUM_THREADS> svslicer <path/to/3d_result.vtu> <path/to/centerline.vtp> <path/to/output_file.vtp>
```

Velocity point arrays (names starting with `velocity`, ignoring case) must
have at least three components; arrays with fewer components cause an error.
When an array has more than three components, the first three are used.

Each participating OpenMP thread slices its own copy of the 3D result, so
memory use grows by about one result mesh per additional thread. The thread
count is capped by the number of centerline points, and grid copies are
allocated only for the team actually created by OpenMP.

## Test

After building, run the regression test from the build directory:

```bash
ctest --output-on-failure
```

It generates a synthetic tube with linear velocity and pressure fields, runs
svSlicer single-threaded and five times with 8 threads, and checks that every
output matches the analytic flow and pressure and agrees with the
single-threaded output to 1e-12 relative. Additional checks cover invalid velocity component
counts, four-component velocities, single-point centerlines, and OpenMP thread
limits. Use `-DSVSLICER_TEST_THREADS=<n>` and
`-DSVSLICER_TEST_REPEATS=<n>` to change the thread count and repeats, or
`-DSVSLICER_BUILD_TESTS=OFF` to skip building the test.

## Estimated runtime

Runtime statistics for a case with 3619 slices of a 3D mesh with 1.6 million
cells for 200 timesteps using two 12-core Intel Xeon Gold 5118 CPUs.

| Number of parallel processes  | Total runtime [s] | Slices per second [1/s] |
| ----------------------------- | ----------------- |  ---------------------- |
| 1                             | 1231              | 3.0                         |
| 2                             | 653               | 5.6                         |
| 4                             | 373               | 9.9                         |
| 8                             | 229               | 16.4                             |
| 12                            | 191               | 19.7                             |
| 16                            | 179               | 21.1                             |
| 20                            | 183               | 21.5                             |
| 24                            | 185               | 20.5                             |

## Development

`clang-format` is used to ensure consistent formatting throughout the code. The
CI pins version 20.1.7 (`pipx install clang-format==20.1.7`); use the same
version locally to avoid formatting differences between clang-format releases.
Apply the right format to your code changes with:

```
clang-format -i --style=file:.clang_format src/main.cpp tests/*.cpp
```

This format is also enforced for pull requests to be allowed.
