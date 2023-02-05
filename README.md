<div align="center">
<h1>svSlicer</h1>

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
manager (e.g. `brew install vtk` on macOS) or build it using:

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

You can build svSlicer with the following script:

```bash
mkdir Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

## Run

```bash
OMP_NUM_THREADS=<NUM_THREADS> ./slicer <path/to/3d_result.vtu> <path/to/centerline.vtp> <path/to/output_file.vtp>
```

> Note: Parallel processing using OpenMP is not supported on macOS

## Estimated runtime

Runtime statistics for a case with 2600 slices and 200 timesteps.

| Number of threads  | Total runtime [min]  | Average runtime per slice [s][^1] | Required RAM [GB] |
| ------------------ | -------------------- |  -------------------------------- | ------------------|
| 8                  | 25                   | 0.54                              |  16               |
| 12                 | 17                   | 0.37                              |  16               |
| 16                 | 14                   | 0.29                              |  32               |
| 20                 | 12                   | 0.24                              |  32               |
| 24                 | 11                   | 0.21                              |  32               |

[^1]: *only slice calculation; time for reading and writing of files not considered*

## Development

`clang-format` is used to ensure consistent formatting throughout the code. Use
the following code to automatically apply the right format to your code changes:

```
clang-format -i --style=file:.clang_format src/main.cpp
```

This format is also enforced for pull requests to be allowed.
