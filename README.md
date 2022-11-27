# Slice extractor

Modified version of slice extractor from
https://github.com/lucapegolotti/cardiovascular/tree/master/slice-vtu-results-centerlines/c++

## Build

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
