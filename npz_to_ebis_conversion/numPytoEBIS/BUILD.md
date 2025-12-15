# Building npz_to_ebis_cpp

This document describes how to compile the `npz_to_ebis_cpp` program.

## Prerequisites

### 1. Build Chombo Libraries

The program requires Chombo 3D libraries with Embedded Boundary support. From `Chombo/lib/`:

```bash
cd ../../Chombo/lib
make all DIM=3 DEBUG=FALSE USE_EB=TRUE
```

This builds the required libraries:
- `libworkshop3d.Linux.64.g++.gfortran.OPT.a`
- `libebamrelliptic3d.Linux.64.g++.gfortran.OPT.a`
- `libebamrtools3d.Linux.64.g++.gfortran.OPT.a`
- `libebtools3d.Linux.64.g++.gfortran.OPT.a`
- `libamrelliptic3d.Linux.64.g++.gfortran.OPT.a`
- `libamrtools3d.Linux.64.g++.gfortran.OPT.a`
- `libboxtools3d.Linux.64.g++.gfortran.OPT.a`
- `libbasetools3d.Linux.64.g++.gfortran.OPT.a`

## Manual Compilation Method

The GNUmakefile has a bug that causes the object file to be linked twice. Use this manual method:

### Step 1: Compile Object File

```bash
cd /home/rochi/work/personal_repos/dan_martin_project/perlmutter_files/geometry/test_build
make clean DIM=3 DEBUG=FALSE
```

This creates: `o/3d.Linux.64.g++.gfortran.OPT/npz_to_ebis_cpp.o`

### Step 2: Link Executable Manually

```bash
g++ -O3 -std=c++11 \
  o/3d.Linux.64.g++.gfortran.OPT/npz_to_ebis_cpp.o \
  -lgfortran -lm -llapack -lblas \
  -L../../Chombo/lib \
  -lworkshop3d.Linux.64.g++.gfortran.OPT \
  -lebamrelliptic3d.Linux.64.g++.gfortran.OPT \
  -lebamrtools3d.Linux.64.g++.gfortran.OPT \
  -lebtools3d.Linux.64.g++.gfortran.OPT \
  -lamrelliptic3d.Linux.64.g++.gfortran.OPT \
  -lamrtools3d.Linux.64.g++.gfortran.OPT \
  -lboxtools3d.Linux.64.g++.gfortran.OPT \
  -lbasetools3d.Linux.64.g++.gfortran.OPT \
  -L/home/rochi/installation_files/hdf5_110_install/lib \
  -lhdf5 -lhdf5_hl -lz \
  -lgfortran -lm \
  -o npz_to_ebis_cpp3d.Linux.64.g++.gfortran.OPT.ex
```

### Output

Executable: `npz_to_ebis_cpp3d.Linux.64.g++.gfortran.OPT.ex`

Verify the build:
```bash
ls -lh npz_to_ebis_cpp3d.Linux.64.g++.gfortran.OPT.ex
file npz_to_ebis_cpp3d.Linux.64.g++.gfortran.OPT.ex
```

## Debug Build

For a debug build, change the flags:

### Step 1: Compile
```bash
make clean DIM=3 DEBUG=TRUE
```

### Step 2: Link
```bash
g++ -g -std=c++11 \
  o/3d.Linux.64.g++.gfortran.DEBUG/npz_to_ebis_cpp.o \
  -lgfortran -lm -llapack -lblas \
  -L../../Chombo/lib \
  -lworkshop3d.Linux.64.g++.gfortran.DEBUG \
  -lebamrelliptic3d.Linux.64.g++.gfortran.DEBUG \
  -lebamrtools3d.Linux.64.g++.gfortran.DEBUG \
  -lebtools3d.Linux.64.g++.gfortran.DEBUG \
  -lamrelliptic3d.Linux.64.g++.gfortran.DEBUG \
  -lamrtools3d.Linux.64.g++.gfortran.DEBUG \
  -lboxtools3d.Linux.64.g++.gfortran.DEBUG \
  -lbasetools3d.Linux.64.g++.gfortran.DEBUG \
  -L/home/rochi/installation_files/hdf5_110_install/lib \
  -lhdf5 -lhdf5_hl -lz \
  -lgfortran -lm \
  -o npz_to_ebis_cpp3d.Linux.64.g++.gfortran.DEBUG.ex
```

Note: Debug libraries must be built first with `DEBUG=TRUE` in the Chombo build step.

## Known Issues

- **GNUmakefile Bug**: The makefile's link command duplicates the object file, causing "multiple definition" errors. This is why manual linking is required.
- **HDF5 Path**: The HDF5 path is hardcoded to `/home/rochi/installation_files/hdf5_110_install/lib`. Adjust if your HDF5 installation is elsewhere.

## Running the Executable

### Usage

```bash
./npz_to_ebis_cpp3d.Linux.64.g++.gfortran.OPT.ex <input.npz> <output.ebis.hdf5> [level_value] [inside_regular]
```

### Arguments

1. **`<input.npz>`** (required): Binary input file containing SDF (Signed Distance Field) data
   - Format: Binary file with header (4 ints: res_x, res_y, res_z, data_size) followed by double-precision data
   - Must be cubic data (res_x == res_y == res_z)

2. **`<output.ebis.hdf5>`** (required): Output HDF5 file for the generated EBIS (Embedded Boundary Index Space)

3. **`[level_value]`** (optional, default: 0.0): The iso-surface level value
   - The zero level-set (surface) is defined where SDF == level_value
   - Typical value: 0.0

4. **`[inside_regular]`** (optional, default: 0): Whether inside or outside is "regular"
   - 0 (false): Outside the surface is regular
   - 1 (true): Inside the surface is regular

### Example

```bash
# Using example geometry data (if available)
./npz_to_ebis_cpp3d.Linux.64.g++.gfortran.OPT.ex ../../npz_files/sphere.npz sphere.ebis.hdf5 0.0 0
```

### Runtime Requirements

**Shared Libraries**:
The executable requires these shared libraries at runtime:
- **HDF5**: `libhdf5.so` and `libhdf5_hl.so`
- **Fortran**: `libgfortran.so`
- **Math**: Standard C math library
- **LAPACK/BLAS**: May be statically linked or require `liblapack.so`, `libblas.so`

**Check Dependencies**:
```bash
ldd npz_to_ebis_cpp3d.Linux.64.g++.gfortran.OPT.ex
```

**Environment Variables**:
If libraries are not in standard paths, set:
```bash
export LD_LIBRARY_PATH=/home/rochi/installation_files/hdf5_110_install/lib:$LD_LIBRARY_PATH
```

### Input File Format

The `.npz` file is actually a **binary format** (not NumPy's compressed archive) with this structure:

```
[Header: 4 integers]
- res_x (int32)
- res_y (int32)
- res_z (int32)
- data_size (int32)

[Data: data_size doubles]
- SDF values (float64/double)
```

Data ordering: `index = i + j*res_x + k*res_x*res_y` (C-order, flattened 3D array)

### Output

The program generates:
- **EBIS HDF5 file**: Contains the Embedded Boundary Index Space that Chombo can use
- **Console output**: Resolution, domain info, data statistics, timing information

### Performance Notes

- Memory usage scales with resolution³
- Uses MPI if compiled with `MPI=TRUE` (this build is serial)
- Includes timing instrumentation via Chombo's CH_TIMER system

## Troubleshooting

If linking fails with "undefined reference to LAPACK functions":
- Ensure LAPACK/BLAS libraries are installed: `sudo apt-get install libblas-dev liblapack-dev`
- Check that `-llapack -lblas` are in the link command

If linking fails with "cannot find -lhdf5":
- Verify HDF5 installation path
- Update the `-L/path/to/hdf5/lib` in the link command

If runtime fails with "cannot open shared object file":
- Check library paths: `ldd ./npz_to_ebis_cpp3d.Linux.64.g++.gfortran.OPT.ex`
- Set `LD_LIBRARY_PATH` to include HDF5 library directory
