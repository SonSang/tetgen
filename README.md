# Build

See build options in CMake file. Then, run following commands.

```bash
mkdir build
cd build
cmake ..
make
```

# Examples (Using Tetgen Executable)

Install following packages before running examples.

* numpy (pip install numpy)
* trimesh (pip install trimesh)

## Delaunay Tetrahedralization of random 3D points

```bash
python _random_points.py
build/tetgen random_points.node
```

## Constrained Delaunay Tetrahedralization

### Plane in Box

In this example, there is a small plane floating inside a unit cube.
We generate quality tetrahedra that conforms to the given PLC.

```bash
python _plane_in_box.py
build/tetgen -pq1.2ak _plane_in_box.poly
```

See documentation for the details of the flags.

# Examples (Using Tetgen Library)

Note that we should bind [build/libtetgen.a].
Therefore, we use "-L./build/" to specify linking directory.

## Default Example from Tetgen

Fixed the original code for using tetgenbehavior struct.

```bash
g++ -o test tetcall.cxx -L./build/ -ltetgen
```