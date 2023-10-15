# Original README

This is TetGen version 1.6.0 (released on August 31, 2020)

Please see the documentation of TetGen for compiling and using TetGen.
It is available at the following link:

            http://www.tetgen.org

For more information on this product, contact :

  Hang Si
  Research Group of Numerical Mathematics and Scientific Computing
  Weierstrass Institute for Applied Analysis and Stochastics
  Mohrenstr. 39
  10117 Berlin, Germany

 EMail: <si@wias-berlin.de>
 Web Site: http://www.wias-berlin.de/~si

------------------- IMPORTANCE NOTICE -----------------------------

BEFORE INTALLING OR USING TetGen(R) READ the 
GENERAL LICENSE TERMS AND CONDITIONS

-------------------------------------------------------------------

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
python examples/exectuable_based/_random_points.py
build/tetgen random_points.node
```

## Constrained Delaunay Tetrahedralization

### Plane in Box

In this example, there is a small plane floating inside a unit cube.
We generate quality tetrahedra that conforms to the given PLC.

```bash
python examples/exectuable_based/_plane_in_box.py
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

# Examples (Using pytetgen)

We need following packages to run these examples.

* pybind11 

```bash
git submodule add -b stable ../../pybind/pybind11 extern/pybind11
git submodule update --init
```

Run following command to build pytetgen.

```bash
pip install -e .
```

Note that -e flag should be used.

* matplotlib (pip install matplotlib)

## Delaunay Tetrahedralization of 3D mesh



```bash
