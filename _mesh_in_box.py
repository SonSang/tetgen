import trimesh
import argparse

import numpy as np

parser = argparse.ArgumentParser(description='Generate a mesh in a box.')
parser.add_argument('--mesh', type=str, default=None, help='Path to the target mesh.')
args = parser.parse_args()

if args.mesh is None:
    # Generate a sphere;
    print("No mesh input, generating a sphere.")
    mesh = trimesh.creation.icosphere(subdivision=1)
else:
    mesh = trimesh.load(args.mesh)

# Get the bounding box;
min_bound = mesh.bounds[0]
max_bound = mesh.bounds[1]

# Define a box that contains the mesh;
box = trimesh.creation.box(
    extents=(max_bound - min_bound) * 1.1,
    transform=trimesh.transformations.translation_matrix((min_bound + max_bound) * 0.5)
)

# Save to .poly;
with open("_mesh_in_box.poly", 'w') as f:
    
    num_box_points = len(box.vertices)
    num_mesh_points = len(mesh.vertices)

    index_start = 0

    # 1. Node list;
    positions = np.concatenate([box.vertices, mesh.vertices], axis=0)
    num_points = len(positions)

    dimension = 3
    num_attrs = 0
    use_bdry = 1
    f.write("# 1. Node list\n")
    f.write("{} {} {} {}\n".format(num_points, dimension, num_attrs, use_bdry))
    for i, position in enumerate(positions):
        f.write("{} {} {} {} {}\n".format(i + index_start, 
                                        position[0], 
                                        position[1], 
                                        position[2], 
                                        1 if i < num_box_points else 2))
    f.write("\n")

    # 2. Facet list;
    num_facets = len(box.faces) + len(mesh.faces)
    f.write("# 2. Facet list\n")
    f.write(f"{num_facets} 1\n")            # num. facets, 1 boundary marker
    for i in range(len(box.faces)):
        facet_nodes = box.faces[i] + index_start
        f.write("1 0 1\n")                  # 1 polygon, no holes, boundary marker 1
        f.write(f"3 {facet_nodes[0]} {facet_nodes[1]} {facet_nodes[2]}\n")   # 3 nodes
    
    for i in range(len(mesh.faces)):
        facet_nodes = mesh.faces[i] + num_box_points + index_start
        f.write("1 0 2\n")                  # 1 polygon, no holes, boundary marker 2
        f.write(f"3 {facet_nodes[0]} {facet_nodes[1]} {facet_nodes[2]}\n")    # 3 nodes
    f.write("\n")

    # 3. Hole list;
    f.write("# 3. Hole list\n")
    f.write("0\n")            # no holes
    f.write("\n")

    # 4. Region list;
    f.write("# 4. Region list\n")
    f.write("0\n")            # no regions
