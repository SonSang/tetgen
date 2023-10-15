import trimesh
import argparse
import os
import time

import numpy as np
import matplotlib.pyplot as plt

import pytetgen

parser = argparse.ArgumentParser(description='Tetrahedralize an arbitrary surface mesh.')
parser.add_argument('--mesh', type=str, default=None, help='Path to the target mesh.')
parser.add_argument('--quiet', action='store_true', help='Suppress tetgen output.')
parser.add_argument('--output', type=str, default='logdir', help='Path to the output.')
args = parser.parse_args()

logdir = os.path.join(args.output, time.strftime("%Y-%m-%d-%H-%M-%S"))
if not os.path.exists(logdir):
    os.makedirs(logdir)

if args.mesh is None:
    # generate a sphere;
    print("No mesh input, generating a sphere.")
    mesh = trimesh.creation.icosphere(subdivision=1)
else:
    mesh = trimesh.load(args.mesh)

# get the bounding box;
min_bound = mesh.bounds[0]
max_bound = mesh.bounds[1]

# define a box that contains the mesh;
box = trimesh.creation.box(
    extents=(max_bound - min_bound) * 1.1,
    transform=trimesh.transformations.translation_matrix((min_bound + max_bound) * 0.5)
)

# define entire mesh;
box_vertices = box.vertices
box_faces = box.faces
mesh_vertices = mesh.vertices
mesh_faces = mesh.faces

vertices = np.concatenate([box_vertices, mesh_vertices], axis=0)
faces = np.concatenate([box_faces, mesh_faces + len(box_vertices)], axis=0)
faces_marker = [1] * len(box_faces) + [2] * len(mesh_faces)

# define tetgenio;
in_mesh = pytetgen.tetgenio()
in_mesh.firstnumber = 0
in_mesh.set_pointlist(vertices)
in_mesh.set_facetlist_from_trimesh(faces)
in_mesh.set_facetmarkerlist(faces_marker)

out_mesh = pytetgen.tetgenio()

# call tetrahedralize;
command = "pq1.414a0.1"
if args.quiet:
    command += "Q"
b = pytetgen.tetgenbehavior()
b.parse_commandline(command)
pytetgen.tetrahedralize(b, in_mesh, out_mesh)

# save the output;
out_verts = out_mesh.get_pointlist()

out_boundary_faces = out_mesh.get_trifacelist()
out_boundary_face_markers = out_mesh.get_trifacemarkerlist()
out_boundary_faces_with_markers = np.concatenate([out_boundary_faces, 
                                                out_boundary_face_markers[:, None]], axis=1)

out_tetra = out_mesh.get_tetrahedronlist()
out_tetra_faces = np.concatenate([out_tetra[:, [0, 1, 2]], 
                                out_tetra[:, [0, 1, 3]],
                                out_tetra[:, [0, 2, 3]],
                                out_tetra[:, [1, 2, 3]]], axis=0)
out_tetra_faces = np.sort(out_tetra_faces, axis=1)
out_tetra_faces = np.unique(out_tetra_faces, axis=0)
out_tetra_face_markers = np.zeros(len(out_tetra_faces), dtype=np.int32)
out_tetra_faces_with_markers = np.concatenate([out_tetra_faces,
                                            out_tetra_face_markers[:, None]], axis=1)

out_faces = np.concatenate([out_boundary_faces_with_markers, out_tetra_faces_with_markers], axis=0)
out_faces = np.unique(out_faces, axis=0)
u_out_faces, u_out_faces_counts = np.unique(out_faces[:, [0, 1, 2]], return_counts=True, axis=0)
u_out_faces_counts_cumsum = np.cumsum(u_out_faces_counts)
out_faces, out_faces_marker = out_faces[u_out_faces_counts_cumsum - 1][:, [0, 1, 2]], out_faces[u_out_faces_counts_cumsum - 1][:, 3]



# render the output;
out_trimesh = trimesh.Trimesh(out_verts, out_faces)
sweep_start = min_bound
sweep_end = max_bound
sweep_count = 200
sweep_direction = np.array([0, 0, 1])
sweep_direction = sweep_direction / np.linalg.norm(sweep_direction)
sweep_right_dir = np.array([1, 0, 0])
sweep_right_dir = sweep_right_dir / np.linalg.norm(sweep_right_dir)
sweep_up_dir = np.cross(sweep_direction, sweep_right_dir)
sweep_up_dir = sweep_up_dir / np.linalg.norm(sweep_up_dir)
sweep_distance = np.dot((sweep_end - sweep_start) / sweep_count, sweep_direction)

render_domain_right_max = np.abs(np.dot(max_bound - sweep_start, sweep_right_dir)).max()
render_domain_up_max = np.abs(np.dot(max_bound - sweep_start, sweep_up_dir)).max()

for i in range(sweep_count):
    # find cross section of the mesh;
    plane_origin = sweep_start + sweep_direction * (sweep_distance * i)
    plane_normal = sweep_direction
    
    cross_section = trimesh.intersections.mesh_plane(
        out_trimesh, plane_normal, plane_origin, return_faces=True
    )

    lines = cross_section[0]
    face_index = cross_section[1]

    # render;
    plt.figure(figsize=(8, 8))
    # plt.xlim(-render_domain_right_max, render_domain_right_max)
    # plt.ylim(-render_domain_up_max, render_domain_up_max)
    # plt.axis("off")
    
    for j in range(len(lines)):
        line = lines[j]
        line_start = line[0]
        line_end = line[1]
        line_start_x = np.dot(line_start - plane_origin, sweep_right_dir)
        line_start_y = np.dot(line_start - plane_origin, sweep_up_dir)
        line_end_x = np.dot(line_end - plane_origin, sweep_right_dir)
        line_end_y = np.dot(line_end - plane_origin, sweep_up_dir)
        color = "red" if out_faces_marker[face_index[j]] == 2 else "black"
        plt.plot([line_start_x, line_end_x], [line_start_y, line_end_y], color=color, linewidth=1.0)
    plt.savefig(f"{logdir}/cross_section_{i}.png", dpi=300)
    plt.close()
