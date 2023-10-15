with open("_plane_in_box.poly", "w") as f:

    # 1. Node list;
    positions = [
        [0, 0, 0],
        [0, 1, 0],
        [0, 1, 1],
        [0, 0, 1],
        [1, 0, 0],
        [1, 1, 0],
        [1, 1, 1],
        [1, 0, 1],
        # floating face in the box;
        [0.25, 0.25, 0.5],
        [0.25, 0.75, 0.5],
        [0.75, 0.75, 0.5],
        [0.75, 0.25, 0.5],
    ]
    num_points = len(positions)
    dimension = 3
    num_attrs = 0
    use_bdry = 1
    f.write("# 1. Node list\n")
    f.write("{} {} {} {}\n".format(num_points, dimension, num_attrs, use_bdry))
    for i, position in enumerate(positions):
        f.write("{} {} {} {} {}\n".format(i + 1, position[0], position[1], position[2], 1 if i < 8 else 2))
    f.write("\n")

    # 2. Facet list;
    f.write("# 2. Facet list\n")
    f.write("7 1\n")          # 7 facets, 1 boundary marker
    f.write("1 0 1\n")        # facet 1, no holes, boundary marker 1
    f.write("4 1 2 3 4\n")    # facet 1, 4 nodes
    f.write("1 0 1\n")        # facet 2, no holes, boundary marker 1
    f.write("4 5 6 7 8\n")    # facet 2, 4 nodes
    f.write("1 0 1\n")        # facet 3, no holes, boundary marker 1
    f.write("4 1 5 6 2\n")    # facet 3, 4 nodes
    f.write("1 0 1\n")        # facet 4, no holes, boundary marker 1
    f.write("4 2 6 7 3\n")    # facet 4, 4 nodes
    f.write("1 0 1\n")        # facet 5, no holes, boundary marker 1
    f.write("4 3 7 8 4\n")    # facet 5, 4 nodes
    f.write("1 0 1\n")        # facet 6, no holes, boundary marker 1
    f.write("4 1 5 8 4\n")    # facet 6, 4 nodes
    f.write("1 0 2\n")        # facet 7, no holes, boundary marker 2
    f.write("4 9 10 11 12\n") # facet 7, 4 nodes
    f.write("\n")

    # 3. Hole list;
    f.write("# 3. Hole list\n")
    f.write("0\n")            # no holes
    f.write("\n")

    # 4. Region list;
    f.write("# 4. Region list\n")
    f.write("0\n")            # no regions

