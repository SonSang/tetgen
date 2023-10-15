import numpy as np

NUM_POINTS = 1_000_000

np.random.seed(1)

# Generate 1000 random points in 3-D
points = np.random.random((NUM_POINTS, 3))

# Save;
with open("random_points.node", 'w') as f:
    f.write("{} 3 0 0\n".format(len(points)))
    for i, point in enumerate(points):
        f.write("{} {} {} {}\n".format(i+1, point[0], point[1], point[2]))