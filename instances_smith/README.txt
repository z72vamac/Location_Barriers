The folder contains 2 subfolders:

- Barriers
- Circles

Each barrier/circle/segment file is labeled with three numbers. The first refers to the number of neighbourhoods to be visited. The second, one of the ten instances that are generated for each number of neighbourhoods. The third label represents the radii of the circles used in the experiment.

To obtain a feasible problem, it is mandatory to take a circle/segment instance with its corresponding barrier instance.

Here, an example Python code to load correctly an Smith instance with 15 circles of radii 0.5 is provided:

# Segments loading
segments = np.genfromtxt('./instances_smith/barriers/barriers15-0-0.5.csv', delimiter = ',')

barriers = []
for segment in segments:
    barriers.append([[segment[0], segment[1]], [segment[2], segment[3]]])

# Neighbourhoods loading
circles = np.genfromtxt('./instances_smith/circles/circles15-0-0.5.csv', delimiter = ',')

neighbourhoods = [Circle(center = [centerx, centery], radii = radius) for centerx, centery, radius in circles]