The folder contains 3 subfolders:

- Barriers
- Circles
- Segments

Each barrier/circle/segment file is labeled with two numbers. The first refers to the number of neighbourhoods to be visited and the second, one of the five instances that are generated for each number of neighbourhoods.

To obtain a feasible problem, it is mandatory to take a circle/segment instance with its corresponding barrier instance.

Here, an example Python code to load correctly an instance with 10 circles to visit is provided:


# Segments loading
segments = np.genfromtxt('./instances_random/barriers/barriers10-2.csv', delimiter = ',')

barriers = []
for segment in segments:
    barriers.append([[segment[0], segment[1]], [segment[2], segment[3]]])

# Neighbourhoods loading
circles = np.genfromtxt('./instances_random/circles/circles10-2.csv', delimiter = ',')

neighbourhoods = [Circle(center = [centerx, centery], radii = radius) for centerx, centery, radius in circles]