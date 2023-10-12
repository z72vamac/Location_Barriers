The folder contains 3 subfolders:

- Barriers
- Circles

Each barrier/circle file is labeled with a number. This number refers to a different zone of the district.

To obtain a feasible problem, it is mandatory to take a circle instance with its corresponding barrier instance.

Here, an example Python code to load correctly an instance for a zone of the neighbourhood:


# Segments loading
segments = np.genfromtxt('./instances_casestudy/barriers/barriers3.csv', delimiter = ',')

barriers = []
for segment in segments:
    barriers.append([[segment[0], segment[1]], [segment[2], segment[3]]])

# Neighbourhoods loading
circles = np.genfromtxt('./instances_casestudy/circles/circles3.csv', delimiter = ',')

neighbourhoods = [Circle(center = [centerx, centery], radii = radius) for centerx, centery, radius in circles]