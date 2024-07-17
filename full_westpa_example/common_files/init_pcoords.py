import mdtraj
import numpy

parent = mdtraj.load('basis.pdb', top='basis.pdb')

distances = {
    "dist1": [[1965, 4239]],
    "dist2": [[4260, 681]]
}
angles = {
    "angle1": [[4150, 4233, 4260]]
}

for name, pairs in distances.items():
    distances[name] = numpy.asarray(mdtraj.compute_distances(parent, pairs, periodic=True)) * 10

for name, triplets in angles.items():
    angles[name] = numpy.degrees(numpy.asarray(mdtraj.compute_angles(parent, triplets, periodic=True)))

for name, data in {**distances, **angles}.items():
    numpy.savetxt(f"{name}.dat", data)
