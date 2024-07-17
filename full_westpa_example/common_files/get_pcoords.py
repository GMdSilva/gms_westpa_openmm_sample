import mdtraj
import numpy

parent = mdtraj.load('parent.pdb', top='parent.pdb')
traj = mdtraj.load('seg.dcd', top='parent.pdb')

distances = {
    "dist1": [[1965, 4239]],
    "dist2": [[4260, 681]]
}
for name, pairs in distances.items():
    distances[name] = numpy.append(mdtraj.compute_distances(parent, pairs, periodic=True),
                                   mdtraj.compute_distances(traj, pairs, periodic=True)) * 10


dihedrals = {
    "dih1": [[4212, 4203, 4220, 4233]],
    "dih2": [[4212, 4203, 4247, 4260]]
}
for name, quadruplets in dihedrals.items():
    dihedrals[name] = (numpy.degrees(numpy.append(mdtraj.compute_dihedrals(parent, quadruplets, periodic=True),
                                                 mdtraj.compute_dihedrals(traj, quadruplets, periodic=True))) + 180)


angles = {
    "angle1": [[4150, 4233, 4260]]
}
for name, triplets in angles.items():
    angles[name] = numpy.degrees(numpy.append(mdtraj.compute_angles(parent, triplets, periodic=True),
                                              mdtraj.compute_angles(traj, triplets, periodic=True)))


for name, data in {**distances, **dihedrals, **angles}.items():
    numpy.savetxt(f"{name}.dat", data)
