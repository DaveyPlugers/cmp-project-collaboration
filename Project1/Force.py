import numpy as np

def DistancePoints(Vector1,Vector2,BoxSize):
    """
    Takes 2 position vectors and calculates distance between them
    :param Vector1: position vector for particle one
    :param Vector2: position vector for particle two
    :return: distance between two particles
    """
    Distance = sum([x * x for x in ((Vector1 - Vector2 + BoxSize/2)%BoxSize - BoxSize/2)])**0.5
    return Distance



def VecPairLJForce(Vector1,Vector2,BoxSize):
    """
    Takes 2 position vectors and returns the Lennard-Jones force on Vector1 due to Vector2
    :param Vector1: position vector for particle one
    :param Vector2: position vector for particle two
    :return: Lennard-Jones force vector between two particles
    """
    Dist = DistancePoints(Vector1,Vector2, BoxSize)
    Constant = 24*(2/Dist**14 - 1/Dist**8)  #Extra 1/r for normalisation Vector1 - Vector2 in next line
    Force = [Constant*x for x in ((Vector1 - Vector2 + BoxSize/2)%BoxSize - BoxSize/2)]
    return Force


def TotalForce(j,tstep,Dimension,ParticleAmount,AllPositions,BoxSize, Particles):
    """
    This calculates the total force vector experienced by particle j
    :param j: index of a particle
    :param tstep: time step with a special case for tstep==-1 for the updated values that aren't saved yet
    :return: Total force on particle j
    """
    TotalForce = np.array([0 for l in range(Dimension)], dtype=np.float64)
    if tstep ==-1: #This is a workaround to get the updated position and thus updated force
        for i in range(ParticleAmount):
            if i !=j:
                TotalForce += VecPairLJForce(Particles[j,0], Particles[i,0], BoxSize)
    else:
        for i in range(ParticleAmount):
            if i !=j:
                TotalForce += VecPairLJForce(AllPositions[tstep][j, 0], AllPositions[tstep][i][0], BoxSize)    #We use the Allpositions instead of Particles since we need the old positions

    return TotalForce
