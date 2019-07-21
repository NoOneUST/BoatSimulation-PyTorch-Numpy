from matplotlib import pyplot as plt
# import multiprocessing
import numpy as np
# from pathos.multiprocessing import ProcessingPool as Pool

from scipy.optimize import root, fsolve
from mpl_toolkits.mplot3d import Axes3D

# cores = multiprocessing.cpu_count()

# x, y = 6, 2
D = 0.2  # meter
L, B = 0.5, 0.5

xRange = 0.305 * 2  # meter
yRange = 0.23 * 2  # meter
zUpperBound = D * 1.2 + 0.5  # meter
zLowerBound = D * 0.2

precision = 0.004  # meter
gravity = 9.8  # m/s^2
waterDensity = 997  # kg/m^3
hullDensity = 110  # kg/m^3
boatThick = 1e-2  # meter

hullMeshDensity = hullDensity * np.power(precision, 3)
waterMeshDensity = waterDensity * np.power(precision, 3)


def axis2MeshPositionZ(axis):
    return list(map(int, (np.array(axis) + zLowerBound) / precision))


class ballast():
    weight = 0.6  # kg
    axisPosition = np.array([0, 0, 0])  # meter
    meshPosition = axis2MeshPositionZ(axisPosition)
    axisShape = np.array([0.06, 0.025, 0.025])  # meter
    meshShape = list(map(int, axisShape / precision))
    meshVolume = meshShape[0] * meshShape[1] * meshShape[2]
    meshDensity = weight / meshVolume


ballast1 = ballast()

MidIndexX = int(xRange / precision + 1)
MidIndexY = int(yRange / precision + 1)
MidIndexZ = int(zLowerBound / precision + 1)

maxZIndex = int((D + zLowerBound) / precision)

xLen = int(2 * xRange / precision + 1)
yLen = int(2 * yRange / precision + 1)
zLen = int((zLowerBound + zUpperBound) / precision + 1)

xArray = np.linspace(-xRange, xRange, xLen)
yArray = np.linspace(-yRange, yRange, yLen)

xGrid, yGrid = np.meshgrid(xArray, yArray)
zGrid = D * (np.power(2 * xGrid / L, 4) + np.power(2 * yGrid / B, 2))

# def f(x, y):
#     return x*x + y*y
#
# pool = Pool(cores)
#
# pool.map(f, range(5), range(5))

hullMesh = np.zeros((xLen, yLen, zLen))

for i in range(xLen):
    for j in range(yLen):
        minZIndex = int((zGrid[j, i] + zLowerBound) / precision)
        hullMesh[i, j, minZIndex:maxZIndex] = 1
        # hullMesh[i, j, minZIndex:minZIndex + int(boatThick / precision)] = 1
        # hullMesh[i, j, maxZIndex - int(boatThick / precision):maxZIndex] = 1


def caculWeight(weightMat):
    M = np.sum(weightMat)
    return M


def calculCOM(weightMat):
    M = caculWeight(weightMat)
    xCOMPre = np.linspace(1, xLen, xLen).dot(np.sum(weightMat, axis=(2, 1)))
    yCOMPre = np.linspace(1, yLen, yLen).dot(np.sum(weightMat, axis=(0, 2)))
    zCOMPre = np.linspace(1, zLen, zLen).dot(np.sum(weightMat, axis=(1, 0)))
    xCOM = xCOMPre / M
    yCOM = yCOMPre / M
    zCOM = zCOMPre / M
    return [xCOM, yCOM, zCOM]


def addBallast(weightMat, ballast1):
    xCOM, yCOM, zCOM = map(int, calculCOM(weightMat))
    shape = ballast1.meshShape
    weightMat[int(xCOM - shape[0] / 2):int(xCOM + shape[0] / 2 + 1),
              int(yCOM - shape[1] / 2):int(yCOM + shape[1] / 2 + 1),
              int(MidIndexZ):int(MidIndexZ + shape[2]) + 1] += ballast1.meshDensity
    return weightMat


hullMesh = hullMeshDensity * hullMesh

hullWeight = caculWeight(hullMesh)

xCOM, yCOM, zCOM = list(map(int, calculCOM(hullMesh)))
print([xCOM, yCOM, zCOM])

hullMesh = addBallast(hullMesh, ballast1)


def addMast(weightMat):
    DiaMeter = 9.5e-3  # meter
    Radius = DiaMeter / 2  # meter
    RadiusMesh = Radius / precision
    Length = 0.5  # meter
    Weight = 96.7e-3  # kg
    Volume = np.pi * Radius * Radius * Length
    Density = Weight / Volume
    DensityMesh = Density * np.power(precision, 3)

    xCOM, yCOM, zCOM = map(int, calculCOM(weightMat))

    for i in np.linspace(
            xCOM - RadiusMesh,
            xCOM + RadiusMesh,
            int(2 * RadiusMesh + 1)):
        for j in np.linspace(
                yCOM - RadiusMesh,
                yCOM + RadiusMesh,
                int(2 * RadiusMesh + 1)):
            if (i - xCOM) * (i - xCOM) + (j - yCOM) * \
                    (j - yCOM) < RadiusMesh * RadiusMesh:
                weightMat[int(i), int(j), int(zLowerBound / precision):int((zLowerBound + Length) / precision)] = DensityMesh

    return weightMat


hullMesh = addMast(hullMesh)

xCOM, yCOM, zCOM = list(map(int, calculCOM(hullMesh)))
print([xCOM, yCOM, zCOM], '\n')

hullWeight = caculWeight(hullMesh)

DisplacementVolumeReal = hullWeight / waterDensity
DisplacementVolumeMesh = DisplacementVolumeReal / pow(precision, 3)


def fliTrans(mat):
    return np.flipud(mat.transpose())


# plt.figure(4)
# plt.matshow(fliTrans(hullMesh[MidIndexX, :, :]))
# plt.show()

def calculDisplacementVolumeMesh(hullMesh, waterAngle, waterOffset, isFinal):
    noWaterMesh = np.ones(
        [hullMesh.shape[0], hullMesh.shape[1], hullMesh.shape[2]])
    waterLineMesh = np.zeros(
        [hullMesh.shape[0], hullMesh.shape[1], hullMesh.shape[2]])

    if (waterAngle >= 0 and waterAngle < np.pi /
            2) or (waterAngle >= 1.5 * np.pi and waterAngle < 2 * np.pi):
        needInverse = 0
    else:
        needInverse = 1
    for i in range(yLen):
        rightXY = (np.tan(waterAngle) * (i * precision - xRange) +
                   yRange) / precision + waterOffset
        rightX1Y = (np.tan(waterAngle) * ((i - 1) * precision -
                                          xRange) + yRange) / precision + waterOffset

        for j in range(zLen):
            if np.logical_xor(j < rightXY, needInverse):
                noWaterMesh[:, i, j] = 0
            if np.logical_xor(
                    j - 1 < rightXY,
                    j < rightXY) or np.logical_xor(
                    j < rightX1Y,
                    j < rightXY):
                waterLineMesh[:, i, j] = 1

    unSubmergedMesh = np.logical_and(noWaterMesh, hullMesh) * hullMesh
    DisplacedWaterWeight = np.sum(
        np.logical_and(
            np.logical_not(noWaterMesh),
            hullMesh)) * waterMeshDensity
    DisplacementVolume = np.sum(np.logical_and(
        np.logical_not(noWaterMesh), hullMesh)) * np.power(precision, 3)
    SubmergedVolume = np.sum(np.logical_and(
        np.logical_not(noWaterMesh),
        hullMesh)) * np.power(precision, 3)
    SubmergedMesh = np.logical_and(
        np.logical_not(noWaterMesh),
        hullMesh) * hullMesh

    if isFinal:
        print('\nSubmergedVolume: ', SubmergedVolume)
        print('Weight of Displaced Water: ', DisplacementVolume * waterDensity)
        print('Weight of Boat: ', hullWeight)

    return [
        unSubmergedMesh,
        SubmergedMesh,
        waterLineMesh,
        DisplacedWaterWeight]

# def calLoss(waterOffset):
#     unSubmergedMesh, SubmergedMesh, waterLineMesh, SubmergedWeight = calculDisplacementVolumeMesh(
#         hullMesh, waterAngle, waterOffset)
#     thisLoss = SubmergedWeight - hullWeight
#     return thisLoss
#
# degreeAngle = 30
#
# waterAngle = np.deg2rad(degreeAngle)
# print('root: ', root(calLoss, [58]))
# #print('fsolve: ', fsolve(calLoss, [55]))


def calculBestWaterOffsetMesh(hullMesh, waterAngle, hullWeight):
    lastLoss = np.sum(hullMesh)
    maxIteration = 5000
    waterOffsetLowerBound = - 2 * zLen
    waterOffsetUpperBound = 3 * zLen
    waterOffset = waterOffsetLowerBound
    lastWaterOffset = waterOffset - 1
    intersect = (waterOffsetUpperBound - waterOffsetLowerBound) / \
        (maxIteration - 1)
    for i in range(maxIteration):
        print('Processed: ' + str((i + 1) / maxIteration * 100) +
              '%' + '     Loss: ' + str(lastLoss))
        unSubmergedMesh, SubmergedMesh, waterLineMesh, DisplacedWaterWeight = calculDisplacementVolumeMesh(
            hullMesh, waterAngle, waterOffset, False)
        thisLoss = DisplacedWaterWeight - hullWeight

        if (abs(thisLoss) - abs(lastLoss) > 0.001 and waterOffset - lastWaterOffset > 0.01 and (lastLoss <
                                                                                                0 or thisLoss < 0)) or i >= maxIteration - 1 or abs(thisLoss / hullWeight * 100) < 1:
            unSubmergedMesh, SubmergedMesh, waterLineMesh, DisplacedWaterWeight = calculDisplacementVolumeMesh(
                hullMesh, waterAngle, lastWaterOffset, True)
            thisLoss = DisplacedWaterWeight - hullWeight
            print(
                '\nDone \nLoss= ',
                thisLoss /
                hullWeight *
                100,
                '%',
                '\nOffset = ',
                waterOffset)
            pass
            return [
                waterOffset,
                unSubmergedMesh,
                SubmergedMesh,
                waterLineMesh,
                thisLoss,
                thisLoss /
                np.sum(hullMesh)]
        else:
            if abs(lastLoss - thisLoss) / \
                    abs(thisLoss) < 0.05 and abs(thisLoss / hullWeight * 100) > 30:
                waterOffset += intersect * 10
            elif abs(thisLoss / hullWeight * 100) < 10:
                waterOffset += intersect * 0.1
            else:
                waterOffset += intersect
            lastLoss = thisLoss
            lastWaterOffset = waterOffset

    #
    # for waterOffset in np.linspace(waterOffsetLowerBound, waterOffsetUpperBound, maxIteration):
    #     print('Processed: ' + str((waterOffset + zLen / 2) / zLen * 100) +
    #           '%' + '     Loss: ' + str(lastLoss))
    #     unSubmergedMesh, SubmergedMesh, waterLineMesh, SubmergedWeight = calculDisplacementVolumeMesh(
    #         hullMesh, waterAngle, waterOffset)
    #     thisLoss = SubmergedWeight - hullWeight
    #     # plt.figure(1)
    #     # plt.matshow(fliTrans(unSubmergedMesh[MidIndexX, :, :]))
    #     # plt.title('unSubmergedMesh')
    #     # plt.show()
    #     # plt.figure(2)
    #     # plt.matshow(fliTrans(SubmergedMesh[MidIndexX, :, :]))
    #     # plt.title('noWaterMesh')
    #     # plt.show()
    #     if abs(thisLoss) > abs(lastLoss + abs(0.1 * lastLoss)) or (waterOffset + zLen / 2) / zLen > (maxIteration - 1) / maxIteration:
    #         unSubmergedMesh, SubmergedMesh, waterLineMesh, SubmergedWeight = calculDisplacementVolumeMesh(
    #             hullMesh, waterAngle, waterOffset - intersect)
    #         thisLoss = SubmergedWeight - hullWeight
    #         print('done, Loss= ', thisLoss / hullWeight *100, '%', '   offset = ', waterOffset)
    #         pass
    #         return [
    #             waterOffset,
    #             unSubmergedMesh,
    #             SubmergedMesh,
    #             waterLineMesh,
    #             thisLoss,
    #             thisLoss /
    #             np.sum(hullMesh)]
    #     else:
    #         lastLoss = thisLoss

    # plt.figure(1)
    # plt.matshow(fliTrans(unSubmergedMesh[MidIndexX, :, :]))
    # plt.title('unSubmergedMesh')
    # plt.show()
    # plt.figure(2)
    # plt.matshow(fliTrans(noWaterMesh[MidIndexX, :, :]))
    # plt.title('noWaterMesh')
    # plt.show()


degreeAngle = 30  # degree

waterAngle = np.deg2rad(degreeAngle)

waterOffset, unSubmergedMesh, SubmergedMesh, waterLineMesh, Loss, lossPercent = calculBestWaterOffsetMesh(
    hullMesh, waterAngle, hullWeight)


# def calculCOB(SubmergedMesh):
#     return calculCOM(SubmergedMesh)


[xCOB, yCOB, zCOB] = map(int, calculCOM(SubmergedMesh))

FBuoyancy = FGravity = hullWeight * gravity

buoyancyTorqueVector = np.cross([(yCOB - yCOM) * precision,
                                 (zCOB - zCOM) * precision,
                                 0],
                                [-FBuoyancy * np.sin(waterAngle),
                                 FBuoyancy * np.cos(waterAngle), 0])

buoyancyTorque = -buoyancyTorqueVector[2]

print(
    '\nBuoyancy Torque under ',
    degreeAngle,
    ' degrees: ',
    buoyancyTorque,
    ' N * M\n')

print('COM: ', [xCOM, yCOM, zLen - zCOM])
print('COB: ', [xCOB, yCOB, zLen - zCOB])

if(buoyancyTorque > 0):
    print('\nIt can recover\n')

# hullMesh = SubmergedMesh

# plt.figure(1)
# plt.matshow(np.array([[int(0), int(0)], [int(1), int(1)]]))
# plt.title('test')
# plt.show()

plt.figure(1)
# plt.matshow(fliTrans(hullMesh[MidIndexX, :, :]))
plt.matshow(fliTrans(SubmergedMesh[MidIndexX, :, :] > 0))
plt.title('Submerged_yz ' + str(degreeAngle))
plt.show()

plt.figure(2)
# plt.matshow(fliTrans(hullMesh[MidIndexX, :, :]))
plt.matshow(fliTrans(hullMesh[MidIndexX, :, :] > 0))
plt.title('yz ' + str(degreeAngle))
plt.show()

plt.figure(3)
plt.matshow(fliTrans(hullMesh[:, MidIndexY, :] > 0))
plt.title('xz')
plt.show()

plt.figure(4)
plt.matshow(fliTrans(hullMesh[:, :, MidIndexZ] > 0))
plt.title('xy')
plt.show()

# fig1 = plt.figure(1)
# ax = Axes3D(fig1)
#
# ax.plot_surface(
#     xGrid,
#     yGrid,
#     zGrid,
#     rstride=1,
#     cstride=1,
#     cmap='rainbow')
# # ax.plot_surface(xGrid, yGrid, DPlane, rstride=1, cstride=1, cmap='rainbow')
#
# plt.show()

pass
