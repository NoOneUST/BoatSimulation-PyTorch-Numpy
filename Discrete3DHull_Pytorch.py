import itertools
from matplotlib import pyplot as plt
import numpy as np
import torch
from mpl_toolkits.mplot3d import Axes3D
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "3"

device = torch.device("cuda")


def showMat(Mat):
    return Mat.cpu().numpy()


for loopIter in range(1):
    print('This is Loop: ', loopIter)
    D = 0.3 * 0.65  # meter
    L, B = 0.3, 0.3

    xRange = 0.20 * 2 * 0.65  # meter
    yRange = 0.20 * 2 * 0.65  # meter
    zUpperBound = np.maximum(D * 1.2, 0.7)  # meter
    zLowerBound = D * 0.2

    precision = 0.0011  # meter
    gravity = 9.8  # m/s^2
    waterDensity = 997  # kg/m^3
    hullDensity = 1.1  # kg/m^3
    boatThick = 1e-2  # meter

    hullMeshDensity = hullDensity * np.power(precision, 3)
    waterMeshDensity = waterDensity * np.power(precision, 3)
    boatThickMesh = int(boatThick / precision)

    class ballast():
        weight = 0.12  # kg
        axisPosition = torch.tensor([0, 0, 0], device=device)  # meter
        meshPosition = ((axisPosition + zLowerBound) / precision).int()
        axisShape = torch.tensor(
            [0.03 * 2, 0.025 * 3, 0.025 * 1.5], device=device)  # meter
        meshShape = (axisShape / precision).int()
        meshVolume = meshShape[0] * meshShape[1] * meshShape[2]
        meshDensity = weight / meshVolume.float()

    ballast1 = ballast()

    MidIndexX = int(xRange / precision + 1)
    MidIndexY = int(yRange / precision + 1)
    MidIndexZ = int(zLowerBound / precision + 1)

    maxZIndex = int((D + zLowerBound + boatThick) / precision)

    xLen = int(2 * xRange / precision + 1)
    yLen = int(2 * yRange / precision + 1)
    zLen = int((zLowerBound + zUpperBound) / precision + 1)

    xArray = torch.linspace(-xRange, xRange, xLen, device=device)
    yArray = torch.linspace(-yRange, yRange, yLen, device=device)

    xGrid, yGrid = torch.meshgrid([xArray, yArray])
    zGrid = D * (torch.pow(2 * xGrid / L, 4) + torch.pow(2 * yGrid / B, 2))

    hullMesh = torch.zeros((xLen, yLen, zLen), device=device)

    minZIndexMat = (
        (torch.t(zGrid) +
         zLowerBound) /
        precision).int()

    zUpperBoundHull = int((D + zLowerBound + boatThick) / precision)

    for index in itertools.product(range(xLen), range(yLen)):
        hullMesh[index[0], index[1], minZIndexMat[index[0], index[1]]:min(maxZIndex, minZIndexMat[index[0], index[1]] + boatThickMesh)] = 1
        if minZIndexMat[index[0], index[1]] > 0:
            hullMesh[index[0], index[1], int((D + zLowerBound) / precision): maxZIndex] = 1


    def caculWeight(weightMat):
        return weightMat.sum()

    def calculCOM(weightMat):
        M = caculWeight(weightMat)
        xCOMPre = torch.dot(torch.linspace(
            1, xLen, xLen, device=device),
            torch.sum(
                weightMat, dim=(
                    2, 1)))
        yCOMPre = torch.dot(torch.linspace(
            1, yLen, yLen, device=device),
            torch.sum(
                weightMat, dim=(
                    0, 2)))
        zCOMPre = torch.dot(torch.linspace(
            1, zLen, zLen, device=device),
            torch.sum(
                weightMat, dim=(
                    1, 0)))
        xCOM, yCOM, zCOM = map(lambda x: x / M, [xCOMPre, yCOMPre, zCOMPre])
        return torch.tensor([xCOM, yCOM, zCOM], device=device).int()

    def addBallast(weightMat, ballast1):
        xCOM, yCOM, zCOM = calculCOM(weightMat)
        shape = ballast1.meshShape
        weightMat[int(xCOM -
                      shape[0] /
                      2):int(xCOM +
                             shape[0] /
                             2 +
                             1), int(yCOM -
                                     shape[1] /
                                     2):int(yCOM +
                                            shape[1] /
                                            2 +
                                            1), int(MidIndexZ +
                                                    ballast1.meshShape[2]):int(MidIndexZ +
                                                                               2 *
                                                                               ballast1.meshShape[2]) +
                  1] += ballast1.meshDensity
        return weightMat

    hullMesh = hullMeshDensity * hullMesh

    hullWeight = caculWeight(hullMesh)

    xCOM, yCOM, zCOM = calculCOM(hullMesh)
    print('COM: ', xCOM.item(), yCOM.item(), zCOM.item())

    print('Weight Before Ballast:', caculWeight(hullMesh).item())

    PureHullMeshNoDeck = torch.zeros(hullMesh.size())
    PureHullMeshNoDeck[:, :, :int((D + zLowerBound) / precision)] = hullMesh[:, :, :int((D + zLowerBound) / precision)]

    PureDeckMesh = torch.zeros(hullMesh.size())
    PureDeckMesh[:, :, int((D + zLowerBound) / precision):] = hullMesh[:, :, int((D + zLowerBound) / precision):]

    hullMesh = addBallast(hullMesh, ballast1)

    PureHullMeshWithDeckAndBallast = hullMesh + 0

    print('Weight After Ballast: ', caculWeight(hullMesh).item())

    def addMast(weightMat):
        DiaMeter = 9.5e-3  # meter
        Radius = DiaMeter / 2  # meter
        RadiusMesh = int(Radius / precision)
        Length = 0.5  # meter
        Weight = 96.7e-3  # kg
        Volume = np.pi * Radius * Radius * Length
        Density = Weight / Volume
        DensityMesh = Density * np.power(precision, 3)

        xCOM, yCOM, zCOM = calculCOM(weightMat)

        iArray = torch.linspace(- RadiusMesh, RadiusMesh, 2 *
                                RadiusMesh + 1, device=device)
        jRadius = torch.sqrt(
            np.power(RadiusMesh, 2) -
            (iArray ** 2)).long()
        iArray = iArray.long()
        for i in range(iArray.shape[0]):
            weightMat[iArray[i] +
                      xCOM, yCOM -
                      jRadius[i]:yCOM +
                      jRadius[i], int(zLowerBound /
                                      precision):int((zLowerBound +
                                                      Length) /
                                                     precision)] = DensityMesh
        return weightMat

    hullMesh = addMast(hullMesh)

    PureMastMesh = hullMesh - PureHullMeshWithDeckAndBallast

    print('Weight After Mast: ', caculWeight(hullMesh).item())

    xCOM, yCOM, zCOM = calculCOM(hullMesh)
    print('COM: ', xCOM.item(), yCOM.item(), zCOM.item(), '\n')

    hullWeight = caculWeight(hullMesh)

    DisplacementVolumeReal = hullWeight / waterDensity
    DisplacementVolumeMesh = DisplacementVolumeReal / np.power(precision, 3)

    def fliTrans(mat):
        return torch.flip(torch.t(mat), [0])

    def calculDisplacementVolumeMesh(
            hullMesh,
            waterAngle,
            waterOffset,
            isFinal,
            needInverse):
        noWaterMesh = torch.ones(
            hullMesh.shape[0],
            hullMesh.shape[1],
            hullMesh.shape[2],
            device=device)
        waterLineMesh = torch.zeros(
            hullMesh.shape[0],
            hullMesh.shape[1],
            hullMesh.shape[2],
            device=device)

        rightXYArray = ((np.tan(waterAngle) *
                         (torch.arange(yLen, device=device).float() *
                          precision -
                          xRange) +
                         yRange) /
                        precision +
                        waterOffset).int()
        rightXYArray = torch.max(
            rightXYArray, torch.zeros(
                yLen, device=device).int())

        if needInverse:
            for i in range(yLen):
                noWaterMesh[:, i, rightXYArray[i]:] = 0
        else:
            for i in range(yLen):
                noWaterMesh[:, i, :rightXYArray[i]] = 0

        hullArea = (hullMesh > 0).int()
        waterArea = (~ noWaterMesh.byte() > 0).int()
        unSubmergedMesh = (noWaterMesh.int() & hullArea).float() * hullMesh
        DisplacedWaterWeight = torch.sum(
            waterArea & hullArea).float() * waterMeshDensity
        DisplacementVolume = torch.sum(
            waterArea & hullArea).float() * np.power(precision, 3)
        SubmergedVolume = torch.sum(
            waterArea & hullArea).float() * np.power(precision, 3)
        SubmergedMesh = (waterArea & hullArea).float() * hullMesh

        if isFinal:
            print('\nSubmergedVolume: ', SubmergedVolume.item())
            print(
                'Weight of Displaced Water: ',
                DisplacementVolume.item() *
                waterDensity)
            print('Weight of Boat: ', hullWeight.item())

        return [
            unSubmergedMesh,
            SubmergedMesh,
            waterLineMesh,
            DisplacedWaterWeight,
            noWaterMesh]

    def calculBestWaterOffsetMesh(hullMesh, waterAngle, hullWeight):
        lastLoss = - 1
        maxIteration = 5000
        waterOffsetLowerBound = - 2 * zLen
        waterOffsetUpperBound = 3 * zLen
        waterOffset = int((waterOffsetLowerBound + waterOffsetUpperBound) / 2)

        if (waterAngle >= 0 and waterAngle < np.pi /
                2) or (waterAngle >= 1.5 * np.pi and waterAngle < 2 * np.pi):
            needInverse = 0
        else:
            needInverse = 1

        for i in range(maxIteration):
            print('Processed: ' +
                  str((i +
                       1) /
                      maxIteration *
                      100) +
                  '%' +
                  '     Loss: ' +
                  str(lastLoss) +
                  '     waterOffsetLowerBound: ' +
                  str(waterOffsetLowerBound) +
                  '     waterOffsetUpperBound: ' +
                  str(waterOffsetUpperBound))
            unSubmergedMesh, SubmergedMesh, waterLineMesh, DisplacedWaterWeight, noWaterMesh = calculDisplacementVolumeMesh(
                hullMesh, waterAngle, waterOffset, False, needInverse)
            thisLoss = (DisplacedWaterWeight - hullWeight) / hullWeight

            if i >= maxIteration - \
                    1 or abs(thisLoss) < 0.01 or waterOffsetLowerBound >= waterOffsetUpperBound - 1:
                unSubmergedMesh, SubmergedMesh, waterLineMesh, DisplacedWaterWeight, noWaterMesh = calculDisplacementVolumeMesh(
                    hullMesh, waterAngle, waterOffset, True, needInverse)
                thisLoss = (DisplacedWaterWeight - hullWeight) / hullWeight
                print(
                    '\nDone \nLoss= ',
                    thisLoss.item() * 100,
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
                    noWaterMesh]
            else:
                if thisLoss < 0 and not(needInverse):
                    waterOffsetLowerBound = waterOffset
                    waterOffset = int(
                        (waterOffset + waterOffsetUpperBound) / 2)
                elif thisLoss < 0 and needInverse:
                    waterOffsetUpperBound = waterOffset
                    waterOffset = int(
                        (waterOffset + waterOffsetLowerBound) / 2)
                elif thisLoss > 0 and not(needInverse):
                    waterOffsetUpperBound = waterOffset
                    waterOffset = int(
                        (waterOffset + waterOffsetLowerBound) / 2)
                else:
                    waterOffsetLowerBound = waterOffset
                    waterOffset = int(
                        (waterOffset + waterOffsetUpperBound) / 2)
                lastLoss = thisLoss.item()

    degreeAngle = 120  # degree

    waterAngle = np.deg2rad(degreeAngle)

    waterOffset, unSubmergedMesh, SubmergedMesh, waterLineMesh, Loss, noWaterMesh = calculBestWaterOffsetMesh(
        hullMesh, waterAngle, hullWeight)

    HullSubmergedMesh = (noWaterMesh.float() < 0.5).float() * SubmergedMesh

    [xCOB, yCOB, zCOB] = calculCOM(HullSubmergedMesh)

    FBuoyancy = FGravity = (hullWeight * gravity).item()

    buoyancyTorqueVector = np.cross([(yCOB - yCOM).item() * precision,
                                     (zCOB - zCOM).item() * precision,
                                     0],
                                    [-FBuoyancy * np.sin(waterAngle),
                                        FBuoyancy * np.cos(waterAngle), 0])

    buoyancyTorque = buoyancyTorqueVector[2]

    print(
        '\nBuoyancy Torque under ',
        degreeAngle,
        ' degrees: ',
        buoyancyTorque,
        ' N * M\n')

    print('COM: ', list(map(lambda x: x.item(), [xCOM, yCOM, zLen - zCOM])))
    print('COB: ', list(map(lambda x: x.item(), [xCOB, yCOB, zLen - zCOB])))

    if(buoyancyTorque > 0):
        print('\nIt can recover\n')

    plt.figure(1)
    plt.matshow(fliTrans(SubmergedMesh[MidIndexX, :, :] > 0).cpu().numpy())
    plt.title('Submerged_yz ' + str(degreeAngle))
    plt.show()

    plt.figure(2)
    plt.matshow(fliTrans(hullMesh[MidIndexX, :, :] > 0).cpu().numpy())
    plt.title('yz ' + str(degreeAngle))
    plt.show()

    plt.figure(3)
    plt.matshow(fliTrans(hullMesh[:, MidIndexY, :] > 0).cpu().numpy())
    plt.title('xz')
    plt.show()

    zWeightArray = torch.sum(hullMesh > 0, dim=(1, 0))
    xNoZeroIndexArray = (torch.sum(hullMesh, dim=(1, 2)) > 0).nonzero()
    yNoZeroIndexArray = (torch.sum(hullMesh, dim=(2, 0)) > 0).nonzero()

    boatHeightIndex = (
        (zWeightArray[:zLen - 2] - 1.1 * zWeightArray[1:zLen - 1]) > 0).nonzero()[-2][0]

    plt.figure(4)
    plt.matshow((fliTrans(hullMesh[:, :, boatHeightIndex]) > 0).cpu().numpy())
    plt.title('xy')
    plt.show()

    print(
        'The height of Boat is: ',
        (boatHeightIndex -
         (zWeightArray > 0).nonzero()[0][0]).item() *
        precision)
    print('The Length of Boat is: ',
          (xNoZeroIndexArray[-1][0] - xNoZeroIndexArray[0][0]).item() * precision)
    print('The Width of Boat is: ',
          (yNoZeroIndexArray[-1][0] - yNoZeroIndexArray[0][0]).item() * precision)

    synthesisMap = noWaterMesh + (hullMesh[MidIndexX, :, :] > 0).float()
    synthesisMap[xCOM - 2:xCOM + 2:, yCOM - 2:yCOM + 2, zCOM - 2:zCOM + 2] += 1
    synthesisMap[xCOB - 2:xCOB + 2, yCOB - 2:yCOB + 2, zCOB - 2:zCOB + 2] += 1

    plt.figure(5)
    plt.matshow(fliTrans(synthesisMap[xCOM, :, :]).cpu().numpy())
    plt.title('synthesisMap')
    plt.show()

    poolingOpt = torch.nn.AdaptiveAvgPool3d(200)
    PureHullMeshNoDeck, PureDeckMesh, PureMastMesh = list(map(lambda x:poolingOpt(x.unsqueeze(
        0)).squeeze().nonzero().cpu().numpy(), [PureHullMeshNoDeck, PureDeckMesh, PureMastMesh]))

    plt.figure(6)
    ax = plt.subplot(111, projection='3d')
    ax.scatter(PureHullMeshNoDeck[:, 0], PureHullMeshNoDeck[:, 1], PureHullMeshNoDeck[:, 2], c='b')
    ax.scatter(PureDeckMesh[:, 0], PureDeckMesh[:, 1], PureDeckMesh[:, 2], c='w')
    ax.scatter(PureMastMesh[:, 0], PureMastMesh[:, 1], PureMastMesh[:, 2], c='y')
    ax.set_zlabel('Z')
    ax.set_ylabel('Y')
    ax.set_xlabel('X')
    plt.show()

    pass
