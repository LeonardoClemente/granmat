import math
import math as mt
import numpy as np
import numpy
import sphere
import granmat
import granmat as gm
import os
import os.path as pa
import scipy.io as sio
import pickle
import glob
import shutil
import time
import numpy.matlib
''' Class to aid in granular materials simulations along with the Sphere library
 By Leonardo Clemente and Eli cuellar. Sphere is made by Anders Damsgaard'''
'''debugging keywords : CAMBIAR (need changes)
    TESTEAR (needs more testing)
'''




class granmatclass:
    def __init__(self, N=32500, r=.0052/2, boxDimensions=np.array([.75, .30, 4])):
        # metadata
        self.simId = 'ricardop08p001' # ricardo+cylinder separation + wall separation
        self.origOutput= '/home/leon/Downloads/sphere-master/output' # Original output
        self.extOutput = "/media/leon/leo1/sphere/output/" # Your output
        self.extInput = "/media/leon/leo1/sphere/input/"
        self.extSave = 1

        # experiment data
        self.N = N  # Number of particles
        self.r = r  # particle radius
        self.particleDiameter = r*2
        self.dx = 2.1*self.r

        self.cylinderSeparation = .08
        self.lWalls = .001
        # Rescaling box dimensions
        self.boxDimensions = boxDimensions
        self.boxDimensions[0] = int(self.boxDimensions[0]/self.dx)*self.dx
        self.boxDimensions[1] = int(self.boxDimensions[1]/self.dx)*self.dx
        self.boxDimensions[2] = self.boxDimensions[2]*self.dx
        self.V = (self.r**3)*mt.pi*4/3  # volume

        # experiment physical parameters *(manually set using setPhysicalParams function)

        self.g = -9.8  # gravity
        self.muS = 0.4  # static friction
        self.muD = 0.4  # dynamic friction
        self.kN = 4084  # normal spring constant
        self.kT = 4084  # tangential spring constant
        self.eps = .5  # restitution coefficient
        self.rho = 14
        self.m = self.V*self.rho  # particle mass
        self.mEff = (self.m**2)/(2*self.m)
        self.s = ((.5*np.log(self.eps)*mt.pi)**2+1)/self.mEff  # parameter to calculate restitution coefficient
        self.gam = (self.kN/self.s)**1/2  # viscosity 0.00385218831922*16
        # booleans
        self.lowerWallsBool = 0
        self.conversionFactorsBool = 0

    def setConversionFactors(self, t):
        if t == 1:
            print('Conversion Factors are set.')
            self.fcK = self.particleDiameter/(self.m*abs(self.g))
            self.fcGamma = mt.sqrt(self.particleDiameter*abs(self.g))
            self.fcD = 1/self.particleDiameter
            self.fcF = 1/(self.m*abs(self.g))
            self.fcT = mt.sqrt(abs(self.g)/self.particleDiameter)
            self.fcG = 1/abs(self.g)
            self.fcM = 1/self.m
            self.conversionFactorsBool = 1
        if t == 0:
            print('Conversion factors deleted. Warning! : Some simulations may not run without adimensionalization of parameters.')
            self.fcK = 1
            self.fcGamma = 1
            self.fcD = 1
            self.fcF = 1
            self.fcT = 1
            self.fcG = 1
            self.fcM = 1
            self.conversionFactorsBool = 0

    def setPhysicalParams(self, g=-9.8, muS=0.4, muD=0.4, kN=4084, kT=4084, eps=.5, rho=14):
        self.g = g  # gravity
        self.muS = muS  # static friction
        self.muD = muD  # dynamic friction
        self.kN = kN  # normal spring constant
        self.kT = kT  # tangential spring constant
        self.eps = eps  # restitution coefficient
        self.rho = rho

    def setLowerWalls(self, nWalls=9, lWalls=None, relativeRes=5, wallSepX=.052):

        if lWalls == None:
            lWalls = self.lWalls
        self.lowerWallsBool = (self.lowerWallsBool+1) % 2
        if self.lowerWallsBool == 1:
            print 'Lower walls option activaded.'
            self.nWalls = nWalls
            self.lWalls = lWalls
            self.relativeRes = relativeRes
            self.wallSepX = wallSepX
        else:
            print 'Lower walls options desactivaded.'

    '''
    def simdata2txt(self,output='', simId=' ', ext=0):  #TESTEAR
        if self.extSave == 1:
            output = self.extOutput
        else:
            output = self.origOutput
        if simId == ' ':
            simId = 'ricardop06p01Sim'

        dataDir = output+'/'+simId
        if pa.isdir(dataDir) == 0:
            os.mkdir(dataDir)

        nSteps = 20000  # CAMBIAR

        # cargar simulacion y entrar al ciclo
        SB = sphere.sim(np=2, sid=simId)
        if ext == 1:
            SB.extsave = 1
        SB.readstep(1)
        #sizes
        n = np.size(SB.x, 0)
        m = np.size(SB.x, 1)

        #Limits for substractBox function TEMPORAL
        limits=self.boxDimensions*self.fcD
        box = [ limits[0]*.6, limits[0]*.7, limits[1], 0, limits[2], 0]
        # vList = []

        vFile = open(self.extOutput+'velocities' + simId + '.csv', 'w')
        xFile = open(self.extOutput+'positions' + simId + '.csv' ,'w')
        fFile = open(self.extOutput+'forces'+simId + '.csv','w')

        for j in range(1, nSteps):
            print j
            # Getting particles
            SB.readstep(j)
            # Writing down on txt files
            for x in range (0,n):
                fFile.write(repr(SB.force[x,0]) + ' , ' + repr(SB.force[x,1]) + ' , ' +repr(SB.force[x,2]) + '\n')
                vFile.write(repr(SB.vel[x,0]) + ' , ' + repr(SB.vel[x,1]) + ' , ' +repr(SB.vel[x,2]) + '\n')
                xFile.write(repr(SB.x[x,0]) + ' , ' + repr(SB.x[x,1]) + ' , ' +repr(SB.x[x,2]) + '\n')
            # Checking for particle velocities in delimited area
            particlesIn = gm.substractBox(vs, box)
            if numpy.size(particlesIn) > 0:
                indices = [index for index, value in enumerate(particlesIn) if value == 1];
                velocities = vs[indices, :]
                vList.append(velocities)
            else:
                vList.append([0, 0])

        # sio.savemat(dataDir+'/fxv.mat', dict(x=x, v=v, f=f))
        fFile.close()
        vFile.close()
        xFile.close()
        vFile = open(dataDir+'/vList.obj', 'w')
        pickle.dump(vList, vFile)
        vFile.close()

    '''

    def simdata2mat(self,output='', simId='ricardop01p', ext=0):  #TESTEAR
        if self.extSave == 1:
            output = self.extOutput
        else:
            output = self.origOutput
        if simId == ' ':
            simId == self.simId

        dataDir = output+'/'+simId
        if pa.isdir(dataDir) == 0:
            os.mkdir(dataDir)

        nSteps = 10000  # CAMBIAR

        # cargar simulacion y entrar al ciclo
        SB = sphere.sim(np=2, sid=simId)
        if ext == 1:
            SB.extsave == 1
        SB.readstep(1)
        #sizes
        n = np.size(SB.x, 0)
        m = np.size(SB.x, 1)

        #Limits for substractBox function TEMPORAL
        limits=self.boxDimensions*self.fcD
        box = [ limits[0]*.6, limits[0]*.7, limits[1], 0, limits[2], 0]
        vList = []

        v = np.zeros([n, m, nSteps])
        x = np.zeros([n, m, nSteps])
        f = np.zeros([n, m, nSteps])
        for j in range(1, nSteps):
            print j
            # Getting particles
            SB.readstep(j)
            x[:, :, j] = SB.x[:, :]
            vs = SB.vel[:, :]
            v[:, :, j] = vs
            f[:, :, j] = SB.force[:, :]

            # Checking for particle velocities in delimited area
            particlesIn = gm.substractBox(vs, box)
            if numpy.size(particlesIn) > 0:
                indices = [index for index, value in enumerate(particlesIn) if value == 1];
                velocities = vs[indices, :]
                vList.append(velocities)
            else:
                vList.append([0, 0])

        sio.savemat(dataDir+'/fxv.mat', dict(x=x, v=v, f=f))

        vFile = open(dataDir+'/vList.obj', 'w')
        pickle.dump(vList, vFile)
        vFile.close()

    def data2txt(self,output='', simId='', ext=0):
        # This function sums the forces from each of the cylinders for every
        # in the simulation and saves them in a  '.m' format matrix.
        if self.extSave == 1:
            output = self.extOutput
        else:
            output = self.origOutput
        if simId == '':
            simId = self.simId

        dataDir = self.extOutput+simId

        if pa.isdir(dataDir) == 0:
            print 'Could not find {0} directory.'.format(simId)
            return


        file = open(dataDir+'/cylinderInfo.obj')
        file2 = open(dataDir+'/simInfo.obj')
        cylinderInfo=pickle.load(file)
        simulationInfo=pickle.load(file2)
        file.close()
        file2.close()
        cylinder1 = open(dataDir+'/cylinder1Info.csv', 'w')
        cylinder2 = open(dataDir+'/cylinder2Info.csv', 'w')
        vFile = open(dataDir+'/velocities.csv', 'w')
        xFile = open(dataDir+'/positions.csv','w')
        # fFile = open(self.extOutput+'forces'+simId + '.csv','w')

        n1=cylinderInfo[0]
        n2=cylinderInfo[1]
        x1=cylinderInfo[2]
        x2=cylinderInfo[3]
        #nSteps=simulationInfo[0]
        nSteps = 19999

        # Entering simulation

        SB = sphere.sim(np=2, sid=simId+'Sim')
        SB.extsave == 1
        SB.readstep(1)

        n = np.size(SB.x, 0)
        m = np.size(SB.x, 1)

        #Limits for substractBox function TEMPORAL
        limits=self.boxDimensions*self.fcD
        box = [ limits[0]*.6, limits[0]*.7, limits[1], 0, limits[2], 0]
        vList = []

        # Looping on every step
        r1=numpy.array([x1[0],x1[1],0])
        r2=numpy.array([x2[0],x2[1],0])
        f1=numpy.zeros([nSteps,3])
        f2=numpy.zeros_like(f1)
        print n1, n2
        try:
            for i in range(1, nSteps):
                SB.readstep(i)
                Fcilinder1 = np.array([0,0,0])
                Fcilinder2 = Fcilinder1
                for j in range(0, n1-1):
                    r1[2] = SB.x[j, 2]
                    rpart = SB.x[j, :]
                    u = r1 - rpart
                    u = u/np.linalg.norm(u)
                    Ftotal = SB.force[j, :]
                    Frad = u*np.dot(Ftotal, u)
                    Fcilinder1 = Fcilinder1 + Frad
                for j in range(n1,n2):
                    r2[2] = SB.x[j, 2]
                    rpart = SB.x[j, :]
                    u = r2 - rpart
                    u = u/np.linalg.norm(u)
                    Ftotal = SB.force[j, :]
                    Frad = u*np.dot(Ftotal, u)
                    Fcilinder2 = Fcilinder2 + Frad
                # Write forces to matrices
                # f1[i, :] = Fcyl1[:]
                # f2[i, :] = Fcyl2[:]
                cylinder1.write(repr( Fcilinder1[0] ) + ' , ' + repr( Fcilinder1[1] ) + ' , ' +repr( Fcilinder1[2]) + '\n')
                cylinder2.write(repr( Fcilinder2[0] ) + ' , ' + repr( Fcilinder2[1] ) + ' , ' +repr( Fcilinder2[2]) + '\n')

                if i % 40 == 0:
                    for x in range (0,n):
                        # fFile.write(repr(SB.force[x,0]) + ' , ' + repr(SB.force[x,1]) + ' , ' +repr(SB.force[x,2]) + '\n')
                        vFile.write(repr(SB.vel[x,0]) + ' , ' + repr(SB.vel[x,1]) + ' , ' +repr(SB.vel[x,2]) + '\n')
                        xFile.write(repr(SB.x[x,0]) + ' , ' + repr(SB.x[x,1]) + ' , ' +repr(SB.x[x,2]) + '\n')
        except error:
            print "An error ocurred. Exiting writing loop and closing files."
        # sio.savemat(dataDir+'/cylinderForces.mat', dict(f1=f1, f2=2))
        # fFile.close()
        vFile.close()
        xFile.close()
        cylinder1.close()
        cylinder2.close()
    def alldata2txtl(self,output='', simId='', ext=0):  # TESTEAR
        # This function sums the forces from each of the cylinders for every
        # in the simulation and saves them in a  '.m' format matrix.
        if self.extSave == 1:
            output = self.extOutput
        else:
            output = self.origOutput
        if simId == '':
            simId = self.simId

        dataDir = self.extOutput+simId

        if pa.isdir(dataDir) == 0:
            print 'Could not find {0} directory.'.format(simId)
            return


        file = open(dataDir+'/cylinderInfo.obj')
        file2 = open(dataDir+'/simInfo.obj')
        cylinderInfo=pickle.load(file)
        simulationInfo=pickle.load(file2)
        file.close()
        file2.close()
        cylinder1 = open(dataDir+'/cylinder1Info.csv', 'w')
        cylinder2 = open(dataDir+'/cylinder2Info.csv', 'w')
        vFile = open(dataDir+'/velocities.csv', 'w')
        xFile = open(dataDir+'/positions.csv','w')
        # fFile = open(self.extOutput+'forces'+simId + '.csv','w')

        n1 = cylinderInfo[0]
        n2 = cylinderInfo[1]
        x1 = cylinderInfo[2]
        x2 = cylinderInfo[3]
        nSteps = simulationInfo[0]
        #nPart = simulationInfo[1]
        nPart = 47062
        nSteps = 20000

        # Entering simulation

        SB = sphere.sim(np=nPart, sid=simId+'Sim')
        SB.extsave == 1
        #SB.readstep(1)
        m = np.size(SB.x, 1)

        #Limits for substractBox function TEMPORAL
        limits=self.boxDimensions*self.fcD
        box = [ limits[0]*.6, limits[0]*.7, limits[1], 0, limits[2], 0]
        vList = []

        # Looping on every step
        r1=numpy.array([x1[0],x1[1],0])
        r2=numpy.array([x2[0],x2[1],0])
        f1=numpy.zeros([nSteps,3])
        f2=numpy.zeros_like(f1)
        try:
            for i in range(1, nSteps):
                fn = "{2}{0}.outputl{1:0=5}.bin".format(simId+"Sim", i,self.extOutput)
                #SB.sid = simId +'Sim'+ ".{:0=5}".format(i)
                SB.readbinl(fn, verbose = False)
                Fcilinder1 = np.array([0,0,0])
                Fcilinder2 = Fcilinder1
                for j in range(0, n1-1):
                    r1[2] = SB.x[j, 2]
                    rpart = SB.x[j, :]
                    u = r1 - rpart
                    u = u/np.linalg.norm(u)
                    Ftotal = SB.force[j, :]
                    Frad = u*np.dot(Ftotal, u)
                    Fcilinder1 = Fcilinder1 + Frad
                for j in range(n1,n2):
                    r2[2] = SB.x[j, 2]
                    rpart = SB.x[j, :]
                    u = r2 - rpart
                    u = u/np.linalg.norm(u)
                    Ftotal = SB.force[j, :]
                    Frad = u*np.dot(Ftotal, u)
                    Fcilinder2 = Fcilinder2 + Frad
                # Write forces to matrices
                # f1[i, :] = Fcyl1[:]
                # f2[i, :] = Fcyl2[:]
                cylinder1.write(repr( Fcilinder1[0] ) + ' , ' + repr( Fcilinder1[1] ) + ' , ' +repr( Fcilinder1[2]) + '\n')
                cylinder2.write(repr( Fcilinder2[0] ) + ' , ' + repr( Fcilinder2[1] ) + ' , ' +repr( Fcilinder2[2]) + '\n')

                if i % 40 == 0:
                    print 'entro a fuerzas'
                    for x in range (0,nPart):
                        # fFile.write(repr(SB.force[x,0]) + ' , ' + repr(SB.force[x,1]) + ' , ' +repr(SB.force[x,2]) + '\n')
                        vFile.write(repr(SB.vel[x,0]) + ' , ' + repr(SB.vel[x,1]) + ' , ' +repr(SB.vel[x,2]) + '\n')
                        xFile.write(repr(SB.x[x,0]) + ' , ' + repr(SB.x[x,1]) + ' , ' +repr(SB.x[x,2]) + '\n')
        except RuntimeError:
            print "An error ocurred. Exiting writing loop and closing files."
        # sio.savemat(dataDir+'/cylinderForces.mat', dict(f1=f1, f2=2))
        # fFile.close()
        vFile.close()
        xFile.close()
        cylinder1.close()
        cylinder2.close()


    def velocityAveragesPerFrame(self,output='', simId='', ext=0):  # TESTEAR #Notdebugged
        box = self.boxDimensions*self.fcD
        nPoints = 2
        boxCenters = numpy.zeros([nPoints,3])
        lcube = self.boxDimensions[2]*self.fcD
        boxCenters[0,:]=box*1/2
        boxCenters[1,:]=[box[0]*.75,box[1]*1/2,box[2]*1/2]

        # This function sums velocity for each of the partcles.
        if self.extSave == 1:
            output = self.extOutput
        else:
            output = self.origOutput
        if simId == '':
            simId = self.simId

        dataDir = self.extOutput+simId

        if pa.isdir(dataDir) == 0:
            print 'Could not find {0} directory.'.format(simId)
            return
        vFile = open(dataDir+'/vAverage.csv', 'w')
        vTot = open(dataDir+'/vAverageTot.csv', 'w')
        # fFile = open(self.extOutput+'forces'+simId + '.csv','w')
        #nSteps=simulationInfo[0]
        nSteps = 400

        # Entering simulation
        SB = sphere.sim(np=2, sid=simId+'Sim')
        SB.extsave == 1
        #SB.readstep(1)

        n = np.size(SB.x, 0)
        m = np.size(SB.x, 1)
        averages = numpy.zeros([nSteps,nPoints])
        point = numpy.zeros([1,3])
        # Looping on every step
        try:
            for i in range(5000, nSteps):
                print i
                fn = "{2}{0}.outputl{1:0=5}.bin".format(simId+"Sim", i,self.extOutput)
                #SB.sid = simId +'Sim'+ ".{:0=5}".format(i)
                SB.readbinl(fn, verbose = False)
                # Write forces to matrices
                # f1[i, :] = Fcyl1[:]
                # f2[i, :] = Fcyl2[:]
                particles = SB.x
                for j in range(0,nPoints):
                    point[0,:]=boxCenters[j,:]
                    indices = gm.getIndices(particles,point,lcube)
                    velocities = SB.vel[indices,:]
                    magnitudes = np.sqrt(np.sum(numpy.power(velocities, 2), axis=1))
                    averages[i,j] = numpy.average(magnitudes)
                    if numpy.isnan(averages[i,j]):
                        averages[i,j]=0
                    vFile.write('{0}, '.format(averages[i,j]))
                vFile.write('\n')



        except RuntimeError:
            print "An error ocurred. Exiting writing loop and closing files."
        for i in range(0,nPoints):
            vTot.write('{0}, \n'.format(numpy.average(averages[:,i])))
        vTot.write('{0}'.format(boxCenters))
        # sio.savemat(dataDir+'/cylinderForces.mat', dict(f1=f1, f2=2))
        # fFile.close()
        vFile.close()
        vTot.close()
    def meanGrainVels(self,output='', simId='', ext=0):  # TESTEAR #Notdebugged
        box = self.boxDimensions*self.fcD
        nPoints = 2
        boxCenters = numpy.zeros([nPoints,3])
        lcube = np.ones([1,3])*self.boxDimensions[2]*self.fcD
        boxCenters[0,:]=box*1/2
        boxCenters[1,:]=[box[0]*.75,box[1]*1/2,box[2]*1/2]

        # This function sums velocity for each of the partcles.
        if self.extSave == 1:
            output = self.extOutput
        else:
            output = self.origOutput
        if simId == '':
            simId = self.simId

        dataDir = self.extOutput+simId

        if pa.isdir(dataDir) == 0:
            print 'Could not find {0} directory.'.format(simId)
            return
        vTot = open(dataDir+'/velocityAverage.txt', 'w')
        # fFile = open(self.extOutput+'forces'+simId + '.csv','w')
        #nSteps=simulationInfo[0]
        #conversionFactors = simulationInfo[2]
        nSteps = 5010
        conversionFactors = 1
        # Entering simulation
        SB = sphere.sim(np=2, sid=simId+'Sim')
        SB.extsave == 1
        #SB.readstep(1)

        n = np.size(SB.x, 0)
        m = np.size(SB.x, 1)
        #averages = numpy.zeros([nSteps, nPoints])
        point = numpy.zeros([1, 3])
        # Looping on every step
        particleCounter = np.zeros([1, nPoints])
        velocitySum = np.zeros([1, nPoints])
        nanCounter = np.zeros([1, nPoints])
        for i in range(5000, nSteps):
            print i
            fn = "{2}{0}.outputl{1:0=5}.bin".format(simId+"Sim", i,self.extOutput)
            #SB.sid = simId +'Sim'+ ".{:0=5}".format(i)
            SB.readbinl(fn, verbose=False)
            # Write forces to matrices
            # f1[i, :] = Fcyl1[:]
            # f2[i, :] = Fcyl2[:]
            particles = SB.x
            for j in range(0, nPoints):
                point[0, :] = boxCenters[j, :]
                indices = gm.getIndices(particles, point, lcube)
                velocities = SB.vel[indices, :]
                magnitudes = np.sqrt(np.sum(numpy.power(velocities, 2), axis=1))
                possibleNan = numpy.isnan(magnitudes)
                for index, value in enumerate(possibleNan):
                    if value == 1:
                        magnitudes[index] = 0
                        nanCounter[0, j] += 1
                        print 'NaN encountered in frame {0}'.format(i)

                velocitySum[0, j] += np.sum(magnitudes)
                particleCounter[0, j] += np.size(indices)

        averages = velocitySum/(particleCounter-nanCounter)

        if conversionFactors == 1:
            averages /= self.fcD
            velocitySum /= self.fcD

        averagesMat = (numpy.matlib.repmat(averages,nPoints,1))**2
        squaredDifference = averagesMat.T-averagesMat
        for i in range(0, nPoints):
            vTot.write('Point {0} \n Velocity sum {1}, \n N particles {2}, Total average {3} \n NaNs  {4} \n'.format(boxCenters[i, :], velocitySum[0, i], particleCounter[0, i], averages[0, i], nanCounter[0, i]))
        vTot.write('Squared difference between points. each element a(i,j) represents the difference  point j - point i. \n point order is the same order from the points vector with size nPoints*3 (input by the user). \n {0}\n'.format(squaredDifference))
        VTot.close()
    def velocityFluctuations(self,output='', simId='', ext=0):  # TESTEAR #Notdebugged
        box = self.boxDimensions*self.fcD
        nPoints = 2
        boxCenters = numpy.zeros([nPoints,3])
        lcube = np.ones([1,3])*self.boxDimensions[2]*self.fcD
        boxCenters[0,:]=box*1/2
        boxCenters[1,:]=[box[0]*.75,box[1]*1/2,box[2]*1/2]

        # This function sums velocity for each of the partcles.
        if self.extSave == 1:
            output = self.extOutput
        else:
            output = self.origOutput
        if simId == '':
            simId = self.simId

        dataDir = self.extOutput+simId

        if pa.isdir(dataDir) == 0:
            print 'Could not find {0} directory.'.format(simId)
            return
        file = open(dataDir+'/fluctuations.txt', 'w')
        pfile = open(dataDir+'/fluctuations.pickle', 'w')
        # fFile = open(self.extOutput+'forces'+simId + '.csv','w')
        #nSteps=simulationInfo[0]
        #conversionFactors = simulationInfo[2]
        nSteps = 5010
        conversionFactors = 1
        # Entering simulation
        SB = sphere.sim(np=2, sid=simId+'Sim')
        SB.extsave == 1
        #SB.readstep(1)

        n = np.size(SB.x, 0)
        m = np.size(SB.x, 1)
        #averages = numpy.zeros([nSteps, nPoints])
        point = numpy.zeros([1, 3])
        # Looping on every step
        particleCounter = np.zeros([1, nPoints])
        ensemble = []
        for i in range(0, nPoints): ensemble.append([])

        nanCounter = np.zeros([1, nPoints])
        fluctuations = numpy.zeros([nPoints,1])

        for i in range(5000, nSteps):
            print i
            fn = "{2}{0}.outputl{1:0=5}.bin".format(simId+"Sim", i,self.extOutput)
            #SB.sid = simId +'Sim'+ ".{:0=5}".format(i)
            SB.readbinl(fn, verbose=False)
            # Write forces to matrices
            # f1[i, :] = Fcyl1[:]
            # f2[i, :] = Fcyl2[:]
            particles = SB.x

            for j in range(0, nPoints):
                vHistory = ensemble[j]
                point[0, :] = boxCenters[j, :]
                indices = gm.getIndices(particles, point, lcube)
                velocities = SB.vel[indices, :]
                magnitudes = np.sqrt(np.sum(numpy.power(velocities, 2), axis=1))
                possibleNan = numpy.isnan(magnitudes)
                for index, value in enumerate(possibleNan):
                    if value == 1:
                        magnitudes[index] = 0
                        nanCounter[0, j] += 1
                        print 'NaN encountered in frame {0}'.format(i)
                vHistory.append(numpy.reshape(magnitudes,[numpy.size(magnitudes), 1]))
                particleCounter[0, j] += np.size(indices)
                # Updating
                ensemble[j] = vHistory

        for j in range(0,nPoints):
            vHistory = numpy.vstack(ensemble[j])
            if conversionFactors == 1:
                vHistory /= self.fcD
            fluctuations[j, 0] = numpy.var(vHistory)

        averagesVar = (numpy.matlib.repmat(fluctuations, 1, nPoints))
        squaredDifference = averagesVar-averagesVar.T
        for i in range(0, nPoints):
            file.write('Point {0} \n  fluctuations {1} \n'.format(boxCenters[i, :], fluctuations[j, 0]))
        file.write('Squared difference between fluctuations. each element a(i,j) represents the difference  point j - point i. \n point order is the same order from the points vector with size nPoints*3 (input by the user). \n {0}\n'.format(squaredDifference))
        file.write(' \n Pickle file contains values for fluctuations and its differences (in that order)')
        pickle.dump([fluctuations,squaredDifference], pfile)
        pfile.close()
        file.close()

    def porosityAverages(self,output='', simId='', resolution=15000,ext=0):  # TESTEAR #Notdebugged
        box = self.boxDimensions*self.fcD
        nPoints = 2
        boxCenters = numpy.zeros([nPoints,3])
        lcube = numpy.ones([1,3])*self.boxDimensions[2]*self.fcD*.8
        boxCenters[0,:]=box*1/2
        boxCenters[1,:]=[box[0]*.3,box[1]*1/2,box[2]*1/2]

        if self.extSave == 1:
            output = self.extOutput
        else:
            output = self.origOutput
        if simId == '':
            simId = self.simId

        dataDir = self.extOutput+simId

        if pa.isdir(dataDir) == 0:
            print 'Could not find {0} directory.'.format(simId)
            return
        pFile = open(dataDir+'/porosityEvolution.csv', 'w')
        pTot = open(dataDir+'/porosityAverageTot.csv', 'w')
        # fFile = open(self.extOutput+'forces'+simId + '.csv','w')
        #nSteps=simulationInfo[0]
        nSteps = 400

        # Entering simulation
        SB = sphere.sim(np=2, sid=simId+'Sim')
        SB.extsave == 1
        #SB.readstep(1)

        n = np.size(SB.x, 0)
        m = np.size(SB.x, 1)
        porosities = numpy.zeros([nSteps,nPoints])
        point = numpy.zeros([1,3])
        # Looping on every step
        try:
            for i in range(0, nSteps):
                print i
                fn = "{2}{0}.outputl{1:0=5}.bin".format(simId+"Sim", i*40,self.extOutput)
                #SB.sid = simId +'Sim'+ ".{:0=5}".format(i)
                SB.readbinl(fn, verbose = False)
                # Write forces to matrices
                # f1[i, :] = Fcyl1[:]
                # f2[i, :] = Fcyl2[:]
                particles = SB.x
                for j in range(0,nPoints):
                    point[0,:]=boxCenters[j,:]
                    indices = gm.getIndices(particles,point,lcube)
                    porosities[i,j] = gm.porosityCube(point,lcube,particles[indices,:],SB.radius[indices],resolution)
                    print porosities[i,j]
                    if numpy.isnan(porosities[i,j]):
                        porosities[i,j]=0
                    pFile.write('{0}, '.format(porosities[i,j]))
                pFile.write('\n')

        except RuntimeError:
            print "An error ocurred. Exiting writing loop and closing files."
        for i in range(0,nPoints):
            pTot.write('{0}, \n'.format(numpy.average(porosities[:,i])))
        pTot.write('{0}'.format(boxCenters))
        # sio.savemat(dataDir+'/cylinderForces.mat', dict(f1=f1, f2=2))
        # fFile.close()
        pFile.close()
        pTot.close()

    def regionAnalysis(self,output='', simId='', resolution=10000,ext=0,nom = ''):  # TESTEAR #Notdebugged

        # Inputs : Temporal while finishing function
        box = self.boxDimensions*self.fcD
        nPoints = 2
        boxCenters = numpy.zeros([nPoints,3])
        lcube = numpy.ones([1,3])*self.boxDimensions[2]*self.fcD*.8
        boxCenters[0,:]=box*1/2
        boxCenters[1,:]=[box[0]*.3,box[1]*1/2,box[2]*1/2]

        #   Summation variables
        pSum = np.zeros([nPoints, 1])
        pCounter = np.zeros([nPoints, 1])
        pNanCounter = np.zeros([nPoints, 1])
        vCounter = np.zeros([nPoints, 1])
        vNanCounter = np.zeros([nPoints, 1])
        velocityAverages = np.zeros([nPoints,1])
        # List that will contain velocity history for each point
        ensemble = []
        fluctuations = numpy.zeros([nPoints,1])

        for i in range(0, nPoints): ensemble.append([])

        # Check directories
        if self.extSave == 1:
            output = self.extOutput
        else:
            output = self.origOutput
        if simId == '':
            simId = self.simId

        dataDir = self.extOutput+simId

        if pa.isdir(dataDir) == 0:
            print 'Could not find {0} directory.'.format(simId)
            return

        # Create data files
        file = open(dataDir+'/regionAnalysis'+nom+'.txt', 'w')
        pFile = open(dataDir+'regionAnalysis'+nom+'.pickle','w')

        # Simulation data
        # nSteps=simulationInfo[0]
        nSteps = 20000
        conversionFactors = 1

        # Entering simulation
        SB = sphere.sim(np=2, sid=simId+'Sim')
        SB.extsave == 1

        n = np.size(SB.x, 0)
        m = np.size(SB.x, 1)
        porosities = numpy.zeros([nSteps,nPoints])
        point = numpy.zeros([1,3])
        # Looping on every step

        try:
            for i in range(5000, nSteps):
                print i

                fn = "{2}{0}.outputl{1:0=5}.bin".format(simId+"Sim", i,self.extOutput)
                #SB.sid = simId +'Sim'+ ".{:0=5}".format(i)
                SB.readbinl(fn, verbose = False)

                particles = SB.x

                for j in range(0, nPoints):
                    vHistory = ensemble[j]
                    point[0, :] = boxCenters[j, :]
                    indices = gm.getIndices(particles,point,lcube)
                    porosity = gm.porosityCube(point, lcube, particles[indices, :], SB.radius[indices], resolution)
                    velocities = SB.vel[indices, :]
                    magnitudes = np.sqrt(np.sum(numpy.power(velocities, 2), axis=1))
                    possibleNan = numpy.isnan(magnitudes)

                    for index, value in enumerate(possibleNan):
                        if value == 1:
                            magnitudes[index] = 0
                            vNanCounter[j, 0] += 1
                            print 'NaN encountered in frame {0} for velocities'.format(i)

                    if numpy.isnan(porosity):
                        porosity = 0
                        pNanCounter[j,0] += 1
                        print 'NaN encountered in frame {0} for porosities'.format(i)
                    vHistory.append(numpy.reshape(magnitudes,[numpy.size(magnitudes), 1]))
                    pSum[j,0] += porosity
                    pCounter[j,0] += 1
                    vCounter[j,0] += np.size(indices)
                    ensemble[j] = vHistory

        except RuntimeError:
            print "An error ocurred. Exiting writing loop and closing files."

        for j in range(0,nPoints):
            vHistory = numpy.vstack(ensemble[j])
            if conversionFactors == 1:
                vHistory /= self.fcD
            fluctuations[j, 0] = numpy.var(vHistory)
            velocityAverages[j,0] = numpy.mean(vHistory)

        porosityAverages = pSum/(pCounter-pNanCounter)

        # Data calc
        vSquared = velocityAverages**2
        vMat = numpy.matlib.repmat(vSquared, 1, nPoints)
        pvMat = numpy.matlib.repmat(vSquared*porosityAverages, 1, nPoints)
        varMat = (numpy.matlib.repmat(fluctuations, 1, nPoints))
        vSquaredDifference = vMat - vMat.T
        pvSquaredDifference = pvMat - pvMat.T
        varDifference = varMat - varMat.T
        # Saving data
        file.write('Region analysis results: \n')
        for i in range(0, nPoints):
            file.write('Point {0} \n center = {1}  \n average = {2}\n'.format(i, boxCenters[i, :], velocityAverages[i, 0]))
            file.write('porosity sum  {0} \n average = {1} \n'.format(pSum[i,0], porosityAverages[i,0]))
            file.write('Fluctuations {0}\n'.format(fluctuations))
        file.write('Squared differences between points. each element a(i,j) represents the difference  point j - point i. \n point order is the same order from the points vector with size nPoints*3 (input by the user). \n v - v\n {0} \n pv - pv\n {1} \n var - var\n{2}\n'.format(vSquaredDifference,pvSquaredDifference,varDifference))
        pickle.dump([velocityAverages, porosityAverages, fluctuations, vSquaredDifference, pvSquaredDifference,varDifference], pFile)

        file.close()
        pFile.close()

    def backupAndClean(self,ext=1,simId='ricardop04p01'):
        # Function that clears all .bin and .vtk files
        if simId == ' ':
            simId = self.simId
        if ext == 0:
            dataDir = self.origOutput
        if ext == 1:
            dataDir= self.extOutput
        os.chdir(dataDir)
        convs = glob.glob('*'+simId+'*-conv.log')
        bins = glob.glob('*'+simId+'*.bin')
        vtks = glob.glob('*'+simId+'*.vtk')

        if  numpy.size(convs) == 0:
            # Checking if dir already exists
            destination=dataDir+simId #CHANGE
            if pa.isdir(destination) == 0:
                os.mkdir(destination)
            # Copying files
            bins.sort()
            n = numpy.size(bins)
            m = numpy.size(vtks)

            shutil.copy(bins[0],destination)
            shutil.copy(bins[n-1],destination)
            #shutil.copy(convs[0],destination)

            # Deleting files
            #os.remove(convs[0])

            for i in range(0,n):
                os.remove(bins[i])
                #os.remove(vtks[i])
            print 'Succesfully removed files. Backup located at {0}'.format(destination)
        else:
            print('Not possible to back up the files. Number of conv files found is {0}.'.format(numpy.size(convs)))
            print('In case number of files is more than one, algorithm is detecting more than one simulation. please backup the one with longer ID first.')
        # ir a directorio

        self.extOutput

    def cylinderForces2txt(self, simId=' '):  # TESTEAR
        # This function sums the forces from each of the cylinders for every
        # in the simulation and saves them in a  '.m' format matrix.


        simId = 'ricardop06p01'

        dataDir = self.extOutput+simId
        if pa.isdir(dataDir) == 0:
            os.mkdir(dataDir)

        if pa.isdir(dataDir) == 0:
            print 'Could not find {0} directory.'.format(simId)
            return
        file = open(dataDir+'/cylinderInfo.obj')
        file2 = open(dataDir+'/simInfo.obj')
        cylinderInfo=pickle.load(file)
        simulationInfo=pickle.load(file2)
        file.close()
        file2.close()
        cylinder1 = open(dataDir+'/cylinder1Info.txt', 'w')
        cylinder2 = open(dataDir+'/cylinder2Info.txt', 'w')
        n1=cylinderInfo[0]
        n2=cylinderInfo[1]
        x1=cylinderInfo[2]
        x2=cylinderInfo[3]
        nSteps=20000
        # Entering simulation

        SB = sphere.sim(np=2, sid=simId+'Sim')
        SB.extsave == 1
        SB.readstep(1)
        # Looping on every step
        r1=numpy.array([x1[0],x1[1],0])*self.fcD
        r2=numpy.array([x2[0],x2[1],0])*self.fcD
        f1=numpy.zeros([nSteps,3])
        f2=numpy.zeros_like(f1)
        for i in range(1, nSteps):
            Fcilinder1 = np.array([0,0,0])
            Fcilinder2 = Fcilinder1
            for j in range(0, n1):
                r1[2] = SB.x[j, 2]
                rpart = SB.x[j, :]
                u = r1 - rpart
                u = u/np.linalg.norm(u)
                Ftotal = SB.force[j, :]
                Frad = u*np.dot(Ftotal, u)
                Fcilinder1 = Fcilinder1 + Frad
            for j in range(n1, n2):
                r2[2] = SB.x[j, 2]
                rpart = SB.x[j, :]
                u = r2 - rpart
                u = u/np.linalg.norm(u)
                Ftotal = SB.force[j, :]
                Frad = u*np.dot(Ftotal, u)
                Fcilinder2 = Fcilinder2 + Frad

            # Write forces to matrices
            # f1[i, :] = Fcyl1[:]
            # f2[i, :] = Fcyl2[:]
            cylinder1.write(repr( Fcilinder1[0] ) + ' , ' + repr( Fcilinder1[1] ) + ' , ' +repr( Fcilinder1[2]) + '\n')
            cylinder2.write(repr( Fcilinder2[0] ) + ' , ' + repr( Fcilinder2[1] ) + ' , ' +repr( Fcilinder2[2]) + '\n')
            SB.readstep(i*50)
        # sio.savemat(dataDir+'/cylinderForces.mat', dict(f1=f1, f2=2))
        cylinder1.close()
        cylinder2.close()

        def backupAndClean(self,ext=1,simId=' '):
            # Function that clears all .bin and .vtk files
            if simId == ' ':
                simId = self.simId
            if ext == 0:
                dataDir = self.origOutput
            if ext == 1:
                dataDir= self.extOutput
            os.chdir(dataDir)
            convs = glob.glob('*'+simId+'*-conv.log')
            bins = glob.glob('*'+simId+'*.bin')
            vtks = glob.glob('*'+simId+'*.vtk')

            if numpy.size(convs) == 1:
                # Checking if dir already exists
                destination=dataDir+simId #CHANGE
                if pa.isdir(destination) == 0:
                    os.mkdir(destination)
                # Copying files
                bins.sort()
                n = numpy.size(bins)
                m = numpy.size(vtks)
                shutil.copy(bins[0],destination)
                shutil.copy(bins[n-1],destination)
                shutil.copy(convs[0],destination)
                # Deleting files
                os.remove(convs[0])
                for i in range(0,n):
                    os.remove(bins[i])
                    os.remove(vtks[i])
                print 'Succesfully removed files. Backup located at {0}'.format(self.destination)
            else:
                print('Not possible to back up the files. Number of conv files found is {0}.'.format(numpy.size(convs)))
                print('In case number of files is more than one, algorithm is detecting more than one simulation. please backup the one with longer ID first.')
            # ir a directorio

            self.extOutput

    def porosityLine(self,simId='',nMeasurements=10,nPoints=10000,xo=0,xf=0,y=0,z=0,lcube=0, pickleFile = 'oneFrame.pickle'):
        if self.lowerWallsBool == 1 :
            box = self.boxDimensions*self.fcD
        else:
            box = self.boxDimensions

        if simId == '':
            simId = self.simId

        dataDir = self.extOutput+simId
        if pa.isdir(dataDir) == 0:
            os.mkdir(dataDir)

        file = open(dataDir+pickleFile)
        frameData = pickle.load(file)
        file.close()
        particles = frameData[0]
        radii = frameData[2]
        xo = box[0]/2
        xf = box*.8
        y = box[1]/2
        z = box[2]/2
        lcube = box[2]
        porosity = numpy.zeros([nMeasurements,1])
        L = (xf-xo-lcube)/n
        sweepPoints[0:n,0]=numpy.arange(1,n+1)*L+xo;
        sweepPoints[0:n,1]=zyCenter[0];
        sweepPoints[0:n,2]=zyCenter[1];
        porosityFile = open(dataDir+'porosityLine.csv','w')

        for i in range(0,nMeasurements):
        	print 'Finding particles for point',sweepPoints[i,:]
        	particlesCell,radiiCell,indices=gm.subCell(cellLengths,sweepPoints[i,:],particles,radii);
        	print 'Done. Particles found:', numpy.size(particlesCell[:,1]),'Calculating porosity.'
        	porosity[i]=gm.porosityCube(sweepPoints[i,:], lcube, particlesCell, radiiCell, nPoints);
        	porosityFile.write('{0}\n'.format(porosity[i]))

        porosityFile.close()

    ''' UNTESTED '''

    def twoCylinderExperiment(self, simId=' ', t=2, dt=.5e-6, dtFrame=.005, cylinderSeparation=.02, cylinderRadii=np.array([.0127, .0127]), nSides=30):
        if simId == ' ':
            simId = self.simId
        # Information printing
        cylinderPositions=np.array([[self.boxDimensions[0]/2, self.boxDimensions[1]/2], [self.boxDimensions[0]/2, self.boxDimensions[1]/2]])
        cylinderPositions[0, 1] = cylinderPositions[0, 1]+(cylinderSeparation+0.0254)/2
        cylinderPositions[1, 1] = cylinderPositions[1, 1]-(cylinderSeparation+0.0254)/2

        particles, radii, minlx, fixVelArr = gm.packYZ(self.boxDimensions[1], self.boxDimensions[2], self.r, self.N)
        print 'particles after packYZ', np.size(particles[:,0])
        if self.lowerWallsBool == 1 :
            print('Lower walls option is activated.')
            particles, radii, fixVelArr = gm.brakeWalls(radii, particles, fixVelArr, self.nWalls, self.lWalls, self.relativeRes, self.wallSepX, self.boxDimensions[0], self.boxDimensions[1], self.boxDimensions[2])
            print 'particles after walls', np.size(particles[:,0])
        else :
            print('Warning! Lower walls option NOT activaded. Setting lower walls off may cause simulation problems. Lower walls can be set using the setLowerWalls() function.')
        particles, radii, fixVelArr, n1 = gm.cylinderYZ(particles, radii, fixVelArr, cylinderRadii[0] , cylinderPositions[0,:], nSides, self.boxDimensions[2])
        particles, radii, fixVelArr, n2 = gm.cylinderYZ(particles, radii, fixVelArr, cylinderRadii[1] , cylinderPositions[1,:], nSides, self.boxDimensions[2])
        print 'particles after cylinders', np.size(particles[:,0])

        print 'Number of particles, particle radius and box dimensions', self.N, self.r, self.boxDimensions[:]
        print 'Initializing 2 cylinder experiment with following parameters'
        print 'Cylinder separation', cylinderSeparation
        print 'Cylinders radius', cylinderRadii[0], cylinderRadii[1]
        print 'Cylinder positions', cylinderPositions[0, :], cylinderPositions[1, :]
        print 'Cylinder indices limits 1 and 2',n1+n2, n1

        # Writing down experiment data
        first = n2
        firstx = cylinderPositions[1, :]
        second = n1+n2
        secondx = cylinderPositions[0, :]
        dataDir = self.extOutput+simId
        # CAMBIAR : convertir en funcion
        if pa.isdir(dataDir) == 0:
            os.mkdir(dataDir)
        file = open(dataDir+'/cylinderInfo.obj', 'w')
        file2 = open(dataDir+'/simInfo.obj', 'w')
        file3 = open(dataDir+'/simInfo.txt','w')

        file3.write('total N of particles (including walls and intruders) = {0} \n'.format(np.size(particles[:,0])))
        file3.write('Number of particles {0} \n Particle radius {1} \n Box dimensions {2} '.format(self.N, self.r, self.boxDimensions[:]))
        file3.write('Cylinder separation {0} \n'.format(cylinderSeparation))
        file3.write('Cylinders radius ({0}, {1}) \n'.format(cylinderRadii[0], cylinderRadii[1]))
        file3.write('Cylinder positions ({0},{1}) \n'.format(cylinderPositions[0, :], cylinderPositions[1, :]))
        file3.write('Cylinder index limits : {0},{1} \n'.format(n1+n2, n1))
        pickle.dump([first, second, firstx, secondx], file)  # index limit and positions for the cylinders
        pickle.dump([int(2*self.fcT/.2)], file2)  # nsteps
        file.close()
        file2.close()
        file3.close()
        # Sphere part
        SBB = sphere.sim(np=np.size(particles, axis=0), sid=simId)
        # Transfering position data to SB
        SBB.x = particles*self.fcD
        SBB.radius = radii*self.fcD
        SBB.fixvel = fixVelArr
        # Material properties
        SBB.g[0] = self.g*self.fcG
        SBB.g[1] = 0.0
        SBB.g[2] = 0.0
        SBB.mu_s[0] = self.muS
        SBB.mu_d[0] = self.muD
        SBB.k_n[0] = self.kN*self.fcK  # 535 #1.0e4 #la constante 535 corresponde a la usada por Ricardo en el codigo de Wassgren
        SBB.k_t[0] = self.kT*self.fcK
        SBB.gamma_n[0] = self.gam*16*self.fcGamma  # Esta constante se calculo a partir del coeficiente de restitucion
        V = ((SBB.radius[0])**3)*mt.pi*4/3

        SBB.rho[0] = (self.m*self.fcM)/V

        minlx = np.amax(particles[:, 0])
        if minlx < self.boxDimensions[0]:
            minlx = self.boxDimensions[0]
        # simulation initialization
        print 'dimensiones', minlx*1.1*self.fcD


        SBB.defineWorldBoundaries(L=np.array([minlx*1.1*self.fcD, self.boxDimensions[1]*self.fcD, self.boxDimensions[2]*1.1*self.fcD]), origo=[0, 0, 0], dx=2.1*self.r*self.fcD)
        SBB.initGrid()  # dx=ddx)
        SBB.initTemporal(total=.7*self.fcT, file_dt=dtFrame*self.fcT, dt=self.fcT*dt)

        # Using a 'dry' run, the sphere main program will display important parameters.
        # sphere will end after displaying these values.
        SBB.run(dry = True)
        # SBB.run()
        SBB.run()
        # SBB.writeVTKall()
        print('Simulation done')
        # Simulation
        # Create a sphere object with two preallocated particles and a simulation ID

        SB = sphere.sim(np=35899, sid=simId)
        SB.readstep(100)
        particles = SB.x[:, :]
        minlx = np.amax(SB.x[:, 0])
        print 'roof',self.boxDimensions[0]*self.fcD
        print 'radius', self.r*self.fcD
        print 'Highest particle', minlx
        print 'Number of particles in simulation', numpy.size(particles[:, 0])

        particlesIn = numpy.less(particles[:, 0], (self.boxDimensions[0]-self.r)*self.fcD)
        print 'number of particles in', numpy.size(particlesIn,0), particlesIn[0:100]
        if numpy.size(particlesIn) > 0:
            indices = [index for index, value in enumerate(particlesIn) if value == 1]
            # Do something
        print numpy.size(particles, 0)-numpy.size(indices, 0)

        # 10282

        lim = numpy.size(indices,0)
        newParticles = particles[indices[0:lim], :]
        newRadii = SB.radius[indices[0:lim]]
        newN = numpy.size(indices[0:lim])
        newfixVel = SB.fixvel[indices[0:lim]]

        # SPHERE

        # Create a sphere object with two preallocated particles and a simulation ID
        print 'entrando a segunda simulacion'
        print numpy.size(indices,0)
        sim = sphere.sim(np=newN, sid=simId+'Sim')

        # Transfering positions data to SB
        sim.x = newParticles
        sim.radius = newRadii
        sim.fixvel = newfixVel

        # Add gravitational acceleration
        sim.g[0] = self.g*self.fcG
        sim.g[1] = 0.0
        sim.g[2] = 0.0

        # SETTING MATERIAL PROPERTIES


        # CAMBIO
        sim.mu_s[0] = SB.mu_s[0]
        sim.mu_d[0] = SB.mu_d[0]
        sim.rho[0] = SB.rho[0]
        sim.k_n[0] = SB.k_n[0]
        sim.k_t[0] = SB.k_t[0]


        sim.periodicBoundariesX()
        minlx=np.amax(particles[:,0])







        sim.defineWorldBoundaries(L=np.array([self.boxDimensions[0]*self.fcD, self.boxDimensions[1]*self.fcD, self.boxDimensions[2]*1.1*self.fcD]), origo=[0, 0, 0], dx=(2.1*self.r*self.fcD))  # quitar el *4
        sim.initGrid()  # dx=ddx
        # Define the temporal parameters, e.g. the total time (total) and the file
        # output interval (file_dt), both in seconds
        sim.initTemporal(total=2*self.fcT, file_dt=.0043, dt=(0.5e-6)*self.fcT)

        # Using a 'dry' run, the sphere main program will display important parameters.
        # sphere will end after displaying these values.
        sim.run(dry=True)
        ### RUNNING THE SIMULATION
        # Start the simulation on the GPU from the sphere program
        #sim.run()
        sim.run()
        #sim.writeVTKall()
    def twoCylinderExperimentl(self, simId=' ', t=2, dt=.5e-6, dtFrame=.005, cylinderSeparation= None, cylinderRadii=np.array([.0127, .0127]), nSides=30):
        if cylinderSeparation == None:
            cylinderSeparation = self.cylinderSeparation
        if simId == ' ':
            simId = self.simId
        # Information printing
        cylinderPositions=np.array([[self.boxDimensions[0]/2, self.boxDimensions[1]/2], [self.boxDimensions[0]/2, self.boxDimensions[1]/2]])
        cylinderPositions[0, 1] = cylinderPositions[0, 1]+(cylinderSeparation+0.0254)/2
        cylinderPositions[1, 1] = cylinderPositions[1, 1]-(cylinderSeparation+0.0254)/2

        particles, radii, minlx, fixVelArr = gm.packYZ(self.boxDimensions[1], self.boxDimensions[2], self.r, self.N)
        print 'particles after packYZ', np.size(particles[:,0])
        if self.lowerWallsBool == 1 :
            print('Lower walls option is activated.')
            particles, radii, fixVelArr = gm.brakeWalls(radii, particles, fixVelArr, self.nWalls, self.lWalls, self.relativeRes, self.wallSepX, self.boxDimensions[0], self.boxDimensions[1], self.boxDimensions[2])
            print 'particles after walls', np.size(particles[:,0])
        else :
            print('Warning! Lower walls option NOT activaded. Setting lower walls off may cause simulation problems. Lower walls can be set using the setLowerWalls() function.')
        particles, radii, fixVelArr, n1 = gm.cylinderYZ(particles, radii, fixVelArr, cylinderRadii[0] , cylinderPositions[0,:], nSides, self.boxDimensions[2])
        particles, radii, fixVelArr, n2 = gm.cylinderYZ(particles, radii, fixVelArr, cylinderRadii[1] , cylinderPositions[1,:], nSides, self.boxDimensions[2])
        print 'particles after cylinders', np.size(particles[:,0])

        print 'Number of particles, particle radius and box dimensions', self.N, self.r, self.boxDimensions[:]
        print 'Initializing 2 cylinder experiment with following parameters'
        print 'Cylinder separation', cylinderSeparation
        print 'Cylinders radius', cylinderRadii[0], cylinderRadii[1]
        print 'Cylinder positions', cylinderPositions[0, :], cylinderPositions[1, :]
        print 'Cylinder indices limits 1 and 2',n1+n2, n1

        # Writing down experiment data
        first = n2
        firstx = cylinderPositions[1, :]
        second = n1+n2
        secondx = cylinderPositions[0, :]
        dataDir = self.extOutput+simId
        # CAMBIAR : convertir en funcion
        if pa.isdir(dataDir) == 0:
            os.mkdir(dataDir)
        file = open(dataDir+'/cylinderInfo.obj', 'w')
        file2 = open(dataDir+'/simInfo.obj', 'w')
        file3 = open(dataDir+'/simInfo.txt','w')

        file3.write('total N of particles (including walls and intruders) = {0} \n'.format(np.size(particles[:,0])))
        file3.write('Number of particles {0} \n Particle radius {1} \n Box dimensions {2} '.format(self.N, self.r, self.boxDimensions[:]))
        file3.write('Cylinder separation {0} \n'.format(cylinderSeparation))
        file3.write('Cylinders radius ({0}, {1}) \n'.format(cylinderRadii[0], cylinderRadii[1]))
        file3.write('Cylinder positions ({0},{1}) \n'.format(cylinderPositions[0, :], cylinderPositions[1, :]))
        file3.write('Cylinder index limits : {0},{1} \n'.format(n1+n2, n1))
        pickle.dump([first, second, firstx, secondx], file)  # index limit and positions for the cylinders
        pickle.dump([int(2*self.fcT/.2), np.size(particles[:, 0]), self.conversionFactorsBool], file2)  # nsteps, #nParticles,
        file.close()
        file2.close()
        file3.close()
        # Sphere part
        SBB = sphere.sim(np=np.size(particles, axis=0), sid=simId)
        # Transfering position data to SB
        SBB.x = particles*self.fcD
        SBB.radius = radii*self.fcD
        SBB.fixvel = fixVelArr
        # Material properties
        SBB.g[0] = self.g*self.fcG
        SBB.g[1] = 0.0
        SBB.g[2] = 0.0
        SBB.mu_s[0] = self.muS
        SBB.mu_d[0] = self.muD
        SBB.k_n[0] = self.kN*self.fcK  # 535 #1.0e4 #la constante 535 corresponde a la usada por Ricardo en el codigo de Wassgren
        SBB.k_t[0] = self.kT*self.fcK
        SBB.gamma_n[0] = self.gam*16*self.fcGamma  # Esta constante se calculo a partir del coeficiente de restitucion
        V = ((SBB.radius[0])**3)*mt.pi*4/3

        SBB.rho[0] = (self.m*self.fcM)/V

        minlx = np.amax(particles[:, 0])
        if minlx < self.boxDimensions[0]:
            minlx = self.boxDimensions[0]
        # simulation initialization
        print 'dimensiones', minlx*1.1*self.fcD


        SBB.defineWorldBoundaries(L=np.array([minlx*1.1*self.fcD, self.boxDimensions[1]*self.fcD, self.boxDimensions[2]*1.1*self.fcD]), origo=[0, 0, 0], dx=2.1*self.r*self.fcD)
        SBB.initGrid()  # dx=ddx)
        SBB.initTemporal(total=.7*self.fcT, file_dt=dtFrame*self.fcT, dt=self.fcT*dt)

        # Using a 'dry' run, the sphere main program will display important parameters.
        # sphere will end after displaying these values.
        SBB.run(dry = True)
        SBB.run()
        #SBB.run()
        # SBB.writeVTKall()
        print('Simulation done')
        # Simulation
        # Create a sphere object with two preallocated particles and a simulation ID
        frame = 100
        SB = sphere.sim(np=35899, sid=simId)
        fn = "{2}{0}.outputl{1:0=5}.bin".format(simId, frame,self.extOutput)
        SB.sid = simId +'Sim'+ ".{:0=5}".format(frame)
        SB.readbinl(fn, verbose = False)
        particles = SB.x[:, :]
        minlx = np.amax(SB.x[:, 0])
        '''
        print 'roof',self.boxDimensions[0]*self.fcD
        print 'radius', self.r*self.fcD
        print 'Highest particle', minlx
        print 'Number of particles in simulation', numpy.size(particles[:, 0])
        '''

        particlesIn = numpy.less(particles[:, 0], (self.boxDimensions[0]-self.r)*self.fcD)
        if numpy.size(particlesIn) > 0:
            indices = [index for index, value in enumerate(particlesIn) if value == 1]
            # Do something
        print numpy.size(particles, 0)-numpy.size(indices, 0)

        # 10282

        lim = numpy.size(indices,0)
        newParticles = particles[indices[0:lim], :]
        newRadii = SB.radius[indices[0:lim]]
        newN = numpy.size(indices[0:lim])
        newfixVel = SB.fixvel[indices[0:lim]]

        # SPHERE

        # Create a sphere object with two preallocated particles and a simulation ID
        print 'entrando a segunda simulacion'
        print numpy.size(indices,0)
        sim = sphere.sim(np=newN, sid=simId+'Sim')

        # Transfering positions data to SB
        sim.x = newParticles
        sim.radius = newRadii
        sim.fixvel = newfixVel

        # Add gravitational acceleration
        sim.g[0] = self.g*self.fcG
        sim.g[1] = 0.0
        sim.g[2] = 0.0

        # SETTING MATERIAL PROPERTIES


        # CAMBIO
        sim.mu_s[0] = SBB.mu_s[0]
        sim.mu_d[0] = SBB.mu_d[0]
        sim.rho[0] = SBB.rho[0]
        sim.k_n[0] = SBB.k_n[0]
        sim.k_t[0] = SBB.k_t[0]


        sim.periodicBoundariesX()
        minlx=np.amax(particles[:,0])







        sim.defineWorldBoundaries(L=np.array([self.boxDimensions[0]*self.fcD, self.boxDimensions[1]*self.fcD, self.boxDimensions[2]*1.1*self.fcD]), origo=[0, 0, 0], dx=(2.1*self.r*self.fcD))  # quitar el *4
        sim.initGrid()  # dx=ddx
        # Define the temporal parameters, e.g. the total time (total) and the file
        # output interval (file_dt), both in seconds
        sim.initTemporal(total=2*self.fcT, file_dt=.0043, dt=(0.5e-6)*self.fcT)

        # Using a 'dry' run, the sphere main program will display important parameters.
        # sphere will end after displaying these values.
        sim.run(dry=True)
        ### RUNNING THE SIMULATION
        # Start the simulation on the GPU from the sphere program
        #sim.run()
        sim.run()
        #sim.writeVTKall()
