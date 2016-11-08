#FUNCTIONS CREATED BY GRANULAR MATERIALS TEAM
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections
matplotlib.rcParams.update({'font.size': 7, 'font.family': 'serif'})
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
from matplotlib.font_manager import FontProperties
import subprocess
import pickle as pl
import random
import math
import numpy as np
import math as mt


'''
granMat BETA.
Authors: Eli cuellar, Leonardo Clemente

Function library made to aid in the research of granular Materials using the Sphere library.

Important data vectors:

Sphere objects have different property variables that must be filled prior to running a simulation.
Most of them being initial conditions and metadata for the simulation. granmat keeps these property
 values on variables of its own in case these are changed prior to initializing the simulation

fixedVelArr:
	fixedVelArr is an Int uni-dimensional array that helps indicate
	which of the particles in the 'particles' array have a fixed velocity and which dont
	(e.g. static obstacles have fixed velocities of 0 in all directions). vector is one of the many
	inputs for the initilization of a sphere simulation.

	particle cases
	-1: Fix horizontal and vertical velocities
	0: No fixe
	1: Fix vertical velocities

radii:

Particles: variable containing information about particles position. all types of particles are

IMPORTANT PREMISES/CONVENTIONS:

	1.-Fixed particles location within arrays are always at the start: Whenever a new obstacle gets inserted within
	the volume, positions and radii for the fixed particles get concatenated at start of the array.
	2.- BETA version doesn't include boundary awareness of  when inserting obstacles. If an obstacle gets inserted outside
	previously set boundaries, sphere algorithm will not run.




Please read sphere documentation for better understanding.

'''



def packYZ(ly,lz,r,N):

	'''
	Fills a rectangular volume with dimensions lx*ly*lz  and  N particles of radius r. lx is defined
	after determining the minimum neccesary length to fill up the volume. ly, lz, r and N are scalars.
	vector implementation for different radius yet to be defined. METHOD MAY VARY.
	r=0.0052
	ly=dx*19
	lz=dx*10
	N=10000
	'''

	#WARNINGS: function wont run if:
		#number of particles is smaller than maximum number of particles that could be contained in one floor.
	D=r*2.15
	Dred=D*.5
	Daug=D*1.1
	dx=r*2.1
	c=0;
	s= 'initParticles'

	#Generating radii
	radii=numpy.zeros(N).reshape(N,1);
	radii=radii+r;

	#Detecting number or particles per floor
	ny=int(math.floor(ly/D));
	nz=int(math.floor(lz/D));
	piso=ny*nz;
	nPisos=int(math.ceil(float(N)/float(piso)));
	rems=N%piso;
	#Generating particle block
	bloque=numpy.zeros((piso,3));

	for i in range(0, ny):
		for j in range(0,nz):
			bloque[c,:]=(0,i*D,j*D);
			c=c+1;
	bloque=bloque+r;

	#Filling up container

	bloqueTemp=numpy.array(bloque)
	bloqueTemp[:,1]=bloque[:,1]+random.random()*Dred
	bloqueTemp[:,2]=bloque[:,2]+random.random()*Dred
	partCont=numpy.array(bloque)

	for k in range(1,nPisos):
		bloqueTemp=numpy.array(bloque)
		bloqueTemp[:,1]=bloque[:,1]+random.random()*Dred
		bloqueTemp[:,2]=bloque[:,2]+random.random()*Dred
		bloqueTemp[:,0]=bloque[:,0]+k*Daug
		partCont=numpy.concatenate((partCont,bloqueTemp),axis=0)

	partCont=partCont[0:N,:];


	#Extra output values
	lx=math.ceil((partCont[N-1,0]+r)/dx)*dx;
	fixedVelArr=numpy.zeros(N).reshape(N,1);

	return partCont, radii, lx, fixedVelArr

def cylinderYZold(particles,radii,fixVelArr,radCyl,c,nSides,lz):
	'''
	#inserts cylinder obstacle within granular material and displaces particles within the obstacle
	geometry to the top of the container. obstacle extends towards z direction with geometry center
	c (cx,cy) and has radius radCyl.  Cylinders are approximated using a regular polygon of 'nSides'
	number of sides with spheres of radius 'radSph' centered on each of the polygon corners. The
	regular polygon is concentric with the cylinder. the spheres circumference make
	contact with the cylinder surface and the neightbour spheres only in one 'point' of space.

		#input parameters
		particles
		radii
		fixVelArr


 	#WARNINGS: If the particles simulating the obstacle are outside boundaries, sphere won't run the simulation.
	radCyl=dx*3
	c=numpy.array([lx/2 ,ly/2]) #center of cylinder
	nSides=15;
	lz=lz;'''

	#WARNINGS: function wont run if:
		#number of particles is smaller than maximum number of particles that could be contained in one floor.


	maxrad=numpy.amax(radii,axis=0)*1.4;

	#determining radius of circle and spheres
	radCir=radCyl/(1+math.sin(math.pi/nSides));
	radSph=radCyl-radCir;
	rp=radCir;


	#Finding particles that may be within the intruders
	ixp=numpy.less(particles[:,0],c[0]+rp+maxrad);
	ixn=numpy.greater(particles[:,0],c[0]-rp-maxrad);
	iyp=numpy.less(particles[:,1],c[1]+rp+maxrad);
	iyn=numpy.greater(particles[:,1],c[1]-rp-maxrad);

	particlesIn=ixp*ixn*iyp*iyn;


	#Generating coordinates
	s=numpy.array(numpy.arange(0,nSides)*(2*math.pi/nSides)); #array with Angle values for each of the
	xCord=radCir*numpy.cos(s).reshape(nSides,1);
	yCord=radCir*numpy.sin(s).reshape(nSides,1);
	zCord=numpy.zeros(nSides).reshape(nSides,1);
	circleCord=numpy.concatenate((xCord,yCord,zCord),axis=1);


	#generating cylinder
	nCuts=lz/(radSph*2);
	nCutsInt=int(math.floor(nCuts));
	freeSpace=(nCuts-nCutsInt)*radSph*2;
	cylinderCord=numpy.array(circleCord);

	for i in range(1,nCutsInt):
		tempBlock=circleCord;
		tempBlock[:,2]=i*radSph*2;
		cylinderCord=numpy.concatenate((cylinderCord,tempBlock),axis=0);


	#Placing cylinder
	cylinderCord[:,0]=cylinderCord[:,0]+c[0];
	cylinderCord[:,1]=cylinderCord[:,1]+c[1];
	cylinderCord[:,2]=cylinderCord[:,2]+freeSpace/2;



	#generating new radii and fixedVelArr
	nCyl=nSides*nCutsInt;
	newRadii=numpy.zeros(nCyl).reshape(nCyl,1);
	newFixVelArr=newRadii;
	newParticlesIn=numpy.zeros(nCyl, dtype=bool);


	newRadii=newRadii+radSph;
	newFixVelArr=newFixVelArr-1;

	#packing all particles together
	particles=numpy.concatenate((cylinderCord,particles),axis=0);
	radii=numpy.concatenate((newRadii,radii),axis=0);
	fixVelArr=numpy.concatenate((newFixVelArr,fixVelArr),axis=0);
	particlesIn=numpy.concatenate((newParticlesIn,particlesIn),axis=0);

	#computing distance of displacement if there is any
	if numpy.size(particles[particlesIn,0],axis=0)> 0:
		highest=numpy.amax(particles[:,0],axis=0);
		lowest=numpy.amin(particles[particlesIn,0],axis=0);
		L=highest-lowest+maxrad*1.01;
		particles[particlesIn,0]=particles[particlesIn,0]+L*1.1;


	return particles,radii,fixVelArr,nCyl


# PURPOSE: to determine weather a point falls inside a sphere. Returns a bool.


def IsInSphere(xp, xs, r):
    if(numpy.linalg.norm(xp-xs) < r):
        return 1
    else:
        return 0


# PURPOSE: To compute the porosity (% volume occupied) of a cubic cell, given a
# series of spheres using montecarlo method.


def porosityCube(xcube, lcube, xspheres, rspheres, Npoints):

    xp = numpy.random.rand(Npoints, 3)
    Nspheres = xspheres.shape[0]

    xlow = float(xcube[0,0]) - float(lcube[0,0]/2)
    xhi = float(xcube[0,0]) + float(lcube[0,0]/2)

    ylow = float(xcube[0,1]) - float(lcube[0,1]/2)
    yhi = float(xcube[0,1]) + float(lcube[0,1]/2)

    zlow = float(xcube[0,2]) - float(lcube[0,2]/2)
    zhi = float(xcube[0,2]) + float(lcube[0,2]/2)
    xp[:, 0] = xp[:, 0]*(xhi-xlow)+xlow
    xp[:, 1] = xp[:, 1]*(yhi-ylow)+ylow
    xp[:, 2] = xp[:, 2]*(zhi-zlow)+zlow

    Ninside = 0

    # For each point, determine wether it falls inside a sphere
    for i in range(0, Npoints):
        ithp = xp[i, :]
        for j in range(0, Nspheres):
            if IsInSphere(ithp, xspheres[j, :], rspheres[j]):
                Ninside = Ninside + 1
                break

    # Return porosity
    return float(Ninside)/float(Npoints)


def subCell(cellLengths,cellPos,particles,radii):

	#Porosity calculation uses a 3D square cell that sweeps through space

	maxRad=numpy.amax(radii);
	c=maxRad+cellLengths/2;

		#Finding particles
	ixp=numpy.less(particles[:,0],cellPos[0]+c[0]);
	ixn=numpy.greater(particles[:,0],cellPos[0]-c[0]);
	iyp=numpy.less(particles[:,1],cellPos[1]+c[1]);
	iyn=numpy.greater(particles[:,1],cellPos[1]-c[1]);
	izp=numpy.less(particles[:,2],cellPos[2]+c[2]);
	izn=numpy.greater(particles[:,2],cellPos[2]-c[2])
	particlesIn=ixp*ixn*iyp*iyn*izp*izn;

	if numpy.size(particlesIn) > 0:
		indices=[index for index,value in enumerate(particlesIn) if value == 1];
		#Do something
	return particles[particlesIn,:],radii[particlesIn],indices

def cylinderYZ(particles,radii,fixVelArr,radCyl,c,nSides,lz):
	'''
	#inserts cylinder obstacle within granular material and displaces particles within the obstacle
	geometry to the top of the container. obstacle extends towards z direction with geometry center
	c (cx,cy) and has radius radCyl.  Cylinders are approximated using a regular polygon of 'nSides'
	number of sides with spheres of radius 'radSph' centered on each of the polygon corners. The
	regular polygon is concentric with the cylinder. the spheres circumference make
	contact with the cylinder surface and the neightbour spheres only in one 'point' of space.

		#input parameters
		particles
		radii
		fixVelArr


 	#WARNINGS: If the particles simulating the obstacle are outside boundaries, sphere won't run the simulation.
	radCyl=dx*3
	c=numpy.array([lx/2 ,ly/2]) #center of cylinder
	nSides=15;
	lz=lz;'''

	#WARNINGS: function wont run if:
		#number of particles is smaller than maximum number of particles that could be contained in one floor.


	maxrad=numpy.amax(radii,axis=0)*1.3;

	#determining radius of circle and spheres
	radCir=radCyl/(1+math.sin(math.pi/nSides));
	radSph=radCyl-radCir;
	rp=radCir;


	#Finding particles that may be within the intruders
	ixp=numpy.less(particles[:,0],c[0]+rp+maxrad);
	ixn=numpy.greater(particles[:,0],c[0]-rp-maxrad);
	iyp=numpy.less(particles[:,1],c[1]+rp+maxrad);
	iyn=numpy.greater(particles[:,1],c[1]-rp-maxrad);

	particlesIn=ixp*ixn*iyp*iyn;


	#Generating coordinates
	s=numpy.array(numpy.arange(0,nSides)*(2*math.pi/nSides)); #array with Angle values for each of the
	xCord=radCir*numpy.cos(s).reshape(nSides,1);
	yCord=radCir*numpy.sin(s).reshape(nSides,1);
	zCord=numpy.zeros(nSides).reshape(nSides,1);
	circleCord=numpy.concatenate((xCord,yCord,zCord),axis=1);


	#generating cylinder
	nCuts=lz/(radSph*2);
	nCutsInt=int(math.floor(nCuts));
	freeSpace=(nCuts-nCutsInt)*radSph*2;
	cylinderCord=numpy.zeros([0,3]);
	spaceDiv=freeSpace/(nCuts+1)

	for i in range(0,nCutsInt):
		tempBlock=circleCord;
		tempBlock[:,2]=i*(radSph*2+spaceDiv)+spaceDiv+radSph;
		cylinderCord=numpy.concatenate((cylinderCord,tempBlock),axis=0);

	#Placing cylinder
	cylinderCord[:,0]=cylinderCord[:,0]+c[0];
	cylinderCord[:,1]=cylinderCord[:,1]+c[1];
	cylinderCord[:,2]=cylinderCord[:,2];



	#generating new radii and fixedVelArr
	nCyl=nSides*nCutsInt;
	newRadii=numpy.zeros(nCyl).reshape(nCyl,1);
	newFixVelArr=newRadii;
	newParticlesIn=numpy.zeros(nCyl, dtype=bool);


	newRadii=newRadii+radSph;
	newFixVelArr=newFixVelArr-1;

	#packing all particles together
	particles=numpy.concatenate((cylinderCord,particles),axis=0);
	radii=numpy.concatenate((newRadii,radii),axis=0);
	fixVelArr=numpy.concatenate((newFixVelArr,fixVelArr),axis=0);
	particlesIn=numpy.concatenate((newParticlesIn,particlesIn),axis=0);

	#computing distance of displacement if there is any
	if numpy.size(particles[particlesIn,0],axis=0)> 0:
		highest=numpy.amax(particles[:,0],axis=0);
		lowest=numpy.amin(particles[particlesIn,0],axis=0);
		L=highest-lowest+maxrad*1.01;
		particles[particlesIn,0]=particles[particlesIn,0]+L*1.1;

	print "Number of particles per cylinder:", nCyl
	return particles,radii,fixVelArr,nCyl


def brakeWalls(radii,particles,fixVelArr,nWalls,lWalls,relativeRes,wallSepX,lx, ly, lz):

	'''
	radii = vector containing information about particles radius
	particles = vector containing particle positions
	fixVelArr = information about particles kind (static/moving)
	nWalls = Number of obstacle walls
	lWalls = length of walls
	relativeRes = sets wall particles size relative to the smallest moving particle in vector particles. e.g. if relativeRes= 		5, particles from walls will be (approximately) 5 times smaller than the smallest moving particle.
	wallSepX = separation from walls at X dimension (heigth)
	lx,ly,lz = container dimensions

	'''


	r=np.amin(radii)
	# Encontrando diametro de esfera
	rRes=(r/relativeRes)
	divisor=lWalls/rRes
	divisorInt=mt.floor(divisor)
	rem=divisor-divisorInt
	if rem>0 :
		diametroExtra=rem*rRes/divisorInt;
		rRes=rRes+diametroExtra

	wallSphD=rRes

	ny=round(lWalls/wallSphD)
	nz=lz/wallSphD
	nzInt=int(mt.floor(nz))
	rem=nz-nzInt


	if rem > 0 :
		dSepZ=(rem*wallSphD)/(nzInt+1) #Distancia que hay entre cada una de las esferas en Z
	else :
		dSepZ=0

	nSph=int(nzInt*ny)
	block=np.zeros((nSph,3))
	tempBlock=block
	row=np.zeros((ny,3))
	for i in range(0,int(ny)):
		row[i,:]=(0,i*wallSphD,0)

	for i in range(0,nzInt):
		row[:,2]=i*(wallSphD)+(i+1)*dSepZ

		block[i*ny:(i+1)*ny,:]=row

	wallSepY = (ly - nWalls*lWalls)/(nWalls+1)
	block[:,1]=block[:,1]+wallSphD/2

	wallsCoord=np.zeros((nWalls*nSph,3))
	newFixVelArr=np.zeros((nWalls*nSph,1))
	newFixVelArr[:,0]=-1;
	newRadii=np.zeros((nWalls*nSph,1))
	newRadii[:,0]=wallSphD;

	l1=wallSepX+wallSphD
	l2=wallSphD/2

	for i in range(0,nWalls) :
		tempBlock=np.zeros((nSph,3))
		tempBlock[:,1]=block[:,1]+(i+1)*(wallSepY)+i*(lWalls)
		if i%2 == 0 :
			tempBlock[:,0]=l1
		else :
			tempBlock[:,0]=l2
		tempBlock[:,2]=block[:,2]
		wallsCoord[nSph*i:nSph*(i+1),:]=tempBlock;
		del tempBlock

	# Translating particles
	L=wallSepX+2*wallSphD
	particles[:,0]=particles[:,0]+L

	# Packing all particles together
	particles=np.concatenate((wallsCoord,particles),axis=0);
	radii=np.concatenate((newRadii,radii),axis=0);
	fixVelArr=np.concatenate((newFixVelArr,fixVelArr),axis=0);

	return particles, radii, fixVelArr


def getIndices(vector,point,lcube): # TESTEAR
    '''Function that substracts values from a 3D vector 'vector' based on limits
    defined by the user in the vector 'box'.'''
    half = lcube/2
    up = numpy.less(vector[:, 0], point[0,0]+half[0,0])
    down = numpy.greater(vector[:, 0], point[0,0]-half[0,0])
    right = numpy.less(vector[:, 1], point[0,1]+half[0,1])
    left = numpy.greater(vector[:, 1], point[0,1]-half[0,1])
    front = numpy.less(vector[:, 2], point[0,2]+half[0,1])
    back = numpy.greater(vector[:, 2], point[0,2]-half[0,1])
    particlesIn = up*down*left*right*front*back

    if numpy.size(particlesIn) > 0:
        indices = [index for index, value in enumerate(particlesIn) if value == 1]
    else:
        indices = numpy.zeros([1,3])

    return indices
