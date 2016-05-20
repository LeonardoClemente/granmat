import numpy as np
import math as mt
import granmat as gm
import sphere
import fractions

# Parametros de experimento
nExperimentos = 1
nameExperiments = 'exp1'
N = 32500
# Parametros particula
r = .0052/2
dx = 2.1*r
# dimensiones contenedor
lx = int(.75/dx)*dx
ly = int(.30/dx)*dx
lz = 4*dx

# Factores de conversion
V = (r**3)*mt.pi*4/3
rho = 14
m = V*rho
fcK = .0052/(m*9.8)
fcGamma = mt.sqrt(.0052*9.8)
fcD = 1/0052
fcF = 1/(9.8*m)
fcT = mt.sqrt(9.8/.0052)

'''
INICIO DE FUNCION
Generando paredes

crea primer bloque de paredes. Las paredes estan hechas de esferas con un
diametro 'WallSph' relativo al radio minimo de las particulas a simular en
Sphere y determinado por el termino 'relativeRes'. relativeRes es un parametro
que indica en aproximacion la proporcion de tamano entre las particulas que
simulan la pared y las particulas con el radio minimo. Una cantidad de 10
indica una relacion aproximada de 10:1.
'''

nWalls = 9
lWalls = 0.005
relativeRes = 5
wallSepX = .052

'''
FUERA DEL CICLO GENERAR LA PRIMERA SIMULACION QUE DEPOSITA LAS PARTICULAS EN UN
CONTENEDOR SIN OBSTACULOS.
PARAMETROS SPHERE.
TIEMPO NECESARIO PARA ESTABILIZACION: 1 segundo?
'''

particles, radii, minlx, fixVelArr = gm.packYZ(ly, lz, r, N)


# Encontrando diametro de esfera
rRes = r/relativeRes
divisor = lWalls/rRes
divisorInt = mt.floor(divisor)
rem = divisor-divisorInt
if rem > 0:
	diametroExtra = rem*rRes/divisorInt
	rRes = rRes+diametroExtra

wallSphD = rRes



#generando bloque
ny=round(lWalls/wallSphD)
nz=lz/wallSphD
nzInt=int(mt.floor(nz))
rem=nz-nzInt


if rem > 0 :
	dSepZ=(rem*wallSphD)/(nzInt+1) #Distancia que hay entre cada una de las esferas en Z
else :
	dSepZ=0




#Generating particle block
nSph=int(nzInt*ny)
block=np.zeros((nSph,3))
tempBlock=block
row=np.zeros((ny,3))
for i in range(0,int(ny)):
	row[i,:]=(0,i*wallSphD,0)




for i in range(0,nzInt):
	row[:,2]=i*(wallSphD)+(i+1)*dSepZ

	block[i*ny:(i+1)*ny,:]=row

print block
print 'DEBUG Z', lz, dSepZ,nzInt,nz
print 'DEBUG Z', nz*wallSphD+dSepZ*nzInt,lz


wallSepY = (ly - nWalls*lWalls)/(nWalls+1)
block[:,1]=block[:,1]+wallSphD/2

wallsCoord=np.zeros((nWalls*nSph,3))
newFixVelArr=np.zeros((nWalls*nSph,1))
newFixVelArr[:,0]=-1;
newRadii=np.zeros((nWalls*nSph,1))
newRadii[:,0]=wallSphD;

l1=wallSepX+wallSphD
l2=wallSphD/2

#print 'QUEPEDO',block


for i in range(0,nWalls) :
	tempBlock=np.zeros((nSph,3))
	tempBlock[:,1]=block[:,1]+(i+1)*(wallSepY)+i*(lWalls)
	#print 'AQUI ESTA EL BEBE',tempBlock[:,1]
	if i%2 == 0 :
		tempBlock[:,0]=l1
	else :
		tempBlock[:,0]=l2
	tempBlock[:,2]=block[:,2]
	wallsCoord[nSph*i:nSph*(i+1),:]=tempBlock;
	del tempBlock

#translating particles
L=wallSepX+2*wallSphD
particles[:,0]=particles[:,0]+L

#packing all particles together
particles=np.concatenate((wallsCoord,particles),axis=0);
radii=np.concatenate((newRadii,radii),axis=0);
fixVelArr=np.concatenate((newFixVelArr,fixVelArr),axis=0);



''' INSERTING CYLINDERS'''

#Parametros cilindro
nIntruders=np.array([2]);
radCylVec=np.array([.0254,.0254]);
cVec=np.array([[lx/2,ly/2+.0254+.0125],[lx/2,ly/2-.0254-.0125]]);
nSides=30;

#for x in range(0,nIntruders): #MODIFIY counter x when added to experiment cycle
particles,radii,fixVelArr,n1 = gm.cylinderYZnew(particles,radii,fixVelArr,radCylVec[0],cVec[0,:],nSides,lz);
particles,radii,fixVelArr,n2 = gm.cylinderYZnew(particles,radii,fixVelArr,radCylVec[1],cVec[1,:],nSides,lz);



'''SPHERE'''

# Create a sphere object with two preallocated particles and a simulation ID
simID = 'factoresNuevos'
SBB = sphere.sim(np = np.size(particles,axis=0), sid = simID);

#Transfering positions data to SB
SBB.x=particles*fcD;
SBB.radius=radii*fcD;
SBB.fixvel=fixVelArr;


# Add gravitational acceleration
SBB.g[0] = -9.8/9.8
SBB.g[1] = 0.0
SBB.g[2] = 0.0


### SETTING MATERIAL PROPERTIES

SBB.mu_s[0] = 0.4
SBB.mu_d[0] = 0.4
SBB.rho[0] = 14
SBB.k_n[0] = 4084 #535 #1.0e4 #la constante 535 corresponde a la usada por Ricardo en el codigo de Wassgren
SBB.k_t[0] = 4084


#calculando gamma
eps=.5
mEff=(m**2)/(2*m)
s=((.5*np.log(eps)*mt.pi)**2+1)/mEff

gam=(SBB.k_n[0]/s)**1/2
#0.00385218831922*16
SBB.gamma_n[0] = gam*16*fcGamma #Esta constante se calculo a partir del coeficiente de restitucion
#CAMBIO
SBB.mu_s[0] = 0.4
SBB.mu_d[0] = 0.4
V=((SBB.radius[0]*fcD)**3)*mt.pi*4/3

m=1
SBB.rho[0] = m/V
SBB.k_n[0] = 4084*fcK #535 #1.0e4 #la constante 535 corresponde a la usada por Ricardo en el codigo de Wassgren
SBB.k_t[0] = 4084*fcK

minlx=np.amax(particles[:,0])

if minlx < lx :
  minlx=lx

ddx= 2.1*r
SBB.defineWorldBoundaries(L=np.array([minlx*1.1*fcD,ly*fcD,lz*1.1*fcD]),origo = [0,0,0], dx=ddx) #quitar el *4
SBB.initGrid()#dx=ddx)
# Define the temporal parameters, e.g. the total time (total) and the file
# output interval (file_dt), both in seconds

print 1*fcT,.005*fcT,.5e-6*fcT
SBB.initTemporal(total = 44, file_dt = .200, dt=2e-6)


# Using a 'dry' run, the sphere main program will display important parameters.
# sphere will end after displaying these values.
SBB.run(dry = True)
### RUNNING THE SIMULATION
# Start the simulation on the GPU from the sphere program
SBB.run()
SBB.writeVTKall()
