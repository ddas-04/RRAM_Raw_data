
def trapz(x,y):
    n=len(y)
    dx=x[1]-x[0]
    index=1
    sum=0
    while index<=n-1:        
        fx=0.5*(y[index-1]+y[index])
        sum=sum+fx*dx
        index=index+1
    
    return sum


from fipy import *
import numpy as np
import scipy
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':15})




################ Set up simulation grid #######################
nx = 100e2;
dx = 1.0e-2;
L = nx * dx
midPoint = L / 2
mesh = Grid1D(nx=nx, dx=dx)
#########################################################
#......... Defining mean and S.D of the Gaussian distribution
mu=midPoint
sigma=L/8

################ Defining Phi variable #######################
phi = FaceVariable(name="solution variable", 
                   mesh=mesh,
                   value=1/L)
'''
valueLeft = 0.0
valueRight = 1.0

phi.constrain(valueLeft, mesh.facesLeft)
phi.constrain(valueRight, mesh.facesRight)
'''
############# Assiging the facecenters to Xgrid #################
Xgrid = mesh.faceCenters[0]

### save Xgrid in the form of numpy array
np.save('Xvalues.npy', Xgrid)
Xvalues=[]
Xvalues=np.load('Xvalues.npy')
#########################################################

################### Convection term ########################
scale_factor=1.0e6
const=1.0/((sigma**3)*numerix.sqrt(2*numerix.pi))
exp_power=-((Xgrid-mu)**2)/(2*sigma**2)

c_value=const*(Xgrid-mu)*numerix.exp(exp_power)
c_value=c_value*scale_factor

cCoeff = FaceVariable(mesh=mesh, value=[c_value,])


### Save c_value in numpy array
np.save('C_values.npy',c_value)
CValues=[]
CValues=np.load('C_values.npy')

######## Calculate Gaussian distribution for ploting purpose only #####
const_F=1.0/((sigma)*np.sqrt(2*numerix.pi))
exp_power=-((Xvalues-mu)**2)/(2*sigma**2)
Fx=1.0-const_F*np.exp(exp_power)

Fx_act=1-Fx

##########################################################
Fx=Fx*scale_factor
'''
############# log rho check for parabola #################
fig = plt.figure(figsize=(10,8))

ax = fig.add_subplot(2, 1, 1)
line, = ax.plot(Xvalues,(Fx_act))
plt.grid()

ax = fig.add_subplot(2, 1, 2)
line, = ax.plot(Xvalues,np.log((Fx_act)))
plt.grid()
plt.xlim([Xvalues[0],Xvalues[-1]])
plt.ylabel('Log')
plt.show()

exit()
'''
######### Plot of F(x) and its derivative ##################
fig = plt.figure(figsize=(12,10))

ax = fig.add_subplot(2, 1, 1)
line, = ax.plot(Xvalues,Fx)
plt.xlabel('Xgrid')
plt.ylabel(r"$(1-F(x))$")
plt.grid()
#plt.show()

ax = fig.add_subplot(2, 1, 2)
line, = ax.plot(Xvalues,CValues)
plt.xlabel('Xgrid')
plt.ylabel(r"$C=\frac{d}{dx}(1-F(x))$")
plt.grid()
plt.show()
fig.savefig('Distribution_and_its_derivative.png')
#exit()


### Diffusion term
Dval= 1.0e3
D = FaceVariable(mesh=mesh, value=Dval)
### create face variable
phi = CellVariable(name=r"$\rho$",
                   mesh=mesh,
                   value=1/L)
'''
phi.faceGrad.constrain([0.], mesh.facesLeft)
phi.faceGrad.constrain([0.], mesh.facesRight)
'''

eqX = TransientTerm() == (DiffusionTerm(coeff=D) + PowerLawConvectionTerm(coeff=cCoeff))

viewer = Viewer(vars=(phi),
                datamin=0., datamax=0.2)

timeStepDuration = 1.0e-6
steps = 300
## Loop to step through time
dexp=-15.0
limit=0.0
incr=0.05
number_of_steps=int((limit-dexp)/incr)+1
t_i=0

while t_i<number_of_steps:
    print("Loop Step Count = "+str(t_i))
    print("dexp = " + str(dexp))
    timeStepDuration=numerix.exp(dexp)
    print('timestep='+str(timeStepDuration))
    eqX.solve(var=phi, dt=timeStepDuration)
    print(max(phi))
    
    if max(phi)>=20:
	break
    print('**********************************************')
    if __name__ == '__main__':
        viewer.plot()
	#viewer.grid()
    dexp=dexp+0.05
    t_i=t_i+1


f=open("Main_result_data.txt","w+")

############ Save steady-state phi values ####################
np.save('phi_steady_state.npy',phi)
phi_steady_state=[]
phi_steady_state=np.load('phi_steady_state.npy')

#print('shape X = ' +str(np.shape(Xvalues)))
#print('shape phi = ' +str(np.shape(phi_staedy_state)))
total_probability=trapz(Xvalues,phi_steady_state)
f.write("Total Probability = %f\r\n" % (total_probability))
print('Total probability = ' + str(total_probability))
############### plot of log(rho) ##########################
xdata=Xvalues[0:len(Xvalues)-1]
print(type(xdata))
#ydata=np.log(phi_staedy_state)
ydata=phi_steady_state

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)
line, = ax.plot(xdata,ydata, linewidth=2)
ax.set_yscale('log')
plt.xlabel('Xgrid')
#plt.ylabel(r"$log(\rho)$")
plt.ylabel(r"$\rho$")

plt.grid()
plt.show()
############### Calculation of CDF ######################


N=len(phi_steady_state)
CDF=np.zeros(N)
i=2
median=[]
while i<=N:
	print('i = '+str(i))
	x=xdata[0:i]
	y=phi_steady_state[0:i]
	CDF[i-1]=trapz(x,y)
	print('CDF='+str(CDF[i-1]))
	print('---------------------------------------------')
	if CDF[i-1]>=0.5:
		median_val=x[-1]
		median=np.append(median,median_val)
	i=i+1
print('Median = '+str(median[0]))
f.write("Median = %f\r\n" % (median[0]))

fig = plt.figure(figsize=(12,10))

ax = fig.add_subplot(2, 1, 1)
line, = ax.plot(xdata,phi_steady_state, linewidth=2)
plt.xlabel('Xgrid')
plt.ylabel(r"$\rho$")
plt.ylim([0,0.2])
plt.grid()
#plt.show()

ax = fig.add_subplot(2, 1, 2)
line, = ax.plot(xdata,CDF, linewidth=2)
plt.grid()
plt.xlabel('Xgrid')
plt.ylabel('CDF')
fig.savefig('CDF.png')
plt.show()

f.close() 
if __name__ == '__main__':
    input("Transient drift-diffusion. Press <return> to proceed...")

