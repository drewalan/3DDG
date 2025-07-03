import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import get_test_data
import h5py as h5
import os
import copy as cp
import matplotlib.ticker as ticker

# PARAMETERS TO CHANGE
plotOverTime = False
plotOverSpace = True
plotBoolArray=[1,1,1,1,1,1,1,1]#use 0 or 1 and NO OTHER VALUES
svg=False
plotOverSpaceStride = 1 #plot only timesteps that are a multiple of this number
Dim1D = 0
DimOmit2D = 2  # which dimension is NOT sampled in 2D
order = 3
valsOverTime = [[],[],[],[],[]]#t, dRho, Ux, Uy, Uz
startTimestep = 0
endTimeStep = 100000#allowed to be an overestimate, but will take slightly longer to run
dt = 0.000001#for converting to ms on labels

folderName = "results"
plotFolderName = folderName+"/plots"
logname = "_log.txt"

if not os.path.exists(plotFolderName):
	os.makedirs(plotFolderName)

numElem = [0,0,0]
numBlocks = [0,0,0]
with open(folderName+"/"+logname) as f:
	for line in f:
		if "NUM_ELEMENTS = " in line:
			substr=line[line.find(" = ")+3:]
			values=substr.split()
			for i in [0,1,2]:
				numElem[i]=int(values[i])
		if "NUM_BLOCKS = " in line:
			substr=line[line.find(" = ")+3:]
			values=substr.split()
			for i in [0,1,2]:
				numBlocks[i]=int(values[i])
	print("Param file read, B=",numBlocks,", E=",numElem)
folderName +="/raw"

numProcs = numBlocks[0] * numBlocks[1] * numBlocks[2]
fixedPlotMax = False #set to false to make dynamic, to a nonzero float to set it to that height
fixedPlotMin = False #set to false to make dynamic, to a nonzero float to set it to that height
fixedPlotXMin = False #set to false to make dynamic, to a nonzero float to set it to that width
fixedPlotXMax = False #set to false to make dynamic, to a nonzero float to set it to that width

def find_nearest(array, value):
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	return idx

def unflatten(oldVals):
	numPts = [numElem[0] * (order+1),numElem[1] * (order+1),numElem[2] * (order+1)]
	vals = np.zeros((numPts[0],numPts[1],numPts[2]))
	for i in range(numPts[0]):
		for j in range(numPts[1]):
			for k in range(numPts[2]):
				idx = i * numPts[2] * numPts[1] + j * numPts[2] + k
				vals[i,j,k] = oldVals[idx]
	return vals


#compute node index offsets from rank
_dx=np.zeros((numProcs),dtype=int)
_dy=np.zeros((numProcs),dtype=int)
_dz=np.zeros((numProcs),dtype=int)
r=0
for i in range(numBlocks[0]):	
	for j in range(numBlocks[1]):	
		for k in range(numBlocks[2]):
			#order, not order+1, b/c redundant pts at elem boundary
			_dx[r]=int(i*numElem[0]*(order+1))
			_dy[r]=int(j*numElem[1]*(order+1))
			_dz[r]=int(k*numElem[2]*(order+1))
			r+=1

#setup empty arrays for read-in
numPts = [numBlocks[0]*numElem[0] * (order+1),numBlocks[1]*numElem[1] * (order+1),numBlocks[2]*numElem[2] * (order+1)]
xVals=np.zeros((numPts[0],numPts[1],numPts[2]))
yVals=np.zeros((numPts[0],numPts[1],numPts[2]))
zVals=np.zeros((numPts[0],numPts[1],numPts[2]))
ba=np.zeros((numPts[0],numPts[1],numPts[2]))
diff=np.zeros((numPts[0],numPts[1],numPts[2]))
rho0=np.zeros((numPts[0],numPts[1],numPts[2]))
c0=np.zeros((numPts[0],numPts[1],numPts[2]))
dRhoVals=np.zeros((numPts[0],numPts[1],numPts[2]))
VxVals=np.zeros((numPts[0],numPts[1],numPts[2]))
VyVals=np.zeros((numPts[0],numPts[1],numPts[2]))
VzVals=np.zeros((numPts[0],numPts[1],numPts[2]))

#Read in space and constants data from HDF5
for rank in range(numProcs):
	xmin=_dx[rank]
	xwidth=numElem[0]*(order+1)
	xmax=_dx[rank]+xwidth
	ymin=_dy[rank]
	ywidth=numElem[1]*(order+1)
	ymax=_dy[rank]+ywidth
	zmin=_dz[rank]
	zwidth=numElem[2]*(order+1)
	zmax=_dz[rank]+zwidth
	with h5.File(folderName+"/coords_"+str(rank)+"_0.h5", "r") as f:
		xVals[xmin:xmax,ymin:ymax,zmin:zmax]=f['coords'][0,0:xwidth,0:ywidth,0:zwidth]
		yVals[xmin:xmax,ymin:ymax,zmin:zmax]=f['coords'][1,0:xwidth,0:ywidth,0:zwidth]
		zVals[xmin:xmax,ymin:ymax,zmin:zmax]=f['coords'][2,0:xwidth,0:ywidth,0:zwidth]
	with h5.File(folderName+"/ba_"+str(rank)+"_0.h5", "r") as f:
		ba[xmin:xmax,ymin:ymax,zmin:zmax]=f['ba'][0:xwidth,0:ywidth,0:zwidth]
	with h5.File(folderName+"/diff_"+str(rank)+"_0.h5", "r") as f:
		diff[xmin:xmax,ymin:ymax,zmin:zmax]=f['diff'][0:xwidth,0:ywidth,0:zwidth]
	with h5.File(folderName+"/rho0_"+str(rank)+"_0.h5", "r") as f:
		rho0[xmin:xmax,ymin:ymax,zmin:zmax]=f['rho0'][0:xwidth,0:ywidth,0:zwidth]
	with h5.File(folderName+"/c0_"+str(rank)+"_0.h5", "r") as f:
		c0[xmin:xmax,ymin:ymax,zmin:zmax]=f['c0'][0:xwidth,0:ywidth,0:zwidth]

#Find array indices of selected point of interest (for 2d and 1d slices)
xInterest=find_nearest(xVals[:,0,0],0.)
yInterest=find_nearest(yVals[0,:,0],0.)
zInterest=find_nearest(zVals[0,0,:],0.)

#Read in all state variable data
for timestep in range(startTimestep, endTimeStep+1):
	skip=False
	try:
		for rank in range(numProcs):
			xmin=_dx[rank]
			xwidth=numElem[0]*(order+1)
			xmax=_dx[rank]+xwidth
			ymin=_dy[rank]
			ywidth=numElem[1]*(order+1)
			ymax=_dy[rank]+ywidth
			zmin=_dz[rank]
			zwidth=numElem[2]*(order+1)
			zmax=_dz[rank]+zwidth
			with h5.File(folderName+"/dRho_"+str(rank)+"_"+str(timestep)+".h5", "r") as f:
				dRhoVals[xmin:xmax,ymin:ymax,zmin:zmax]=f['rho'][0:xwidth,0:ywidth,0:zwidth]
			with h5.File(folderName+"/vel_"+str(rank)+"_"+str(timestep)+".h5", "r") as f:
				VxVals[xmin:xmax,ymin:ymax,zmin:zmax]=f['velocity'][0,0:xwidth,0:ywidth,0:zwidth]
				VyVals[xmin:xmax,ymin:ymax,zmin:zmax]=f['velocity'][1,0:xwidth,0:ywidth,0:zwidth]
				VzVals[xmin:xmax,ymin:ymax,zmin:zmax]=f['velocity'][2,0:xwidth,0:ywidth,0:zwidth]
	except:
		skip=True #read in for this timestep did not occur (expected for most steps if printFreq>1)

	if not skip:
		print("T=",timestep)
		print("Max DRho: ",dRhoVals.max())
		if plotOverTime:
			valsOverTime[0].append(timestep)
			valsOverTime[1].append(dRhoVals[xInterest][yInterest][zInterest])
			valsOverTime[2].append(VxVals[xInterest][yInterest][zInterest])
			valsOverTime[3].append(VyVals[xInterest][yInterest][zInterest])
			valsOverTime[4].append(VzVals[xInterest][yInterest][zInterest])

		if plotOverSpace and timestep % plotOverSpaceStride == 0:
			xCent=int(len(xVals[:,0,0])*.5)
			yCent=int(len(xVals[0,:,0])*.5)
			zCent=int(len(xVals[0,0,:])*.5)

			# --------------------------------------------------------------------------
			NUMPLOTS=np.sum(plotBoolArray)
			if NUMPLOTS==1:
				NUMPLOTROWS=1
				NUMPLOTCOLS=1
			elif NUMPLOTS==2:
				NUMPLOTROWS=1
				NUMPLOTCOLS=2
			elif NUMPLOTS==3:
				NUMPLOTROWS=1
				NUMPLOTCOLS=3
			elif NUMPLOTS==4:
				NUMPLOTROWS=2
				NUMPLOTCOLS=2
			elif NUMPLOTS==5:
				NUMPLOTROWS=2
				NUMPLOTCOLS=3
			elif NUMPLOTS==6:
				NUMPLOTROWS=2
				NUMPLOTCOLS=3
			elif NUMPLOTS==7:
				NUMPLOTROWS=2
				NUMPLOTCOLS=4
			elif NUMPLOTS==8:
				NUMPLOTROWS=2
				NUMPLOTCOLS=4

			#Lots of index magic here to prevent an enormous web of nested if statements
			#Essentially slice1D and slice2D function as a *range* of indices:
			#	e.g. instead of myArray[:,0,:] or myArray[:,:,0], both can be represented as myArray[slice2D]
			slice2D=[slice(None), slice(None), slice(None)]#'4d plot' must have one dimension flattened
			slice1D=[xCent, yCent, zCent]#'0d plot' must have one dimension expanded
			Dims2D=[0,1,2]
			Dims2D.remove(DimOmit2D)
			Dims2D = tuple(Dims2D)
			slice2D[DimOmit2D]=slice1D[DimOmit2D]#set one dimension to just be the center
			slice1D[Dim1D]=(slice(None))#set one dimension to be the whole span
			slice1D=tuple(slice1D)
			slice2D=tuple(slice2D)
			physicalDimensions=[xVals,yVals,zVals]#allows us to select a dimension with an integer variable as an index
			dimNames=["x", "y", "z"]#allows us to select a dimension with an integer variable as an index


			fig = plt.figure(figsize=(NUMPLOTCOLS*4,NUMPLOTROWS*3))
			axes=[0]*(NUMPLOTCOLS*NUMPLOTROWS)
			plotIdx=-1

			# --------------------------------------------------------------------------
			if plotBoolArray[0]:
				plotIdx+=1
				axes[plotIdx] = fig.add_subplot(NUMPLOTROWS, NUMPLOTCOLS, plotIdx+1, projection='3d')
				axes[plotIdx].set_title(r'$\delta$$\rho$')
				axes[plotIdx].set(xlabel=dimNames[Dims2D[0]], ylabel=dimNames[Dims2D[1]])

				surf = axes[plotIdx].plot_surface(physicalDimensions[Dims2D[0]][slice2D], physicalDimensions[Dims2D[1]][slice2D], dRhoVals[slice2D], rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)

				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_zlim(fixedPlotMin, fixedPlotMax)
				# fig.colorbar(surf, shrink=0.5, aspect=10)
			# --------------------------------------------------------------------------
			if plotBoolArray[1]:
				plotIdx+=1
				axes[plotIdx] = fig.add_subplot(NUMPLOTROWS, NUMPLOTCOLS, plotIdx+1, projection='3d')
				axes[plotIdx].set_title('$u_{x}$')
				axes[plotIdx].set(xlabel=dimNames[Dims2D[0]], ylabel=dimNames[Dims2D[1]])
				surf = axes[plotIdx].plot_surface(physicalDimensions[Dims2D[0]][slice2D], physicalDimensions[Dims2D[1]][slice2D], VxVals[slice2D], rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_zlim(fixedPlotMin, fixedPlotMax)
			# --------------------------------------------------------------------------
			if plotBoolArray[2]:
				plotIdx+=1
				axes[plotIdx] = fig.add_subplot(NUMPLOTROWS, NUMPLOTCOLS, plotIdx+1, projection='3d')
				axes[plotIdx].set_title('$u_{y}$')
				axes[plotIdx].set(xlabel=dimNames[Dims2D[0]], ylabel=dimNames[Dims2D[1]])
				surf = axes[plotIdx].plot_surface(physicalDimensions[Dims2D[0]][slice2D], physicalDimensions[Dims2D[1]][slice2D], VyVals[slice2D], rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_zlim(fixedPlotMin, fixedPlotMax)
			# --------------------------------------------------------------------------
			if plotBoolArray[3]:
				plotIdx+=1
				axes[plotIdx] = fig.add_subplot(NUMPLOTROWS, NUMPLOTCOLS, plotIdx+1, projection='3d')
				axes[plotIdx].set_title('$u_{z}$')
				axes[plotIdx].set(xlabel=dimNames[Dims2D[0]], ylabel=dimNames[Dims2D[1]])
				surf = axes[plotIdx].plot_surface(physicalDimensions[Dims2D[0]][slice2D], physicalDimensions[Dims2D[1]][slice2D], VzVals[slice2D], rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_zlim(fixedPlotMin, fixedPlotMax)
			# --------------------------------------------------------------------------
			if plotBoolArray[4]:
				plotIdx+=1
				axes[plotIdx] = fig.add_subplot(NUMPLOTROWS, NUMPLOTCOLS, plotIdx+1)
				#ax.set_title("DRho Slice in "+dimNames[Dim1D])
				axes[plotIdx].set(xlabel=dimNames[Dim1D], ylabel=r'$\delta$$\rho$')
				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_ylim(fixedPlotMin, fixedPlotMax)
				if fixedPlotXMax!=False and fixedPlotXMin!=False:
					axes[plotIdx].set_xlim(fixedPlotXMin, fixedPlotXMax)
				axes[plotIdx].yaxis.set_label_position("right")
				axes[plotIdx].yaxis.tick_right()
				axes[plotIdx].plot(physicalDimensions[Dim1D][slice1D], dRhoVals[slice1D],color='k')

				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_ylim(fixedPlotMin, fixedPlotMax)
			# --------------------------------------------------------------------------
			if plotBoolArray[5]:
				plotIdx+=1
				axes[plotIdx] = fig.add_subplot(NUMPLOTROWS, NUMPLOTCOLS, plotIdx+1)
				#ax.set_title("Ux Slice in "+dimNames[Dim1D])
				axes[plotIdx].set(xlabel=dimNames[Dim1D], ylabel='$u_{x}$')
				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_ylim(fixedPlotMin, fixedPlotMax)
				axes[plotIdx].yaxis.set_label_position("right")
				axes[plotIdx].yaxis.tick_right()
				axes[plotIdx].plot(physicalDimensions[Dim1D][slice1D], VxVals[slice1D],color='k')

				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_ylim(fixedPlotMin, fixedPlotMax)
			# --------------------------------------------------------------------------
			if plotBoolArray[6]:
				plotIdx+=1
				axes[plotIdx] = fig.add_subplot(NUMPLOTROWS, NUMPLOTCOLS, plotIdx+1)
				#ax.set_title("Uy Slice in "+dimNames[Dim1D])
				axes[plotIdx].set(xlabel=dimNames[Dim1D], ylabel='$u_{y}$')
				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_ylim(fixedPlotMin, fixedPlotMax)
				axes[plotIdx].yaxis.set_label_position("right")
				axes[plotIdx].yaxis.tick_right()
				axes[plotIdx].plot(physicalDimensions[Dim1D][slice1D], VyVals[slice1D],color='k')

				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_ylim(fixedPlotMin, fixedPlotMax)
			# --------------------------------------------------------------------------
			if plotBoolArray[7]:
				plotIdx+=1
				axes[plotIdx] = fig.add_subplot(NUMPLOTROWS, NUMPLOTCOLS, plotIdx+1)
				#ax.set_title("Uz Slice in "+dimNames[Dim1D])
				axes[plotIdx].set(xlabel=dimNames[Dim1D], ylabel='$u_{z}$')
				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_ylim(fixedPlotMin, fixedPlotMax)
				axes[plotIdx].yaxis.set_label_position("right")
				axes[plotIdx].yaxis.tick_right()
				axes[plotIdx].plot(physicalDimensions[Dim1D][slice1D], VzVals[slice1D],color='k')

				if fixedPlotMax!=False and fixedPlotMin!=False:
					axes[plotIdx].set_ylim(fixedPlotMin, fixedPlotMax)
			# --------------------------------------------------------------------------

			#plt.subplots_adjust(left=0.1,bottom=0.2,right=0.8,top=0.9,wspace=0.5,hspace=0.2)#good for mix of 2d & 2d
			plt.subplots_adjust(left=0.05,bottom=0.05,right=0.9,top=0.95,wspace=0.2,hspace=0.2)#good for 2d

			
			#labels must be set AFTER all 3d plots are drawn to override rotation. (only works on last 3d plot for unknown reason)
			labels=[r'   $\delta$$\rho$','   $u_{x}$','   $u_{y}$','   $u_{z}$']
			for idx in range(0,NUMPLOTCOLS*NUMPLOTROWS):
				try:
					axes[idx].zaxis.set_rotate_label(False)
					axes[idx].zaxis.set_tick_params(labelsize=8)
					axes[idx].set_zlabel(labels[idx],rotation=0)
					axes[idx].yaxis.set_rotate_label(False)
				except:
					pass
				axes[idx].xaxis.set_major_locator(ticker.MultipleLocator(.25*(xVals[-1][0][0]-xVals[0][0][0])))#control axes tick frequency
				axes[idx].yaxis.set_major_locator(ticker.MultipleLocator(.25*(yVals[0][-1][0]-yVals[0][0][0])))#control axes tick frequency
				axes[idx].xaxis.set_tick_params(labelsize=8)
				axes[idx].yaxis.set_tick_params(labelsize=8)

				ylabel=axes[idx].get_ylabel()
				axes[idx].set_ylabel(ylabel,rotation=0)
			

			plt.savefig(plotFolderName+"/_MasterPlot"+str(timestep))
			if svg:
				plt.savefig(plotFolderName+"/_MasterPlot"+str(timestep),format='svg')
			plt.clf()
			plt.close()

			old=dRhoVals #keep track to enable "change from previous timestep" mode


