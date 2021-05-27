"""
CityGen - A Maya Python script used to generate a city.
Created as a 1st year programming project for 
Computer Animation Technical Arts at Bournemouth University.
Made by Peter Dravecky
"""

for name in dir():
    if not name.startswith('_'):
        try:
            del globals()[name]
        except:
            pass

import maya.cmds as mc
import maya.api.OpenMaya as om
import maya.mel as mel

import math
import random
from functools import partial
import os

projectDirectory = mc.workspace(q=True, rd=True)
if os.path.exists(projectDirectory+'scripts/mayaCityGen.py'):
    scriptPath = projectDirectory+"scripts"

def remainder(a,b):
    """Returns the remainder of a divided by b*  
    *Python uses different modulo calculation than C++ - this is the C++ one

    a - (float), b - (float)

    return - (float)
    """
    return a - int(a/b) * b

def ppdistance(a,b):
    """Returns the squared distance between the two 2D points.

    a - [float,float], b - [float,float]

    return - (float)
    """
    return (b[0]-a[0])**2+(b[1]-a[1])**2

def dotproduct(a,b):
    """Returns the dot product of two 2D vectors starting at origin.

    a - [float,float], b - [float,float]

    return - (float)
    """
    u = (a[0]*b[0]) + (a[1]*b[1])
    return u

def closestpoint(a1,a2,b):
    """Finds the closest point on line (a1,b2) to point b.

    a1 - [float,float], a2 - [float,float], b = [float,float]

    return - [float,float] - point on line (a1,b1)

    Source: https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
    """
    if a1 == a2:
        return [a1[0],a2[1]]

    ax = a2[0] - a1[0]
    ay = a2[1] - a1[1]

    norm = ax*ax + ay*ay

    u =  dotproduct( [b[0]-a1[0], b[1]-a1[1]] , [ax, ay]) / float(norm)
    
    if u > 1:
        u = 1
    elif u < 0:
        u = 0
    
    #print(u)

    x = a1[0] + u * ax
    y = a1[1] + u * ay
    #print(x)
    #print(y)
    return [x,y]

def lpdistance(a1,a2,b):
    """Returns the distance between a line (a1,a2) and b in 2D space.

    a1 - [float,float], a2 - [float,float], b - [float,float]

    return - (float) - distance
    """
    dist = ppdistance(closestpoint(a1,a2,b), b)

    return dist

def normalize(a):
    """Returns a normalized version of vector a.

    a - [float, float]

    return - [float, float]
    """
    length = (a[0]**2+a[1]**2)**0.5
    if(length> 0 ):
        return [a[0]/length, a[1]/length]
    else:
        print("normalize error")
        return a

def lcosangle(a,b):
    """returns cosine of the angle between 2 vectors starting at origin.

    a - [float,float], b - [float, float]

    return - (float) between -1.0 and 1.0
    """
    dotProduct = dotproduct(a, b)
    dotDiv = abs((a[0]*a[0] + a[1]*a[1])**0.5)*abs((b[0]*b[0] + b[1]*b[1])**0.5)

    if(dotDiv!=0):
        cosAngle = dotProduct/dotDiv
    else:
        cosAngle = 0

    if(cosAngle>1.0):
        cosAngle = 1.0
    elif(cosAngle<-1.0):
        cosAngle = -1.0
        
    return cosAngle

class Surface: 
    """Contains information and functions 
    connected to the plane the city is generated on.
    
    __init__ is used to pass in the UI class, prepareSurface() should be used to fully initialize it
    """
    name = ""
    ready = False
    subd = 60

    def __init__(self, menu):
        """menu - (UI Class)
        """
        self.ui = menu
        self.name = self.ui.surfaceName

    def addPaintAttr(self, attrName):
        """Creates a new paintable attribute and assigns it to the mesh connected to the Surface class.

            attrName - (string) - name of the attribute
        """
        mc.select(self.name, hi=True)
        mc.select(self.name, d=True)
        if mc.attributeQuery(attrName, node = self.shape, exists = True) == False:
            mc.addAttr(self.shape, ln = attrName, nn = attrName, dt = 'doubleArray')
        mc.makePaintable('mesh', attrName, at = 'doubleArray')
        mc.select(cl=True)
        mc.setToolTo('selectSuperContext')

    def removePaintAttr(self, attrName):
        """Removes given paintable attribute from the mesh connected to this Surface class.

            attrName - (string) - name of the attribute
        """
        mc.select(self.name, hi=True)
        mc.select(self.name, d=True)
        if mc.attributeQuery(attrName, node = self.shape, exists = True) == True:
            mc.deleteAttr(self.shape, at=attrName)
        mc.select(cl=True)
        mc.setToolTo('selectSuperContext')

    def getAttrAtPos(self, pos, attrArray):
        """Returns the attribute value from the mc.getAttr() matrix at given position relative to the surface.

            pos - [float, float], attrArray - [float, float...]

        return - (float) - attribute value at given pos
        """
        if(attrArray == None):
            return 0

        self.updateBounds()
        subd = self.subd + 1

        width = self.maxx - self.minx
        Xoffset = pos[0] - self.minx
        Xpos = int((Xoffset*subd)/width)


        height = self.maxz - self.minz
        Yoffset = self.maxz - pos[1]
        Ypos = int((Yoffset*subd)/height)-1
        index = Ypos*subd + Xpos

        try:
            attrIndex = attrArray[index]
        except:
            try:
                attrIndex = attrArray[index-1]
            except:
                attrIndex = 0.0

        return attrIndex

    def updateBounds(self):
        """Re-checks the world bounding box of the surface objects and updates the bound properties.
        
            bound properties - self.(bounds, minx, miny, minz, maxx, maxy, maxz)
        """
        self.bounds = mc.exactWorldBoundingBox(self.name)
        bounds = self.bounds
        #print("updating bounds: " + str(bounds))

        self.minx = bounds[0]
        self.miny = bounds[1]
        self.minz = bounds[2]

        self.maxx = bounds[3]
        self.maxy = bounds[4]
        self.maxz = bounds[5]

    def isInBounds(self, pos):
        """Checks if the objects 2D x and z coordinates lay within the bounds of the surface.

            pos - [float, float]

        return - (bool)
        """
        return (pos[0]<=self.maxx and pos[0]>=self.minx and pos[1]<=self.maxz and pos[1]>=self.minz)

    def isOnBounds(self, pos):
        """Checks if the objects 2D x and z coordinates lay on the bounds of the surface.

            pos - [float, float]

        return - (bool)
        """
        return (pos[0]==self.minx or pos[0]==self.maxx) or (pos[1]==self.minz or pos[1]==self.maxz)

    def findBoundOverlap(self, pos, vector):
        """Finds the position of the point on vector where vector intersects the bounds of the surface.

            pos - [float, float] - Position within the bounds

            vector - [float, float] - Vector which overlaps with the bounds if added to the pos

        return - [float, float] - The position on bounds
        """
        newPos = [0,0]
        #print(pos)

        outPos = [0,0]
        outPos[0] = pos[0] + vector[0]
        outPos[1] = pos[1] + vector[1]
        
        x = False
        z = False
        ratio = 0

        if not self.isInBounds(outPos):

            if((outPos[0]<self.minx or outPos[0]>self.maxx) and vector[0]!=0):
                x = True
                if(outPos[0]<self.minx):
                    xbounds = self.minx
                else:
                    xbounds = self.maxx
                xratio = (xbounds-pos[0])/vector[0]

            if((outPos[1]<self.minz or outPos[1]>self.maxz) and vector[1]!=0):
                z = True
                if(outPos[1]<self.minz):
                    zbounds = self.minz
                else:
                    zbounds = self.maxz
                zratio = (zbounds-pos[1])/vector[1]

            if(x and z):
                if(xratio<zratio):
                    ratio = xratio
                else:
                    ratio = zratio
            elif(x):
                ratio = xratio
            elif(z):
                ratio = zratio
            else:
                ratio = 1

            #print(ratio)

            newPos[0] = pos[0] + vector[0]*ratio
            newPos[1] = pos[1] + vector[1]*ratio
            #print(newPos)

        return newPos

    def prepareSurface(self):
        """Assigns the surface object to the Surface class. Should be executed at the beginning of the program.
        """
        ui = self.ui
        newSurface = self.name
        if not(mc.objExists(self.name)):
            newSurface = mc.polyPlane(n = self.name, w = self.subd, h = self.subd, sx = self.subd, sy = self.subd)
        mc.select(newSurface, hi=True)

        self.shape = mc.ls(sl=True, l=True)[1]

        _curSel = om.MGlobal.getActiveSelectionList()
        self.mesh = om.MFnMesh(_curSel.getDagPath(0))

        self.addPaintAttr(ui.streetAttrName)

        self.updateBounds()
        
class TensorField:
    """Contains information and functions
    connected to the Tensor Field, a field 
    made up of vectors whose rotation is influenced by the different types of Tensor Points.

    __init__ passes in the Surface and UI classes, no other initialization needed
    """
    pointArray = []
    pointNameArray = []
    dirArray = []
    radArray = []

    locatorArray = []

    def __init__(self, surface):
        """surface - (Surface Class) that the TensorField should work with
        """
        self.surface = surface
        self.ui = self.surface.ui

    def getTensor(self, pos):
        """Returns the Tensor information at the current position

            pos - [float, float]
        
        return - (TensorField.Tensor Class)
        """
        return self.Tensor(self, pos, True)

    def getDirTensor(self, _pos):
        """Returns the Tensor information at the current position, when the Tensor is only influence by directional Tensor Points

            pos - [float, float]
        
        return - (TensorField.Tensor Class)
        """
        return self.Tensor(self, _pos, False)

    class Tensor:
        """Contains information and functions
        connected to one specific point on the Tensor Field

        Creating a new tensor should be made through TensorField's getTensor() and getDirTensor()
        """
        pos = [0,0]
        major = 0
        minor = 90
        isRad = False   
        radPoints = []
        useRad = True

        def __init__(self, tensorField, pos, useRad):
            """ tensorField - (TensorField Class)

                pos - [float, float]

                useRad - (bool) - False if the tensor should ignore radial Tensor Points
            """
            self.pos = pos
            self.tensorField = tensorField
            self.useRad = useRad
            self.setRotation()
        
        def setRotation(self):
            """ Sets the rotation of the Tensor's vectors according to the Tensor Points that influence it.

                tensorField - (TensorField Class)

                pos - [float, float]
                
                useRad - (bool) - False if the tensor should ignore radial Tensor Points
            """
            self.radPoints = []
            tempRadArray = []
            tempAngle = 0
            distAngle = 0

            closestDist = -1

            #Get directional Tensor Points
            for eachPoint in self.tensorField.dirArray:
                distance = ppdistance(self.pos, eachPoint.pos)
                if distance<closestDist or closestDist == -1:
                    closestDist = distance
                    distAngle = eachPoint.angle

            if(self.useRad):
                #Get radial Tensor Points
                for eachPoint in self.tensorField.radArray:
                    tensorVector = self.pos
                    pointVector = eachPoint.pos
                    distance = ppdistance(self.pos, eachPoint.pos)
                    size = (eachPoint.radius)**2
                    blendRatio = 1
                    if(distance<=size):
                        distVector = [tensorVector[0]-pointVector[0] , tensorVector[1]-pointVector[1]]
                        
                        if(distVector[0]<0):
                            distVector[1] *= -1
                            distVector[0] *= -1
                        
                        tempAngle = mc.angleBetween(euler=True, v1=(1,0,0), v2=(distVector[0], 0, distVector[1]))[1]
                        if(distance/size > 1 - eachPoint.decay):
                            self.isRad == False
                            blendRatio = (1-distance/size)/eachPoint.decay
                            tempAngle *= blendRatio
                            tempAngle += (1-blendRatio)*distAngle
                        else:
                            self.isRad = True
                            self.radPoints.append(eachPoint)

                        tempRadArray.append((eachPoint, blendRatio, tempAngle))

            self.radArray = tempRadArray
            if len(tempRadArray)>1:
                for x in range(len(tempRadArray)-1):
                    tempAngle += tempRadArray[x][2]
                tempAngle /= len(tempRadArray)
            elif len(tempRadArray)==0:
                tempAngle = distAngle

            #Clamps the major angle between 45 and -45
            tempAngle = remainder(tempAngle, 90)
            if(tempAngle<-45): tempAngle +=90
            elif(tempAngle>45): tempAngle -=90
            
            self.major = tempAngle          #horizontal road angle - angles between -45 and 45
            self.minor = tempAngle + 90     #vertical road angle - angles between 45 and 135
        
    def addPoint(self, _type, *pArgs):
        """Creates a TensorPoint class object and adds it to the list of all TensorPoints

        return - (TensorField.TensorPoint Class)
        """
        if not (mc.objExists(self.ui.tensorPointGroupName)):
            mc.group(em = True, name = self.ui.tensorPointGroupName)
        newPoint = self.TensorPoint(self, _type)
        self.pointArray.append(newPoint)
        self.pointNameArray.append(newPoint.object)
        if(newPoint.type == "dir"):
            self.dirArray.append(newPoint)
        else: self.radArray.append(newPoint)

        return newPoint

    def updatePoints(self):
        """Checks the tensor group and reimports all the tensor points.
        """
        if(mc.objExists(self.ui.tensorPointGroupName)):
            mc.select(self.ui.tensorPointGroupName, hi=True)
            objects = mc.ls(sl=True, tr=True)
            objectsLong = mc.ls(sl=True, tr=True, l=True)
            x=0
            while x < len(objects):
                if(objects[x][0] == "d"):
                    newPoint = self.addPoint('dir')
                    mc.matchTransform(newPoint.object, objectsLong[x])
                elif(objects[x][0] =="r"):
                    newPoint = self.addPoint('rad')
                    mc.matchTransform(newPoint.object, objectsLong[x])
                    mc.matchTransform(newPoint.decayTransform.fullPathName(), objectsLong[x+1])
                    x+=1
                x+=1
            mc.delete(objectsLong[1::])
        deletedPoints = []
        for eachPoint in self.pointArray:
            if not mc.objExists(eachPoint.object):
                if eachPoint.type=="dir":
                    self.dirArray.remove(eachPoint)
                else: self.radArray.remove(eachPoint)
                deletedPoints.append(eachPoint)
                self.pointNameArray.remove(eachPoint.object)
            else:    
                eachPoint.update()
        for eachPoint in deletedPoints:
            self.pointArray.remove(eachPoint)

    class TensorPoint:
        """Contains the information and functions connected to
        a specific Tensor Point.

        Creating a new TensorPoint object should only be done through TensorField's addPoint() function.
        """
        type = "rad"
        angle = 0
        radius = 10
        decay = 0.2
        pos = [0,0]

        def __init__(self, tensorField, type):
            """ tensorField - (TensorField Class)

                type - (string) - type of the Tensor Point - either "dir" or "rad"
            """
            self.tensorField = tensorField
            if type not in ['rad','dir']:
                print("Warning: Unknown tensor point type.")
            self.type = type
            self.visualizePoint()

        def visualizePoint(self):
            """Visualizes the TensorPoint information in the workspace as a movable control.
            """
            if(self.type == "rad"):
                print("created a rad point")
                newPoint = mc.circle(n="radPoint", nr=(0,1,0), c=(0,0,0), r=1)

                _curSel = om.MGlobal.getActiveSelectionList()
                mainTransform = om.MFnTransform().setObject(_curSel.getDependNode(0))

                decayTransform = mainTransform.duplicate()
                mainTransform.addChild(decayTransform)
                decayTransform = om.MFnTransform().setObject(decayTransform)
                decayName = decayTransform.fullPathName()
                mc.select(decayName)
                mc.rename("fallofPoint")

                decayTransform.scaleBy((0.8,0.8,0.8))
                mainTransform.scaleBy((10,10,10))
                
                self.mainTransform = mainTransform
                self.decayTransform = decayTransform

                mc.setAttr(decayTransform.fullPathName()+".translate",l=True)

                self.decayTransform.setLimit(1, 1)
                self.decayTransform.setLimit(5, 1)
                self.object = newPoint[0]

            if(self.type == "dir"):
                print("created a dir point")
                newPoint = mc.curve(n="dirPoint", d=1, p=[(0,0,4),(-2,0,1.5),(-1,0,1.5),(-1,0,-1.5),(-2,0,-1.5),(0,0,-4),(2,0,-1.5),(1,0,-1.5),(1,0,1.5),(2,0,1.5),(0,0,4)])
                _curSel = om.MGlobal.getActiveSelectionList()
                self.mainTransform = om.MFnTransform().setObject(_curSel.getDependNode(0))
                for limitType, amount in zip([0,1,2,3,4,5], [1,1,1,1,1,1]):
                    self.mainTransform.setLimit(limitType, amount)
                self.object = newPoint

            for limitType, amount in zip([20,21,12,13,16,17], [0,0,0,0,0,0]):
                self.mainTransform.setLimit(limitType, amount)
            
            mc.select(self.object)
            mc.parent(self.object, self.tensorField.ui.tensorPointGroupName, a=True)
            self.object = mc.ls(sl=True)[0]
            print(self.object)

        def update(self):
            """Updates the TensorPoint properties according to the transforms of the control connected to the TensorPoint
            """
            self.pos = [self.mainTransform.translation(1).x, self.mainTransform.translation(1).z]
            self.transform = (self.pos[0], 0, self.pos[1])
            if self.type == "rad":
                self.streetSizes = []
                self.radius = self.mainTransform.scale()[0]
                self.decay = 1 - self.decayTransform.scale()[0]
            self.angle = math.degrees(self.mainTransform.rotation(1).y)

    def visualize(self):
        """Visualizes the Tensor Field as evenly spaced group of Locators that rotate according to the tensor vectors.
        """
        if(mc.objExists(self.surface.ui.tensorGroupName)):
            mc.delete(self.surface.ui.tensorGroupName)

        self.locatorArray = []
        locatorNameArray = []
        _progress = self.ui.Progress('Visualizing tensor field...')

        mainLocator= mc.spaceLocator()
        mc.select(mainLocator, hi=True)
        mc.select(mainLocator, d=True)
        _curSel = om.MGlobal.getActiveSelectionList()
        locatorObject = _curSel.getDependNode(0)

        xRange = range(int(self.surface.minx), int(self.surface.maxx))[2::3]
        yRange = range(int(self.surface.minz), int(self.surface.maxz))[2::3]

        self.visBounds = self.surface.bounds[:]
        for x in xRange:
            for y in yRange:
                #print(str(x) + ", " + str(y))
                tempObject = om.MFnTransform().create()
                tempTransform = om.MFnTransform().setObject(tempObject)
                tempTransform.addChild(locatorObject, index=0, keepExistingParents=True)
                tempTransform.translateBy(om.MVector(x, 0,y), 1)
                locatorNameArray.append(tempTransform.fullPathName())
                self.locatorArray.append(tempTransform)
                _progress.add((len(xRange)*len(yRange))/100.0)
                

        mc.delete(mainLocator)
        mc.select(locatorNameArray)
        if not locatorNameArray == []:
            mc.group(n=self.surface.ui.tensorGroupName)

        _progress.finish()

    def updatePreVis(self, *pArgs):
        """Updates the rotation of locators made in the visualize() function and hides them in places where the obstructionMap is painted.

            *pArgs - Maya UI command arguments
        """
        self.surface.updateBounds()
        if self.visBounds != self.surface.bounds:
            mc.select(self.surface.ui.tensorGroupName)
            mc.delete()
            self.visualize()

        self.updatePoints()
        obstructionMap = mc.getAttr(self.surface.shape + '.' + self.ui.streetAttrName)
        for eachLocator in self.locatorArray:
            lPos = [eachLocator.translation(1).x, eachLocator.translation(1).z]

            tempTensor = self.getTensor(lPos)
            eachLocator.setRotation(om.MEulerRotation(0,math.radians(tempTensor.major),0),1)
            obstr = self.surface.getAttrAtPos(lPos, obstructionMap)
            if(obstr>0.5):
                mc.hide(eachLocator.fullPathName())
            else:
                mc.showHidden(eachLocator.fullPathName())

    def showPreVis(self, *pArgs):
        """Shows the locators after they've been hidden.
        """
        if(self.locatorArray == []):
            self.visualize()
            mc.button(self.ui.updateTensor, e=True, en = True, vis = True)
            self.updatePreVis()
        else:
            self.updatePreVis()
            mc.showHidden(self.surface.ui.tensorGroupName)

class StreetMap:
    """Contains properties and functions connected to
    street generation.

    __init__ passes in the other working classes, createMap() executes the main function of this class
    """
    streetArray = []
    radArray = []
    branchArray = []

    dstep = 1.0     #Step between creating a new point on the road

    def __init__(self, tensorField):
        """tensorField - (TensorField Class) you want the StreetMap to get tensor information from
        """
        self.tensorField = tensorField
        self.ui = self.tensorField.surface.ui
        self.basedist = 5
        self.updateDist(self.basedist)

    def updateDist(self, num):
        """Updates the properties avgdist, maxdist, mindist according to the num value

            num - (float)
        """
        self.avgdist = num
        self.maxdist = self.avgdist*2
        self.mindist = self.avgdist/2

    class Street:
        """Contains properties and functions connected to the street currently being generated.
        __init__ assigns the type and passes in the StreetMap class, getRadStreet() and getDirStreet() are used for singular street generation
        """
        def __init__(self, sm, startType):
            """ sm - (StreetMap Class)
                startType - (string) - type of the street being generated, either "minor", "major" or "rad"
            """
            self.type = startType
            self.streetMap = sm
            self.streetList = []
            self.pstreet = []
            self.bounds = []
            self.endEarly = False

        def getNoise(self):
            """Gets the noise amount from the menu and makes it either negative or positive.

            return - (float) - new noise amount
            """
            noise = 0
            noise_amount = self.streetMap.noise_amount
            if(noise_amount > 0):
                noise = int(random.choice([-1,1]))*random.random()*noise_amount
            return noise

        def isPointObstructed(self, curPos):
            """Checks if the value of the obstructionMap is higher than 0.5.

                curPos - [float,float]

            return - (bool)
            """
            sm = self.streetMap
            obstr = sm.tensorField.surface.getAttrAtPos(curPos, sm.obstructionMap)
            if(obstr > 0.5):
                return True
            else:
                return False

        def getClosestPointInfo(self, street, curPos):
            """Gets the information about the closest point on the street to curPos.

                street - ([float,float, float], [float,float, float]...)

                curPos - [float, float]

            Return:

                [float,float] - coordinates of the closest point

                (int) - index of the end point of a street segment that the closestPoint lies on.

                (float) - distance from curPos to the closestPoint

                (float) - cosine of the angle between the current line segment and the line between the curPos and closestPoint
            
            """
            #Finds closest point to bounds
            dist = 100000
            closestPoint = [0,0]
            closestPoint_index = 20
            cosAngle = -2

            for p in range(len(street)-1):
                p1 = [street[p][0], street[p][2]]
                p2 = [street[p+1][0], street[p+1][2]]

                tempClosestPoint = closestpoint(p1,p2, curPos)
                tempDist = ppdistance(tempClosestPoint, curPos)
                tempCosAngle = lcosangle([p2[0]-p1[0], p2[1]-p1[1]], [tempClosestPoint[0]-curPos[0], tempClosestPoint[1] - curPos[1]])

                if(tempDist <= dist):
                    closestPoint_index = p+1
                    dist = tempDist
                    closestPoint = tempClosestPoint
                    cosAngle = tempCosAngle

            return closestPoint, closestPoint_index, dist, cosAngle

        def findEarlyEndPoint(self, street, cp, cpi, pos, vector):
            """Gets the point on the nearby streets that will be the end point for the currently generated street.

                street - ([float,float, float], [float,float, float]...) - nearby street

                cp - [float, float] - closest point on street to pos

                cpi - (int - ) - index of the end point of a street segment that the cp lies on

                pos - [float, float] - current position

                vector - [float, float] - current street direction

            return - [float, float] - end point
            """
            p = [street[cpi][0]-cp[0], street[cpi][2]-cp[1]]
            cosAngle = lcosangle(vector,p)
            if(cosAngle>0.99 or cosAngle<=0):
                return cp

            pdist = ppdistance([0,0], p)**0.5
            cdist = ppdistance(pos, cp)**0.5
            tanAngle = math.tan(math.acos(cosAngle))

            newDist = cdist/tanAngle
            if newDist<pdist:
                np = [p[0]/pdist, p[1]/pdist]
                p = [np[0]*newDist, np[1]*newDist]

            endpoint = [cp[0]+p[0], cp[1]+p[1]]

            return endpoint

        def findBestConnectPoint(self, street, cp, cpi, pos, vector):
            """Calculates and returns a connect point to the nearby street, for which 
            the angles between the line from the connect point to the current point and the two streets 
            are more even than the angles, if the line was between the closest point and the current point.

                street - ([float,float, float], [float,float, float]...) - nearby street

                cp - [float, float] - closest point on street to pos

                cpi - (int - ) - index of the end point of a street segment that the cp lies on

                pos - [float, float] - current position

                vector - [float, float] - current street direction

            return - [float, float] - connect point
            """
            p = [street[cpi][0]-cp[0], street[cpi][2]-cp[1]]

            cosAngle = lcosangle(vector,p)
            print(cosAngle)

            angle = (math.pi - math.acos(cosAngle)) / 2

            pdist = ppdistance([0,0], p)**0.5
            cdist = ppdistance(pos, cp)**0.5
            tanAngle = math.tan(angle)

            newDist = cdist/tanAngle
            if newDist<pdist:
                np = [p[0]/pdist, p[1]/pdist]
                p = [np[0]*newDist, np[1]*newDist]

            endpoint = [cp[0]+p[0], cp[1]+p[1]]

            return endpoint

        def finishStreetPart(self, pos, streetPoints):
            """If executed, adds pos to streetPoints, appends streetPoints to the
            classes' steetList and clears the streetPoints list

                pos - [float, float] - position to finish the street at

                streetList - ([float,float, float], [float,float, float]...) - list of previous points

            return - []
            """
            streetPoints.append((pos[0],0,pos[1]))
            if(len(streetPoints)>1):
                self.streetList.append(streetPoints[:])
            return []

        def getSmoothVector(self, startTensor, curVector, snum):
            """This function smooths the forward vector by checking the tensors in front
            of it snum times.

                startTensor - (TensorField.Tensor Class) - tensor at current position

                curVector - [float, float] -  current forward vector

                snum - (int) - number of smoothing iterations (0 == no smoothing)

            return - [float, float] - smoothed vector
            """
            sm = self.streetMap
            tf = sm.tensorField

            curPos = startTensor.pos
            curTensor = tf.getTensor(curPos)
            angle = getattr(curTensor, self.type) + self.getNoise()

            for i in range(snum):

                curVector = [math.cos(math.radians(angle))*sm.dstep, -math.sin(math.radians(angle))*sm.dstep]

                newPos = [curPos[0]+curVector[0], curPos[1]+curVector[1]]
            
                newTensor = tf.getTensor(newPos)
                angle = (angle + getattr(newTensor, self.type) + self.getNoise())/2.0

                curPos = newPos
            
            curVector = [math.cos(math.radians(angle))*sm.dstep, -math.sin(math.radians(angle))*sm.dstep]
            return curVector

        def getSafeRadStreet(self, tempTensor, curPos, curVector, streetPoints, prevObstructed):
            """This function makes the directional Streets behave properly when encountering radial Tensor Points.

                tempTensor - (TensorField.Tensor Class) - tensor at current position

                curPos - [float, float] - current position

                curVector - [float, float] - current vector

                streetPoints - ([float,float, float], [float,float, float]...) - previous points on the street

                prevObstructed - (bool) - value checking if the point before current position was obstructed

            Return:

                (TensorField.Tensor Class) - tensor at finished position

                [float, float] - finished position

                [float, float] - previous position

                ([float,float, float], [float,float, float]...) - new street points
            """
            sm = self.streetMap
            tf = sm.tensorField
            newTensor = None
            new_pstreet = self.bounds
            print("getsaferad?")
            radPos = [0,0]
            radCenter = [0,0,0]
            endEarly = False

            for eachRad in tempTensor.radPoints:

                if(eachRad not in self.usedRads):
                    self.usedRads.append(eachRad)
                    rsSize = eachRad.streetSizes[-1]
                    radCenter = eachRad.transform

                    vec = [curPos[0]-radCenter[0], curPos[1]-radCenter[2]]
                    cosAngle = lcosangle(vec, curVector)

                    if(cosAngle<-0.9):

                        print("***"+str(cosAngle))

                        isObstructed = self.isPointObstructed(curPos)

                        if(isObstructed == False):
                            streetPoints.append((curPos[0],0,curPos[1]))

                        elif(isObstructed == True and prevObstructed == False):
                            prevPos = curPos
                            prevPos[0] = curPos[0]-curVector[0]
                            prevPos[1] = curPos[1]-curVector[1]
                            streetPoints = self.finishStreetPart(prevPos, streetPoints)

                        prevObstructed = isObstructed

                        new_pstreet.append((curPos[0],0,curPos[1])) 

                        normvec = normalize(vec)
                        
                        radPos1 = [radCenter[0]+normvec[0]*rsSize, radCenter[2]+normvec[1]*rsSize]

                        if not (tf.surface.isInBounds(radPos1)):
                            radPos1 = tf.surface.findBoundOverlap(curPos, [radPos1[0] - curPos[0], radPos1[1] - curPos[1]])

                        if not (self.isPointObstructed(radPos1)):
                            streetPoints.append((radPos1[0], 0, radPos1[1]))

                        new_pstreet.append((radPos1[0], 0, radPos1[1]))

                        streetPoints = self.finishStreetPart(radPos1, streetPoints)

                        if(self.type == "major"):
                            vec[0] *= -1
                            normvec[0] *= -1
                        else:
                            vec[1] *= -1
                            normvec[1] *= -1

                        radPos2 = [radCenter[0]+normvec[0]*rsSize, radCenter[2]+normvec[1]*rsSize]

                        if(tf.surface.isInBounds(radPos2) == False or endEarly == True):
                            newTensor = tf.getTensor(radPos1)
                            return newTensor, radPos1, curPos, streetPoints

                        if not (self.isPointObstructed(radPos2)):
                            streetPoints.append((radPos2[0], 0, radPos2[1]))

                        new_pstreet.append((radPos2[0], 0, radPos2[1]))

                        newPos = [radCenter[0]+vec[0], radCenter[2]+vec[1]]
                        curPos = newPos

                        if (tf.surface.isInBounds(curPos) == False):
                            curPos = tf.surface.findBoundOverlap(radPos2, [ curPos[0] - radPos2[0] , curPos[1] - radPos2[1] ])                      
                        
                        newTensor = tf.getTensor(curPos)
                    else:
                        print("not going in")
                        newTensor = tf.getDirTensor(curPos)
                else:
                    print("nah")
                    newTensor = tf.getDirTensor(curPos)

                return newTensor, curPos, [radCenter[0], radCenter[2]], streetPoints
            return newTensor, curPos, [curPos[0]-curVector[0], curPos[1]-curVector[1]], streetPoints

        def createDirStreet(self, startPos, pstreet, r):
            """Generates a directional street.

                startPos - [float, float] - starting position

                pstreet - ([float,float, float], [float,float, float]...) - previous street bounds

                r - (int) - recursion number
            """
            sm = self.streetMap
            tf = sm.tensorField
            r +=1

            if(r>20 or self.type not in ["minor","major"]):
                print("WARNING: Maximum recursion reached or invalid street type... r="+str(r))
                return

            isObstructed = False
            prevObstructed = False

            streetPoints =  []
            streetList = self.streetList

            curPos = startPos[:]
            curType = self.type
            curVector = [0,0]
            smooth_amount = sm.smooth_amount

            self.pstreet = pstreet
            new_pstreet = self.bounds
            
            self.usedRads = []

            prevPos = startPos[:]
            prevVector = [0,0]
            prevDist = -1

            endEarly = self.endEarly
            
            while(tf.surface.isInBounds(curPos) and (endEarly == False)):

                tempTensor = tf.getTensor(curPos) 

                if(tempTensor.isRad==True):
                    print("isRad")
                    tempTensor, curPos, prevPos, streetPoints = self.getSafeRadStreet(tempTensor, curPos, curVector, streetPoints, prevObstructed)

                isObstructed = self.isPointObstructed(curPos)

                if(pstreet != []):
                    closestPoint, closestPoint_index, dist, cosAngle = self.getClosestPointInfo(pstreet, curPos)

                    if prevDist == -1:
                        prevDist = dist
                    
                    if dist != prevDist and closestPoint != []: 

                        if(dist>sm.maxdist**2 and tempTensor.isRad == False and (cosAngle<0.7 and cosAngle>-0.7)):
                            closestPoint = self.findBestConnectPoint(pstreet, closestPoint, closestPoint_index, curPos, curVector)
                            branchPos = [(curPos[0]+closestPoint[0])/2, (curPos[1]+closestPoint[1])/2]
                            branchTensor = tf.getTensor(branchPos)

                            if(branchTensor.isRad == False):
                                old_pstreet = pstreet[:]
                                branchStreet = sm.Street(sm , curType)
                                branchStreet.createDirStreet(branchPos, pstreet, r)
                                pstreet = branchStreet.pstreet
                                
                                if(pstreet!= old_pstreet):
                                    sm.branchArray.append(branchStreet)

                        if(dist<sm.mindist**2 and closestPoint != pstreet[closestPoint_index]):
                            p1 = [pstreet[closestPoint_index-1][0],pstreet[closestPoint_index-1][2]]
                            p2 = [pstreet[closestPoint_index][0],pstreet[closestPoint_index][2]]
                            if not ((p1[1]==p2[1] and (p1[1]==tf.surface.minz)) or (p1[0]==p2[0] and (p1[0]==tf.surface.minx))):
                                endEarly = True

                curVector = self.getSmoothVector(tempTensor, curVector, smooth_amount)

                if(isObstructed == True and prevObstructed == False):
                    streetPoints = self.finishStreetPart(prevPos, streetPoints)

                if((curVector != prevVector or endEarly==True or (isObstructed == False and prevObstructed==True))):
                    prevVector = curVector
                    if(isObstructed == False):
                        streetPoints.append((curPos[0], 0, curPos[1]))
                    new_pstreet.append((curPos[0], 0, curPos[1]))

                prevObstructed = isObstructed

                prevPos = curPos
                curPos[0] += curVector[0]
                curPos[1] += curVector[1]

            prevPos[0] = curPos[0] - curVector[0]
            prevPos[1] = curPos[1] - curVector[1]


            if(endEarly == False):
                prevPos = tf.surface.findBoundOverlap(prevPos, curVector)

                streetPoints.append((prevPos[0], 0, prevPos[1]))
                new_pstreet.append((prevPos[0], 0, prevPos[1]))

            elif(endEarly == True and tempTensor.isRad == False and sm.do_branch == True):
                prevPos = self.findEarlyEndPoint(pstreet, closestPoint, closestPoint_index, prevPos, prevVector)
                streetPoints.append((prevPos[0], 0, prevPos[1]))
                new_pstreet.append((prevPos[0], 0, prevPos[1]))

            if(pstreet!=[]):
                i = closestPoint_index
                while(i < len(pstreet)):
                    new_pstreet.append(pstreet[i])
                    i += 1
            
            if(len(new_pstreet) <= 1 or (len(new_pstreet)==2 and new_pstreet[0]==new_pstreet[1])):
                self.pstreet = pstreet
                return

            if(len(streetPoints)>1):
                streetList.append(streetPoints[:])
            x=0
            for eachStreet in streetList:
                if(len(eachStreet)>2 or (len(eachStreet)==2 and eachStreet[0]!=eachStreet[1])):
                    newStreet = mc.curve(n=str(r)+"_"+str(x)+"roadTry", d=1, p=eachStreet)
                    x+=1
                    sm.streetArray.append(newStreet)
                    print("new Road created: " + str(newStreet) + str(eachStreet))
            mc.refresh()
            
            self.pstreet =  new_pstreet

        def createRadStreet(self, pos, radius):
            """Generates a radial (circular) street.

                pos - [float, float] - starting position

                radius - (float) - radius of the street
            """
            sm = self.streetMap
            tf = sm.tensorField
            streetPoints = []
            streetList = self.streetList
            isObstructed = False

            radNum = int((2*radius*math.pi)/sm.dstep)

            if(radNum<6):
                radNum = 6

            angleStep = (2*math.pi)/radNum


            curAngle = 0 #radians
            curPos = [0,0]
            prevBounds = None
            startPos = None

            for x in range(radNum):
                prevObstructed = isObstructed
                prevPos = curPos[:]
                
                curPos[0] = pos[0] + math.cos(curAngle)*radius + self.getNoise()/100
                curPos[1] = pos[2] + math.sin(curAngle)*radius + self.getNoise()/100
                if(prevPos == [0,0]):
                    prevPos = curPos[:]

                isInBounds = tf.surface.isInBounds(curPos)

                if(prevBounds == None):
                    prevBounds = isInBounds
                
                if(isInBounds== True and prevBounds == False):
                    curPos = tf.surface.findBoundOverlap(curPos, [prevPos[0] - curPos[0], prevPos[1] - curPos[1]])
                    
                elif(prevBounds== True and isInBounds == False):
                    curPos = tf.surface.findBoundOverlap(prevPos, [curPos[0] - prevPos[0], curPos[1] - prevPos[1]])
                    streetPoints = self.finishStreetPart(curPos, streetPoints)

                isObstructed = self.isPointObstructed(curPos)
                if(isObstructed == True and prevObstructed == False):
                    streetPoints = self.finishStreetPart(prevPos, streetPoints)
                if(isObstructed == False and isInBounds == True):
                    streetPoints.append((curPos[0],0,curPos[1]))
                    if(curAngle == 0):
                        startPos = curPos[:]

                prevBounds = isInBounds
                curAngle += angleStep

            if(startPos != None and isInBounds == True and isObstructed == False):
                streetPoints.append((startPos[0],0,startPos[1]))

            if(len(streetPoints)>1):
                streetList.append(streetPoints[:])
                
            for eachStreet in streetList:
                if(len(eachStreet)>2 or (len(eachStreet)==2 and eachStreet[0]!=eachStreet[1])):
                    newStreet = mc.curve(n="roundRoad", d=1, p=eachStreet)
                    sm.streetArray.append(newStreet)
                    print("new Road created: " + str(newStreet) + str(eachStreet))

            mc.refresh()

    def getStartPointsAtBound(self, startx, endx, y, axis, angleRange, tensorType):
        """Returns an array of starting points for streets along a specified bound of the surface.

            startx - (float) - starting position

            endx - (float) - ending position

            y - (float) - second coordinate for both positions

            axis - (string) - either "v" - vertical or "h" - horizontal

            angleRange - (float, float) - creates a starting point only when the tensor at the location has an angle within this range

            tensorType - (string) - either "minor" or "major"

        return - ([float, float], [float, float]...) - array of starting points
        """
        tf = self.tensorField
        point_array = []
        curPos = [0,0]
        prevAngle = 360
        

        if startx == endx:
            return 0
        
        tf = self.tensorField
        sign = (endx-startx)/abs(endx-startx)
        x = startx
        nextX = startx

        if(axis == "v"):
            i_axis = 0
            i_point = 1
            axisAngle = 90

        elif(axis == "h"):
            i_axis = 1
            i_point = 0
            axisAngle = 0

        else:
            print("invalid axis")

        while x*sign < endx*sign:
            curPos[i_point] = x
            curPos[i_axis] = y
            tempTensor = tf.getDirTensor(curPos)
            tensorAngle = getattr(tempTensor, tensorType)
            isWithinBounds = (tensorAngle > angleRange[0] and tensorAngle <= angleRange[1] and x*sign < endx*sign)
            stepRatio = math.sin(math.radians(axisAngle - tensorAngle))

            if(isWithinBounds and (tensorAngle != prevAngle or x*sign >= nextX*sign) and stepRatio != 0):
                point_array.append(curPos[:])
                prevAngle = tensorAngle
                step = abs(self.avgdist / stepRatio)
                nextX = x + step*sign

            x += self.dstep*sign
        return point_array

    def getStartPointMap(self, type):
        """Executes the getStartPointAtBound() function along for both types of street along the relevant bounds.
        
            type - (string) - either "minor" or "major"

        return - ([float, float], [float, float]...) - array of starting points
        """
        tf = self.tensorField
        start_array = []
        prevAngle = 360

        if(type == "major"):
            side1_array = self.getStartPointsAtBound(tf.surface.maxx, tf.surface.minx, tf.surface.minz, "h", (-45, 0), type)
            main_array = self.getStartPointsAtBound(tf.surface.minz, tf.surface.maxz, tf.surface.minx, "v", (-45, 45), type)
            side2_array = self.getStartPointsAtBound(tf.surface.minx, tf.surface.maxx, tf.surface.maxz, "h", (0, 45), type)

        elif(type == "minor"):
            side1_array = self.getStartPointsAtBound(tf.surface.minz, tf.surface.maxz, tf.surface.minx, "v", (45, 90), type)
            main_array = self.getStartPointsAtBound(tf.surface.minx, tf.surface.maxx, tf.surface.maxz, "h", (45, 135), type)
            side2_array = self.getStartPointsAtBound(tf.surface.maxz, tf.surface.minz, tf.surface.maxx, "v", (90, 135), type)
        else:
            print("invalid type")
            return 0

        for x in side1_array:
            start_array.append(x[:])
        for x in main_array:
            start_array.append(x[:])
        for x in side2_array:
            start_array.append(x[:])

        return start_array

    def createMap(self, *pArgs):
        """Main function of the StreetMap Class. Creates a street map.

            *pArgs - Maya UI arguments
        """
        tf = self.tensorField
        ui = self.ui
        self.streetArray = []

        self.majorArray = []    #unused
        self.minorArray = []    #unused

        tf.updatePoints()
        tf.surface.updateBounds()

        if(mc.objExists(ui.streetGroupName)):
            mc.delete(ui.streetGroupName)

        streetDensity = mc.floatSliderGrp(ui.streetDensity, q=True, v=True)
        self.updateDist(self.basedist/streetDensity)

        self.smooth_amount = mc.intField(ui.smoothAmount, q=True, v=True)
        self.noise_amount = mc.floatSliderGrp(ui.noiseAmount, q=True, v=True)
        self.do_branch = mc.checkBox(ui.doEarlyBranch, q=True, v=True)
        self.dstep = mc.floatField(ui.stepSize, q=True, v=True)


        self.streetArray = []
        self.majorArray = []
        self.minorArray = []
        self.radArray = []

        prevSelection = mc.ls(sl=True)

        self.obstructionMap = mc.getAttr(tf.surface.shape + '.' + self.ui.streetAttrName)

        majorStartArray = self.getStartPointMap("major")
        minorStartArray = self.getStartPointMap("minor")
        _progress = ui.Progress('Generating streets...')

        for eachRad in tf.radArray:
            size = eachRad.radius-eachRad.decay*eachRad.radius
            tempSize = size-(1-eachRad.decay)**8*self.avgdist
            pos = [eachRad.pos[0]-5, eachRad.pos[1]-5]
            while(tempSize>self.avgdist/4):
                street = self.Street(self, "rad")
                street.createRadStreet(eachRad.transform, tempSize)
                eachRad.streetSizes.append(tempSize)
                tempSize -= self.avgdist
            _progress.add(10.0/len(tf.radArray))

        pstreet = [(tf.surface.minx,0, tf.surface.minz), (tf.surface.maxx,0, tf.surface.minz)]
        
        for i in majorStartArray:
            street = self.Street(self, "major")  #start horizontally
            street.createDirStreet(i, pstreet, 0)
            pstreet = street.pstreet
            _progress.add(45.0/len(majorStartArray))
        
        pstreet = [(tf.surface.minx, 0, tf.surface.maxz), (tf.surface.minx, 0, tf.surface.minz)]
        
        for i in minorStartArray:
            #print(i)
            street = self.Street(self, "minor")  #start horizontally
            street.createDirStreet(i, pstreet, 0)
            pstreet = street.pstreet
            _progress.add(int(45.0/len(minorStartArray)))

        mc.group(self.streetArray, name = ui.streetGroupName)
        _progress.finish()
        mc.select(prevSelection)

class City:
    """Contains properties and functions connected to building generation.

    __init__ passes in the StreetMap, generateCity() executes all else.
    """

    def __init__(self, sm):
        """ sm - (StreetMap Class) - streetMap that the city generation should use
        """
        self.ui = sm.ui

    def safeProjectCurve(self, offCurve, keepOriginal=True):
        """Projects curve onto the object only if the whole curve lies within the bounds of said object.

            offCurve - (string) - name of the transform node of the curve to project

            keepOriginal - (bool) - if False, deletes the curve after projeting

        Return - [string, string] - transform and shape nodes  of the projected curve
        """
        try:
            projectedOffCurve = mc.polyProjectCurve(offCurve, self.pSurfaceName, d=(0,1,0))[0]
            if(keepOriginal==False):
                mc.delete(offCurve)
            mc.select(projectedOffCurve, hi=True)
            projectedOffCurve = mc.ls(sl=True)
            offCurve = mc.parent(projectedOffCurve[1], w=True)
            self.emptyProjectGroups.append(projectedOffCurve[0][:])
        except:
            mc.warning(str(offCurve)+" could not be projected onto the surface - make sure lays under the whole street map")
        return offCurve

    def prepareStreets(self):
        """Prepares the curves in the group for street generation - 
        if enabled, projects curves and fixes periodic curves.
        
        Returns two lists:

        streetList - a list of object names

        fnStreetList - a list of mFnNurbsCurve objects
        """
        ui = self.ui

        mc.select(ui.streetGroupName, hi=True)
        streetNameList = mc.ls(sl=True, l=True)[1::2]
        #print(streetNameList)

        streetList = []
        fnStreetList = []
        
        if(self.projectOnSurface):
            for street in streetNameList:
                tempStreet = self.safeProjectCurve(street)
                streetList.append(tempStreet[0])
        else:
            streetList = streetNameList[:]

        mc.select(streetList, hi=True)
        streetList = mc.ls(sl=True, l=True)[1::2]
        #print(streetList)
        mc.select(streetList)
        _curSel = om.MGlobal.getActiveSelectionList()

        for x in range(len(streetList)):
            fnStreet = om.MFnNurbsCurve().setObject(_curSel.getDependNode(x))

            if fnStreet.form==3:
                #print("fixing...")
                fixedCurve = mc.offsetCurve(streetList[x], d=0, n=streetList[x], ch=False, sd=0)[0]
                mc.delete(mc.listRelatives(streetList[x], p=True))
                streetList[x] = fixedCurve
                mc.select(streetList[x], hi=True)
                tempSel = om.MGlobal.getActiveSelectionList()
                fnStreet = om.MFnNurbsCurve().setObject(tempSel.getDependNode(1))
            fnStreetList.append(fnStreet)

        return streetList, fnStreetList

    def getStreetSegments(self, fnOffCurve, fnStreetList):
        """Splits the street into segments by finding intersections with other streets

            fnOffCurve - (OpenMaya.MFnNurbsCurve) - a curve offsetted from the street on which buildings are to be generated

            fnStreetList - list of all the other streets
        
        Return - ([float, float], [float, float]...)
            returns a list containing each segment in a format [start length, end length*]
            
            *length along the curve
        """
        lenList = []
        offset = self.offset
        l = 0
        lend = fnOffCurve.length()
        
        startLen = l
        endLen = l

        prevLen = l
        prevDist = l
        prevClosestStreet = []

        step  = offset/10

        while l < lend:
            
            p = fnOffCurve.findParamFromLength(l)
            cvPos = fnOffCurve.getPointAtParam(p)
            dist = 100000
            closestStreet = []

            for eachFnStreet in fnStreetList:
                tempDist = round(eachFnStreet.distanceToPoint(cvPos, space=2),5)
                if(tempDist<dist):
                    dist = tempDist
                    closestStreet = eachFnStreet
            #print(dist)
            
            if(closestStreet == prevClosestStreet):
                if(dist <= offset and prevDist > offset):
                    if(dist== offset):
                        dist-=step
                    dLen = prevDist - dist
                    incr = (prevDist - offset)/dLen
                    endLen = round(prevLen + dLen*incr,5)
                    lenList.append((startLen, endLen))
                    
                        
                elif(dist >= offset and prevDist < offset):
                    if(dist== offset):
                        dist+=step
                    dLen = dist-prevDist
                    incr = (offset - prevDist)/dLen
                    startLen = (prevLen + dLen*incr)
                
            if(l+step>lend and endLen<=startLen):
                lenList.append((startLen,lend))

            prevClosestStreet = closestStreet
            prevDist = dist
            prevLen = l
            l += step

        return lenList

    def getCoordsFromSegments(self, fnOffCurve, lenList, buildingCoords):
        """This function creates coordinates from passed in line segments, 
        so that the buildings are evenly spaced on every street.

        fnOffCurve - (OpenMaya.MFnNurbsCurve) - a curve offsetted from the street on which buildings are to be generated

        lenList - ([float, float], [float, float]...) - a list containing street segments in a format [start length, end length*]

        buildingCoords (OpenMaya.MPoint, ... ) - a list to which the function appends the coordinates

        *length along the curve

        """
        for eachLen in lenList:
            diff = eachLen[1] - eachLen[0]
            numOfHouses = int(diff/(self.buildingSize+self.buildingGap))
            if(numOfHouses>0):
                stepSize = diff/numOfHouses
                l = eachLen[0]
                while l <= eachLen[1]+stepSize/2:
                    p = fnOffCurve.findParamFromLength(l)
                    tempPos = fnOffCurve.getPointAtParam(p)
                    buildingCoords.append(tempPos)
                    l += stepSize

    def collisionFix(self,buildingCoords):
        """This function 
        
        - uses the City's construction dictionary to compare the
        current building coordinates to all the previously made building coordinates.

        - checks the obstruction paint attribute

        - checks if the coordinates are in bounds of the surface

        If any of these are True it removes the coordinates from the buildingCoords list.
        
        buildingCoords - (OpenMaya.MPoint, ... ) - a list of coordinates on which the function operates
        """
        surface = self.ui.surface
        obstructionMap = mc.getAttr(self.ui.surface.shape + '.' + self.ui.streetAttrName)
        removeList = []

        for fnPos in buildingCoords:
            toRemove = False
            buildPos = [fnPos.x, fnPos.z]
            obst = surface.getAttrAtPos(buildPos, obstructionMap)

            if(obst>0.5 or surface.isInBounds(buildPos)==False or surface.isOnBounds(buildPos)==True):
                toRemove = True
            else:
                for eachStreet in self.constructionDict.keys():
                    if(toRemove): break
                    for xPos in self.constructionDict[eachStreet]:
                        usedPos = [xPos.x, xPos.z]
                        sqdist = ppdistance(buildPos, usedPos)
                        #print(sqdist)
                        if sqdist < self.buildingSize**2:
                            toRemove = True
                            break
            if(toRemove):
                removeList.append(fnPos)

        for fnPos in removeList:
            buildingCoords.remove(fnPos)
        
    def getBuildPositions(self, buildingCoords, fnOffCurve, fnStreet):
        """this function is used to execute all the relevant functions for getting
        correct coordinates for all the buildings that are to be placed.

        buildingCoords - (OpenMaya.MPoint, ... ) - a list, which the function fills with the new coordinates

        fnOffCurve - (OpenMaya.MFnNurbsCurve) - a curve offsetted from the street on which buildings are to be generated

        fnStreet - (OpenMaya.MFnNurbsCurve) - the street buildings are generated along
        """
        houseNum = int(fnOffCurve.length()/(self.buildingSize+self.buildingGap))

        fnTempList = self.fnStreetList[:]
        fnTempList.remove(fnStreet)

        lenList = self.getStreetSegments(fnOffCurve, fnTempList)

        self.getCoordsFromSegments(fnOffCurve, lenList, buildingCoords)
        
        self.collisionFix(buildingCoords)

    def generateBuildings(self):
        """Uses the StreetMap's construction dictionary to instance correct buildings to the build possitions.
        """
        if(mc.objExists(self.ui.buildingGroupName)):
            mc.delete(self.ui.buildingGroupName)

        surface = self.ui.surface
        groupDict = self.ui.getGroupDict()

        buildingList = []

        _progress = self.ui.Progress('Generating Buildings...')

        for eachStreet, fnStreet in zip(self.streetList, self.fnStreetList):
            for cvPos in self.constructionDict[eachStreet]:
                buildingSelection = []
                for eachGroup in groupDict.keys():
                    val = surface.getAttrAtPos([cvPos.x, cvPos.z], groupDict[eachGroup].map)
                    rdVal = int(val*10)
                    for i in range(rdVal):
                        buildingSelection.append((random.choice(groupDict[eachGroup].members), groupDict[eachGroup]))
                do = True
                #If the paint value is lower than 10, there is a chance nothing will be placed
                if(len(buildingSelection)<10):
                    randNum = random.randint(1,10)
                    if(randNum > len(buildingSelection)):
                        do = False
                if(do==True):
                    choice = random.randint(0, len(buildingSelection)-1)

                    buildingObj = buildingSelection[choice][0]
                    tempBuilding = mc.instance(buildingObj)
                    buildingList.append(tempBuilding[0])
                    cvStreet = fnStreet.closestPoint(cvPos, tolerance=self.offset)[0]
                    streetNormal = [cvStreet.x-cvPos.x, cvStreet.z-cvPos.z]
                    cosAngle = lcosangle(streetNormal, [0,1])
                    
                    angle = math.degrees(math.acos(cosAngle))
                    if(streetNormal[0]<0):
                        angle *= -1

                    buildGrp = buildingSelection[choice][1]
                    isRandom = mc.checkBox(buildGrp.randomRotation, q=True, v=True)
                    if(isRandom):
                        angle += random.choice([-90,0,90,180])
                    mc.select(tempBuilding)
                    mc.showHidden(tempBuilding)
                    buildingList.append(tempBuilding[0])
                    s = self.buildingSize*self.sizeMult
                    mc.scale(s, s, s)
                    mc.xform(ws=True, t=(cvPos.x, cvPos.y, cvPos.z), ro=(0, angle, 0))
            _progress.add(int(100.0/len(self.constructionDict.keys())))

        if(buildingList!=[]):
            mc.group(buildingList, n=self.ui.buildingGroupName, w=True)
        _progress.finish()

    def generateCity(self, *pArgs):
        """Main City class function. Generates buildings along the streets.

            *pArgs - Maya UI arguments
        """
        mc.setToolTo('selectSuperContext')
        print("Generating city..")
        ui = self.ui
        surface = ui.surface
        sm = ui.streetMap
        tf = ui.tensorField

        self.buildingGap = mc.floatSliderGrp(ui.buildingGap, q=True, v=True)/10
        self.streetWidth = mc.floatSliderGrp(ui.streetWidth, q=True, v=True)/10
        self.buildingSize = mc.floatSliderGrp(ui.buildingSize, q=True, v=True)
        self.projectOnSurface = mc.checkBox(ui.projectCheck, q=True, v=True)
        self.sizeMult = mc.floatSliderGrp(ui.buildSizeMult, q=True, v=True)

        self.offset = self.streetWidth + self.buildingSize
        offset = self.offset
        self.constructionDict = {}
        sideStreetList = []
        

        if(mc.objExists(ui.sideStreetName)):
            mc.delete(ui.sideStreetName)

        if not (mc.objExists(ui.streetGroupName)):
            mc.warning("No streets found to generate on - should be in "+str(ui.streetGroupName))
            return

        if(self.projectOnSurface==True):
            self.emptyProjectGroups = []
            try:
                self.pSurfaceName = mc.textField(ui.projectSurfaceName, q=True, text=True)
                mc.select(self.pSurfaceName)
                mc.select(cl=True)
            except:
                mc.warning("Project to Surface is on but no valid surface was provided")
                return

        buildingList = []
        
        self.streetList, self.fnStreetList = self.prepareStreets()
        #print("Prepared Streets - " +str(self.streetList))
        fnStreetList = self.fnStreetList
        streetList = self.streetList
        _progress = self.ui.Progress('Getting Building Coordinates')

        for (street, fnStreet) in zip(streetList, fnStreetList):

            tempList = streetList[:]
            tempList.remove(street)

            #print("Working on street "+street)

            buildingCoords = []

            for j in [offset, -offset]:

                offCurve = mc.offsetCurve(street, d=j, n= street+"_O", ch=False, cl=True, sd=0)
                #print(offCurve)
                if(self.projectOnSurface):
                    offCurve = self.safeProjectCurve(offCurve, False)
                sideStreetList.append(offCurve[0])

                mc.select(offCurve, hi=True)
                #print(mc.ls(sl=True))
                tempSel = om.MGlobal.getActiveSelectionList()
                fnOffCurve = om.MFnNurbsCurve().setObject(tempSel.getDependNode(1))

                self.getBuildPositions(buildingCoords, fnOffCurve, fnStreet)
                _progress.add(100.0/len(streetList)*2)

            self.constructionDict[street] = buildingCoords[:]
            
        _progress.finish()

        self.generateBuildings() 
        if(sideStreetList!=[]):
            mc.group(sideStreetList, n=ui.sideStreetName)
            mc.hide(ui.sideStreetName)
        if(self.projectOnSurface):
            mc.group(streetList, n=ui.projectedStreetGroupName)
            mc.delete(self.emptyProjectGroups)
            hideObjects([ui.surfaceName,ui.streetGroupName])
        hideObjects([ui.tensorGroupName, ui.tensorPointGroupName, ])

def hideObjects(array, *pArgs):
    """Checks if each object in array exists and the hides it.

        array - [string, string, ...] - array of objects

        *pArgs - Maya UI arguments
    """
    for eachObject in array:
        if(mc.objExists(eachObject)):
            mc.hide(eachObject)

def showObjects(array, *pArgs):
    """Checks if each object in array exists and the shows it if it's hidden.

        array - [string, string, ...] - array of objects

        *pArgs - Maya UI arguments
    """
    for eachObject in array:
        if(mc.objExists(eachObject)):
            mc.showHidden(eachObject)

def deleteObjects(array, *pArgs):
    """Checks if each object in array exists and the deletes it.

        array - [string, string, ...] - array of objects

        *pArgs - Maya UI arguments
    """
    for eachObject in array:
        if(mc.objExists(eachObject)):
            mc.delete(eachObject)

def changePaintAttr(object,type,*pArgs):
    """Sets the tool to attribute paint and switches the current attribute to the one given.

        object - (string) - object name to execute the function on

        type - (string) - name of the paintable attribute

        *pArgs - Maya UI arguments
 
    """
    mc.select(object, hi=True)
    mel.eval('source "artAttrCreateMenuItems.mel"')
    mel.eval('artSetToolAndSelectAttr( "artAttrCtx", "'+str('mesh.'+mc.ls(sl=True, l=True)[1]+'.'+type)+'")')
    mc.select(object)

class UI:
    """Class for handling the UI
    """
    def __init__(self):
        """Sets the names for every single group object the program generates.
        """
        self.surfaceName = "CG_TensorSurface"
        self.tensorGroupName = "CG_TensorField"
        self.tensorName = "Tensor"
        self.tensorPointGroupName = "CG_TensorPoints"
        self.streetGroupName = "CG_StreetCurves"
        self.buildingGroupName = "CG_Buildings"
        self.sideStreetName = "CG_SideStreets"
        self.projectedStreetGroupName = "CG_ProjectedStreets"
        self.buildAttrName = "buildingGroup"

        self.streetAttrName = "obstructionPaint"
        self.buildAttrName = "buildingPaint"
            
    class Progress:
        """A subclass that used for handling executing the progress window easily.
        """
        amount = 0
        name = "Progress"

        def __init__(self, name):
            self.name = name
            mc.progressWindow(title=" ",
                        progress=self.amount,
                        status=self.name,
                        isInterruptable=True )

        def add(self, add):
            """Adds the percentage to the progress amount.
                add - (float) - amount to add
            """
            self.amount += add
            mc.progressWindow( edit=True, progress=self.amount)

        def finish(self):
            """Closes the progress window.
            """
            mc.progressWindow(endProgress=1)

    def saveObjects(self, groupName, *pArgs):
        """Renames the passed group and hides it.
            groupName - (str) - group name

            *pArgs - Maya UI arguments
        """
        if(mc.objExists(groupName)):
            mc.select(groupName)
            mc.rename(groupName+"_S")
            mc.hide()

    def streetGenMenu(self):
        """Street Generation Tab
        """
        surface = self.surface
        tf = self.tensorField
        sm = self.streetMap

        self.streetGenColumn = mc.columnLayout(adj = True, rs=5)
        mc.columnLayout(rs=20)
        mc.setParent("..")
        #----
        mc.rowLayout(nc=2, cat = [(1,'left',10)])
        #-----
        mc.columnLayout(rs = 5, adj = True, w=160)
        mc.iconTextButton( style='iconAndTextHorizontal', image1='putty.png', label='Paint Obstructions', dcc=mc.toolPropertyWindow, c=partial(changePaintAttr, self.surfaceName, self.streetAttrName),ebg = True)
        mc.button(label = "Add Directional Constraint", command=partial(tf.addPoint,'dir'))
        mc.button(label = "Add Radial Constraint", command=partial(tf.addPoint,'rad'))
        mc.setParent('..')

        #-----
        mc.columnLayout()
        mc.rowLayout(nc=2)
        mc.text(l="Street Density ", al='right', w=90)
        self.streetDensity = mc.floatSliderGrp(adj=True, minValue=0.1, maxValue=3.0, value=1.0, s=0.1, w=113, f=True, fieldMaxValue = 10.0)
        mc.setParent('..')
        mc.rowLayout(nc=2)
        mc.text(l="Noise Amount ", al='right', w=90)
        self.noiseAmount = mc.floatSliderGrp(adj=True, minValue=0.0, maxValue=20.0, value=0.0, s=1.0, w=113, f=True, fieldMaxValue = 100.0)
        mc.setParent('..')
        mc.rowLayout(nc=3, h=25)
        mc.text(l="Show Guides ", al='right', w=110)
        mc.checkBox(label = "", ofc = partial(hideObjects, [self.tensorGroupName]), onc = tf.showPreVis, v=False, w=20)
        self.updateTensor = mc.button(label = "Update", command=tf.updatePreVis, en = False, h=15)
        mc.setParent('..')
        mc.rowLayout(nc=2)
        mc.text(l="Show Constraints ", al='right', w=110)
        mc.checkBox(label = "", ofc = partial(hideObjects, [self.tensorPointGroupName]), onc = partial(showObjects, [self.tensorPointGroupName]), v=True, h=20)
        mc.setParent('..')
        mc.setParent('..')

        mc.setParent('..')

        #---
        mc.separator()
        mc.button(label = "Create Streets", command=sm.createMap)
        mc.rowLayout(nc=5)
        mc.button(w=100, h=21, label = "Delete Streets", command=partial(deleteObjects, [self.streetGroupName, self.sideStreetName]))
        mc.button(w=100, h=21, label = "Save Streets", command=partial(self.saveObjects, self.streetGroupName))
        mc.text(l="", al='left', w=20)
        mc.checkBox(l="", ofc = partial(hideObjects, [self.streetGroupName]), onc = partial(showObjects, [self.streetGroupName]), v=True, h=20)
        mc.text(l="Show Streets", al='left', w=110)
        mc.setParent('..')
        mc.separator()

        
        #----
        mc.frameLayout(label = "Advanced", cll=True, cl = True, mw=50, mh=5, cc = "mc.window('"+self.windowID+"', edit = True, rtf=True)")
        
        mc.rowLayout(nc=2, cat = [(1,'left',0), (2, 'left',50)])
        mc.rowLayout(nc=2, cat = [(1,'left',0), (2, 'right',0)])
        mc.text(label = "Smooth Amount")
        self.smoothAmount = mc.intField(minValue=0, maxValue=10, s=1,value=2, w=20)
        mc.setParent('..')
        mc.rowLayout(nc=2, cat = [(1,'left',0), (2, 'right',0)])
        mc.text(label = "Step Size")
        self.stepSize = mc.floatField(minValue=0.01, maxValue=10, value=1.0, w=50)
        mc.setParent('..')
        mc.setParent('..')
        self.doEarlyBranch = mc.checkBox(label = "Connect Early Streets", v=True,)
        
        mc.setParent('..')  #end Advanced
        #---
        mc.setParent('..') #end streetGen
        #--

    def getGroupDict(self):
        """Gets all usable building groups and updates them.

        Return - (dictionary) - dictionary of all usable building groups
        """
        activeGroupDict = {}
        for eachGroup in self.buildGroupDict.keys():
            self.buildGroupDict[eachGroup].map = mc.getAttr(self.surface.shape + '.' + self.buildGroupDict[eachGroup].attrName)
            if(self.buildGroupDict[eachGroup].group != None):
                mc.select(self.buildGroupDict[eachGroup].group, hi=True)
                self.buildGroupDict[eachGroup].members = mc.ls(sl=True)[1::2]

                #print(self.buildGroupDict[eachGroup].members)
                activeGroupDict[eachGroup] = self.buildGroupDict[eachGroup]
            else:
                mc.warning(str(eachGroup)+" couldn't be use because it has no group assigned.")
        return activeGroupDict

    class BuildGroupInfo:
        """A subclass used to create, delete and query values and the side menu for the different building groups.
        """
        def __init__(self, ui, name = None, group = None):
            """ui - (UI Class), name - (string)
            
            group - (string) - object name of the group containing the buildings
            """
            self.name = name
            self.ui = ui
            self.group = group
            self.map = None
            self.members = []

            tag = ui.buildAttrName
            i=0
            while tag in ui.buildGroupDict.keys():
                i += 1
                tag = ui.buildAttrName+str(i)

            self.utg = tag
            
            if(name!=None):
                mc.textScrollList(ui.buildGroupList, e=True, a=name, utg=tag)
                self.attrName = name[:]
            else:
                mc.textScrollList(ui.buildGroupList, e=True, a=tag, utg=tag)
                self.name = tag[:]
                self.attrName = tag[:]
            ui.surface.addPaintAttr(self.attrName)
            
            self.createMenu()

            ui.buildGroupDict[tag] = self
            mc.textScrollList(ui.buildGroupList, e=True, sut=tag)
            ui.selectBuildGroup()

        def updateName(self, *pArgs):
            """Updates the name of the building group in the scrollList.
            """
            newName = mc.textField(self.nameField, q=True, text=True)
            mc.textScrollList(self.ui.buildGroupList, e=True, sut=self.utg)
            oldName = mc.textScrollList(self.ui.buildGroupList, q=True, si = True)
            mc.textScrollList(self.ui.buildGroupList, e=True, ri = oldName)
            mc.textScrollList(self.ui.buildGroupList, e=True, a=newName, utg=self.utg)
            mc.textScrollList(self.ui.buildGroupList, e=True, sut=self.utg)

        def updateGroup(self, *pArgs):
            """Updates the object group name.
            """
            newGroup = mc.textField(self.groupField, q=True, text=True)
            self.group = newGroup
        
        def addToGroup(self, *pArgs):
            """Parents currently selected objects under the group assigned to the current BuildGroupInfo.
            """
            objects = mc.ls(sl=True)

            if(self.group == None):
                mc.warning("No name for group specified!")
                return
            if not(mc.objExists(self.group)):
                mc.group(em=True, n=self.group)

            for eachObject in objects:
                try:
                    newObj = mc.parent(eachObject, w=True)
                except:
                    newObj = eachObject
                print(newObj)
                mc.parent(newObj, self.group)
                

        def visible(self, bool):
            """Changes the visibility of the current building group menu.

                bool - (bool)
            """
            mc.layout(self.column, e=True, vis=bool)

        def cleanUp(self):
            """Deletes the paint attribute associated to this BuildGroupInfo and hides the menu.
            """
            self.ui.surface.removePaintAttr(self.attrName)
            self.visible(False)
            
        def createMenu(self):
            """Menu layout for BuildGroupInfo.
            """
            mc.setParent(self.ui.buildGroupColumn)
            self.column = mc.columnLayout(vis=False, adj=True)
            mc.rowLayout(nc=2)
            mc.text(l="Name: ", w=50)
            self.nameField = mc.textField(w=120, text=self.name, cc=self.updateName)
            mc.setParent('..')
            mc.rowLayout(nc=2, h=30)
            mc.text(l="Group: ", w=50)
            self.groupField = mc.textField(w=120, text=self.group, cc=self.updateGroup)
            mc.setParent('..')
            mc.button(l="Add Selection to Group", c=self.addToGroup, )
            mc.columnLayout(h=5)
            mc.setParent('..')
            mc.iconTextButton( style='iconAndTextHorizontal', 
                image1='putty.png', label='Paint Building Density', dcc=mc.toolPropertyWindow, 
                c=partial(changePaintAttr, self.ui.surfaceName, self.attrName),ebg = True)
            mc.columnLayout(h=10)
            mc.setParent('..')
            mc.rowLayout(nc=2)
            mc.text(l="Random Building Rotation ", w=155, al='right')
            self.randomRotation = mc.checkBox(l="")
            mc.setParent('..')
    
    def addBuildGroup(self, *pArgs):
        """Creates a new building group.
        """
        newBuildGroup = self.BuildGroupInfo(self)

    def removeBuildGroup(self, *pArgs):
        """Removes a building group from the list and cleans up.
        """
        tag = mc.textScrollList(self.buildGroupList, q=True, sut=True)

        selectedItem = mc.textScrollList(self.buildGroupList, q=True, si=True)
        mc.textScrollList(self.buildGroupList, e=True, ri=selectedItem)

        tag = str(tag).split("'")[1]
        self.buildGroupDict[tag].cleanUp()

        self.buildGroupDict.pop(tag)
        if(self.curSelectedItem==tag):
            self.curSelectedItem = None

    def selectBuildGroup(self, *pArgs):
        """Switches the menu and selected attribute according to the selected item in the scrollLayout.
        """
        if(self.curSelectedItem != None):
            self.buildGroupDict[self.curSelectedItem].visible(False)
        selectedItem = mc.textScrollList(self.buildGroupList, q=True, sut=True)
        tag = str(selectedItem).split("'")[1]

        self.buildGroupDict[tag].visible(True)

        self.curSelectedItem = tag

        changePaintAttr(self.surfaceName, self.buildGroupDict[tag].attrName)

    def boolSurfaceProjection(self, bool, *pArgs):
        """Un/greys out the text field for surface projection name.

            bool - (bool)
        """
        mc.textField(self.projectSurfaceName, edit=True, ed=bool)
    
    def cityGenMenu(self):
        """City Generation Tab
        """
        surface = self.surface
        tf = self.tensorField
        sm = self.streetMap
        city = self.city

        self.buildGroupDict = {}
        self.curSelectedItem = None

        self.cityGenColumn = mc.columnLayout(adjustableColumn = True,  columnAttach = ('both', 5),rs=5)
        mc.columnLayout(rs=20)
        mc.setParent("..")

        mc.columnLayout(rs=5)
        self.buildingSize = mc.floatSliderGrp(adj = True, w = 340, l="Building Size ",minValue=0.1, maxValue=5.0, fieldMaxValue = 100.0, s=0.2,value=1.0, f = True)
        self.buildingGap = mc.floatSliderGrp(adj = True, w = 340, l="Building Gap ",minValue=0.0, maxValue=5.0, fieldMaxValue = 100.0, s=0.2,value=2.0, f = True)
        self.streetWidth = mc.floatSliderGrp(adj = True, w = 340, l="Street Width ",minValue=0.0, maxValue=5.0, fieldMaxValue = 100.0, s=0.2,value=2.0, f = True)
        mc.setParent('..')
        mc.rowLayout(nc=3)
        mc.text(l="Project On Surface ", al='right', w=191)
        self.projectCheck = mc.checkBox(label = "", v=False, w=20, onc=partial(self.boolSurfaceProjection, True), ofc=partial(self.boolSurfaceProjection, False))
        self.projectSurfaceName = mc.textField(w=120, ed=False)
        mc.setParent('..')
        mc.separator()
        mc.paneLayout( configuration='vertical2', h = 150)
        mc.columnLayout(w=60, h = 120, adj=True)
        mc.rowLayout(nc=2)
        mc.button(l="Add Group",  c=self.addBuildGroup, w=80)
        mc.button(l="Remove Group",  c=self.removeBuildGroup, w=100)
        mc.setParent('..')
        self.buildGroupList = mc.textScrollList( sc = self.selectBuildGroup)
        mc.setParent('..')

        self.buildGroupColumn = mc.columnLayout(w=60, adj=True)
        mc.setParent('..')

        mc.setParent('..')
        mc.separator()
        mc.button(label = "Generate", command=city.generateCity)
        mc.rowLayout(nc=5)
        mc.button(w=110, h=21, label = "Delete Buildings", command=partial(deleteObjects, [self.buildingGroupName, self.projectedStreetGroupName, self.sideStreetName]))
        mc.button(w=100, h=21, label = "Save Buildings", command=partial(self.saveObjects, self.buildingGroupName))
        mc.text(l="", al='left', w=20)
        mc.checkBox(l="", ofc = partial(hideObjects, [self.buildingGroupName]), onc = partial(showObjects, [self.buildingGroupName]), v=True, h=20)
        mc.text(l="Show Buildings", al='left', w=110)
        mc.setParent('..')
        mc.separator()


        mc.frameLayout(label = "Advanced", cll=True, cl = True, mw=30, mh=5, cc = "mc.window('"+self.windowID+"', edit = True, rtf=True)")
        
        mc.rowLayout(nc=2, cat = [(1,'left',0), (2, 'left',0)])
        mc.text(label = "Mesh Scale Multiplier: ")
        self.buildSizeMult = mc.floatSliderGrp(minValue=0.01, maxValue=2.0,fieldMaxValue=10.0, s=0.01,value=1.0,w = 200, field=True)
        mc.setParent('..')

        mc.setParent('..')

        mc.columnLayout(rs=20)
        mc.setParent('..')
        mc.setParent("..")
        mc.setParent('..')

    def prepareProgram(self, *pArgs):
        """This function
        - checks if the test building groups exist and appends them to the scrollLayout
        - checks for tensor points from previous sessions of the program and initializes them
        - prepares the surface plane for street generation
        """
        surface = self.surface
        tf = self.tensorField
        sm = self.streetMap
        city = self.city

        surface.prepareSurface()
        mc.setToolTo('selectSuperContext')
        try:
            try1 = mc.file(scriptPath+"/CG_buildings.fbx", i=True, mnc=True)
            print("Building file successfuly loaded!")
        except:
            pass
            
        if(mc.objExists("CG_smallBuildings")):
            mc.hide("CG_smallBuildings")
            newBuildGroup = self.BuildGroupInfo(self, "smallBuildings", "CG_smallBuildings")
        if(mc.objExists("CG_midBuildings")):
            mc.hide("CG_midBuildings")
            newBuildGroup = self.BuildGroupInfo(self, "midBuildings", "CG_midBuildings")
        if(mc.objExists("CG_tallBuildings")):
            mc.hide("CG_tallBuildings")
            newBuildGroup = self.BuildGroupInfo(self, "tallBuildings", "CG_tallBuildings")
        
        mc.select(cl=True)
        

        mc.layout(self.tabs, edit=True, en = True)
        mc.button(self.start, edit= True, en = False)

    def mainMenu(self):
        """Main Window Layout
        """
        self.surface = Surface(self)
        self.tensorField = TensorField(self.surface)
        self.streetMap = StreetMap(self.tensorField)
        self.city = City(self.streetMap)

        surface = self.surface
        tf = self.tensorField
        sm = self.streetMap
        city = self.city


        self.window = mc.window(self.windowID, title="City Generator 1.0")
        mc.window(self.window, edit = True, s=False, rtf=True, w=350, h=200)

        mc.columnLayout(adj = True,  columnAttach = ('both', 5),rs=10)

        #-
        self.mainImage = mc.image( image = scriptPath + '/CG_thumbnail.png')
        self.start = mc.button(label = "Start City Generation", command=self.prepareProgram)

        #--
        self.tabs = mc.tabLayout(innerMarginWidth=5, innerMarginHeight=20, w=350, en=False)

        self.streetGenMenu()

        self.cityGenMenu()

        mc.tabLayout( self.tabs, edit=True, tabLabel=((self.streetGenColumn, "Street Generation"), (self.cityGenColumn, "Building Generation")) )

        mc.setParent('..')  

        mc.showWindow(self.window)

    def deleteStart(self, exObjects, *pArgs):
        """Deletes all the groups from the previous session and opens the main menu.

            exObjects - [string, string, ...] - previous objects
        """

        for eachObject in exObjects:
            if(mc.objExists(eachObject)):
                mc.delete(eachObject)
        mc.deleteUI(self.checkID)
        self.mainMenu()

    def fixStart(self, exObjects, *pArgs):
        """Renames all the groups from previous session, hides them and opens the main menu.

            exObjects - [string, string, ...] - previous objects
        """
        dupObjects = []
        for eachObject in exObjects:
            dupObject = mc.duplicate(eachObject)[0]
            dupObjects.append(dupObject)
        deleteObjects(exObjects)
        print(dupObjects)
        mc.hide(dupObjects)
        mc.deleteUI(self.checkID)
        self.mainMenu()

    def keepStart(self, *pArgs):
        """Keeps all the objects from the previous session and opens the main menu.
        """
        mc.deleteUI(self.checkID)
        self.mainMenu()

    def checkForPrevRun(self):
        """Checks if there are any objects that match the naming of the different groups this program creates.
        """
        exObjects = []
        nameArray = [self.projectedStreetGroupName, self.sideStreetName, self.buildingGroupName, self.surfaceName, self.tensorGroupName, self.tensorName, self.tensorPointGroupName, self.streetGroupName]
        for eachName in nameArray:
            if(mc.objExists(eachName)):
                exObjects.append(eachName)

        if(exObjects != []):
                window2 = mc.window(self.checkID, s=False, title = "CityGen 1.0")
                mc.window(window2, edit=True, widthHeight=(300,100))
                mc.columnLayout(adjustableColumn = True,  columnAttach = ('both', 10))
                mc.text(label = "\nFound objects from a previous run of this program.\n Save a copy or delete?")
                mc.columnLayout(h=10)
                mc.setParent('..')
                mc.separator()
                mc.columnLayout(h=10)
                mc.setParent('..')
                mc.rowLayout(nc=3)
                mc.button(w=95, label = "Delete", command=partial(self.deleteStart, exObjects))
                mc.button(w=95, label = "Save a Copy", command=partial(self.fixStart, exObjects))
                mc.button(w=95, label = "Keep", command=self.keepStart)
                mc.setParent('..')
                mc.showWindow(window2)
        else:
            self.mainMenu()

    def start(self):
        """Main Program Function
        """
        self.windowID = "cityGenWindow"
        if mc.window(self.windowID, exists=True):
            mc.deleteUI(self.windowID)
        self.checkID = "checkWindow"
        if mc.window(self.checkID, exists=True):
            mc.deleteUI(self.checkID)

        self.checkForPrevRun()

def Start():

    if mc.artAttrCtx('cityPaintCtx', query=True, exists=True) == False:
        mc.artAttrCtx('cityPaintCtx')
        mc.select(cl = True)

    newMenu = UI()
    newMenu.start()

#main program
if __name__ =="__main__":
    Start()   
