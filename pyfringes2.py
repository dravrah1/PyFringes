import sys
import numpy as np
import math
import time
from PyQt5 import QtCore, QtGui, uic, QtWidgets
import zernike as Z

qtCreatorFile = "pyfringesGUI-R1.ui" # Enter file here.
img = QtGui.QImage()
zernarrays = np.zeros((2,2,2))
glbmask = np.zeros((2,2))
havezernarray = False
priorvalue = []

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class MyApp(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.GenFringes.clicked.connect(self.GenerateFringeImage)
        self.piston.editingFinished.connect(self.coeffchanged)
        self.tiltx.editingFinished.connect(self.coeffchanged)
        self.tilty.editingFinished.connect(self.coeffchanged)
        self.power.editingFinished.connect(self.coeffchanged)
        self.astigx.editingFinished.connect(self.coeffchanged)
        self.astigy.editingFinished.connect(self.coeffchanged)
        self.comax.editingFinished.connect(self.coeffchanged)
        self.comay.editingFinished.connect(self.coeffchanged)
        self.spherical.editingFinished.connect(self.coeffchanged)   
        self.numPix.editingFinished.connect(self.Nchanged)
        self.shearplateCheck.stateChanged.connect(self.shearcheck)
        self.shearplatecomboBox.currentIndexChanged.connect(self.shearcombo)
        self.UpdateShear.clicked.connect(self.GenerateShearImage)
        self.platethickness.editingFinished.connect(self.shearparamchanged)
        self.platewedge.editingFinished.connect(self.shearparamchanged)
        self.plateindex.editingFinished.connect(self.shearparamchanged)
        self.platesize.editingFinished.connect(self.shearparamchanged)
        self.beamsize.editingFinished.connect(self.shearparamchanged)
#        self.frame_4.setHidden(True)
        self.frame_5.setHidden(True)
    
    def resizeEvent(self, event):
        if havezernarray: 
            self.resizeimgs()

    def resizeimgs(self):
        GView = self.graphicsView
        GView.setScene(self.scene)
 #       self._pixmapHandle = None
        GView.aspectRatioMode = QtCore.Qt.KeepAspectRatio
        GView.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        GView.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        GView.zoomStack = []
        GView.canZoom = True
        GView.canPan = True
        
        self.updateViewer(GView)

        GView2 = self.graphicsView2
        GView2.setScene(self.scene2)
 #       self._pixmapHandle2 = None
        GView2.aspectRatioMode = QtCore.Qt.KeepAspectRatio
        GView2.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        GView2.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        GView2.zoomStack = []
        GView2.canZoom = True
        GView2.canPan = True
        
        self.updateViewer(GView2)

    def coeffchanged(self):
        coeffedit = self.sender()
#        print self.sender().objectName()
#        print self.sender().isModified()
        if coeffedit.isModified():
            self.GenerateFringeImage()
            coeffedit.setModified(False)

    def shearparamchanged(self):
        shearedit = self.sender()
        if shearedit.isModified():
            self.shearcheck()
            shearedit.setModified(False)
        
    def Nchanged(self):
        global havezernarray
        havezernarray = False
        self.coeffchanged()

    def shearcheck(self):
        if self.shearplateCheck.isChecked():
            self.graphicsView2.setHidden(False)
            self.GenerateFringeImage()
            self.GenerateShearImage()
        else:
            self.graphicsView2.setHidden(True)
            self.GenerateFringeImage()
            self.resizeimgs()
   
    def shearcombo(self):
        select = self.shearplatecomboBox.currentIndex()
#        print select
#        print "selected" 
        thick_o = float(self.platethickness.text())
        wedge_o = float(self.platewedge.text())
        index_o = float(self.plateindex.text())
        size_o = float(self.platesize.text())
        beam_o = float(self.beamsize.text())
        beam = beam_o
        if select == 0:
            thick = 0.75
            wedge = 117
            index = 1.457
            size = 5.28
            if beam_o > 3:
                beam = 3
        elif select == 1:
            thick = 1.3
            wedge = 83
            index = 1.457
            size = 7.94
            if beam_o > 5:
                beam = 5
        elif select == 2:
            thick = 2.6
            wedge = 40
            index = 1.457
            size = 15.6
            if beam_o > 10:
                beam = 10
        elif select == 3:
            thick = 6.35
            wedge = 18
            index = 1.457
            size = 38
            if beam_o > 25.4:
                beam = 25.4
        elif select == 4:
            thick = 13
            wedge = 10
            index = 1.457
            size = 78
            if beam_o > 50:
                beam = 50
        self.platethickness.setText(str(thick))
        self.platewedge.setText(str(wedge))
        self.plateindex.setText(str(index))
        self.platesize.setText(str(size))
        self.beamsize.setText(str(beam))
        self.GenerateShearImage()
                            
#blah = QtWidgets.QLineEdit
#blah.setText()
#blah = QtWidgets.QProgressBar
#blah.setValue

#    def framechange(self,self.frame.changeEvent):
#        QtWidgets.QMessageBox.information(self,"Information!","Frame has been changed...")

#    def resizeEvent(self,resizeEvent):
#        QtWidgets.QMessageBox.information(self,"Information!","Window has been resized...")
        
    def CircleMask(self,N,d):
        mask = np.zeros((N,N))
        for x in range(N):
            for y in range(N):
                r=math.sqrt((float(x)-N/2)**2+(float(y)-N/2)**2)
                if r <= d/2.0:
                    mask[x,y]=1
        mask = np.uint8(mask)
        return mask

   
    def GeneratePhase(self):
        global zernarrays
        global havezernarray
        
        nterms=int(self.numterms.text())
        wavelength = float(self.wavelength.text())*(10.0**(-6))
        zernterms = Z.genZern(nterms)
        N = int(self.numPix.text())
        N = int(int(float(N)/4.0)*4)     
        if havezernarray == False:
            zernarrays = Z.genZernArrays(nterms,N,zernterms,self.progressBar)
            havezernarray = True
        zcoeff = np.empty(nterms)
        try:
            zcoeff[0] = float(self.piston.text())
            zcoeff[1] = float(self.tiltx.text())
            zcoeff[2] = float(self.tilty.text())
            zcoeff[3] = float(self.power.text())
            zcoeff[4] = float(self.astigx.text())
            zcoeff[5] = float(self.astigy.text())
            zcoeff[6] = float(self.comax.text())
            zcoeff[7] = float(self.comay.text())
            zcoeff[8] = float(self.spherical.text())
        except ValueError as e:
            mbox = QtWidgets.QMessageBox()
            buttonreply = mbox.question(self, 'PyFringes Error Message','Bad Coefficent Value', mbox.Ok)
            return np.zeros(1)
        mask = self.CircleMask(N,N)
        global glbmask
        glbmask = mask
        phsimg = np.zeros((N,N))
        for i in range(nterms):
            phsimg = phsimg + zcoeff[i] * zernarrays[i,:,:] 
        phsimgraw = phsimg.transpose()
        phsimg = phsimg.transpose() * mask

        phsmean = phsimg.sum()/mask.sum()
        phsimgtmp = (phsimg - phsmean) * mask
        
        PV = phsimgtmp.max()-phsimgtmp.min()
        self.txtPV.setText(str(PV))
        RMS = phsimgtmp.std() * math.sqrt(float(N**2)/mask.sum())
        self.txtRMS.setText(str(RMS))
        return (phsimg*2*math.pi)

    def GenerateShearImage(self):
        N = int(self.numPix.text())
        N = int(int(float(N)/4.0)*4)
        wavelength = float(self.wavelength.text())*(10.0**(-6))
        zernphs = self.GeneratePhase()
        if zernphs.size == 1: return
#        print zernphs.sum()
        try:
            thick = float(self.platethickness.text())
            wedge = float(self.platewedge.text())/3600.0*math.pi/180.0
            size = float(self.platesize.text())
            index = float(self.plateindex.text())
            beamsz = float(self.beamsize.text())
        except ValueError as e:
            mbox = QtWidgets.QMessageBox()
            buttonreply = mbox.question(self, 'PyFringes Error Message','Bad Shear Parameter Value', mbox.Ok)
            return 0       
         
        alpha = index*math.sin(2*wedge)
        dist = size/2/math.sqrt(2)*1.25
        shear = 2**(0.5) * thick * math.tan(1/index * math.sin(math.pi/4))
        pix = beamsz/N
        shearpix=int(shear/pix)
#        mask = self.CircleMask(N,N)
        global glbmask
        mask = glbmask
        N2=N+shearpix
        pad=np.zeros((N, shearpix))
        I = np.zeros((N,N2))
        phs1 = np.append(zernphs,pad,axis=1)
        phs2 = np.append(pad,zernphs,axis=1)
        int1 = np.append(mask,pad,axis=1)
        int2 = np.append(pad,mask,axis=1)
        tanalphapix = math.tan(alpha)*(pix*10**(-3))/wavelength * 2 * math.pi
        for x in range(N):
            for y in range(N2): 
                I[x,y] = int1[x,y] + int2[x,y] + (2 * math.sqrt(int1[x,y] * int2[x,y]) * math.cos(phs1[x,y]-phs2[x,y]+tanalphapix*x ))
                self.progressBar.setValue(int(round(float(x*N2+y)/float(N2*N)*100.0)))
                self.progressBar.update()
        I2 = I/4*255
        #make sure image is 32 bit aligned
        if round(float(N2)/4.0 + 0.5) != float(N2)/4.0:
            N3=round(float(N2)/4.0 + 0.5)*4-N2
            pad2 = np.zeros((N, int(N3)))
            I2 = np.append(I2,pad2,axis=1)
        # add shear plate line
        I2[N/2,:]=0
        I3 = np.uint8(I2)
        imgdata = I3
        wid = int(imgdata.shape[1])
        ht = int(imgdata.shape[0])
        im = QtGui.QImage(imgdata.data, wid, ht, QtGui.QImage.Format_Grayscale8)
        pixmap = QtGui.QPixmap.fromImage(im)

        GView = self.graphicsView2
    
        self.scene2 = QtWidgets.QGraphicsScene()
        GView.setScene(self.scene2)
        self._pixmapHandle2 = None
        
        GView.aspectRatioMode = QtCore.Qt.KeepAspectRatio

        # Scroll bar behaviour.
        #   Qt.ScrollBarAlwaysOff: Never shows a scroll bar.
        #   Qt.ScrollBarAlwaysOn: Always shows a scroll bar.
        #   Qt.ScrollBarAsNeeded: Shows a scroll bar only when zoomed.
        GView.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        GView.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)

        # Stack of QRectF zoom boxes in scene coordinates.
        GView.zoomStack = []

        # Flags for enabling/disabling mouse interaction.
        GView.canZoom = True
        GView.canPan = True
        
        self._pixmapHandle2 = self.scene2.addPixmap(pixmap)

        GView.setSceneRect(QtCore.QRectF(pixmap.rect()))

        self.updateViewer(GView)
                    
    
     
    def GenerateFringeImage(self):
        GView = self.graphicsView
        # generate zernike terms
        nterms=int(self.numterms.text())
        wavelength = float(self.wavelength.text())*(10.0**(-6))
#        zernterms = Z.genZern(nterms)
#        zcoeff = np.empty(nterms)
#        zcoeff[0] = float(self.piston.text())
#        zcoeff[1] = float(self.tiltx.text())
#        zcoeff[2] = float(self.tilty.text())
#        zcoeff[3] = float(self.power.text())
#        zcoeff[4] = float(self.astigx.text())
#        zcoeff[5] = float(self.astigy.text())
#        zcoeff[6] = float(self.comax.text())
#        zcoeff[7] = float(self.comay.text())
#        zcoeff[8] = float(self.spherical.text())
        N = int(self.numPix.text())
        N = int(int(float(N)/4.0)*4)
        zernimg = np.zeros((N,N))
        zernimg = np.uint8(zernimg)   
#        mask = self.CircleMask(N,N)
        phase = self.GeneratePhase()
        global glbmask
        mask = glbmask        
#        for x in range(N):
#            for y in range(N):
#               if mask[x,y] != 0:
#                    xx=float(x)-N/2.0
#                    yy=float(y)-N/2.0
#                    p = ((xx/(N/2))**2 + (yy/(N/2))**2)**0.5
#                    theta = math.atan2(yy,xx)
#                    phase = Z.calcZern(zernterms,zcoeff,nterms,p,theta)
#                    fringe = (math.cos(phase*2*3.14)+1)/2*255
#                    val2=int(fringe)
#                    zernimg[x,y]=val2

        if phase.size == 1: return
        fringe = (np.cos(phase)+1)/2*255
        fringe = mask * fringe
        zernimg = np.uint8(fringe)
        #print zernimg
#       zernimg2 = np.multiply(zernimg,mask)
        im = QtGui.QImage(zernimg.data, N, N, QtGui.QImage.Format_Grayscale8)
        pixmap = QtGui.QPixmap.fromImage(im)
    
        self.scene = QtWidgets.QGraphicsScene()
        GView.setScene(self.scene)
        self._pixmapHandle = None
        
        GView.aspectRatioMode = QtCore.Qt.KeepAspectRatio

        # Scroll bar behaviour.
        #   Qt.ScrollBarAlwaysOff: Never shows a scroll bar.
        #   Qt.ScrollBarAlwaysOn: Always shows a scroll bar.
        #   Qt.ScrollBarAsNeeded: Shows a scroll bar only when zoomed.
        GView.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        GView.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)

        # Stack of QRectF zoom boxes in scene coordinates.
        GView.zoomStack = []

        # Flags for enabling/disabling mouse interaction.
        GView.canZoom = True
        GView.canPan = True
        
        self._pixmapHandle = self.scene.addPixmap(pixmap)

        GView.setSceneRect(QtCore.QRectF(pixmap.rect()))

        self.updateViewer(GView)

        if self.shearplateCheck.isChecked():
            self.GenerateShearImage()
    

    def updateViewer(self,GV):
        if not self.hasImage():
            return
        if len(GV.zoomStack) and GV.sceneRect().contains(GV.zoomStack[-1]):
              GV.fitInView(GV.zoomStack[-1], QtCore.Qt.IgnoreAspectRatio)  # Show zoomed rect (ignore aspect ratio).
        else:
             GV.zoomStack = []  # Clear the zoom stack (in case we got here because of an invalid zoom).
             GV.fitInView(GV.sceneRect(), GV.aspectRatioMode)  # Show entire image (use current aspect ratio mode).
            

    def hasImage(self):
        """ Returns whether or not the scene contains an image pixmap.
        """
        return self._pixmapHandle is not None

    def clearImage(self):
        """ Removes the current image pixmap from the scene if it exists.
        """
        if self.hasImage():
            self.scene.removeItem(self._pixmapHandle)
            self._pixmapHandle = None

    def pixmap(self):
        """ Returns the scene's current image pixmap as a QPixmap, or else None if no image exists.
        :rtype: QPixmap | None
        """
        if self.hasImage():
            return self._pixmapHandle.pixmap()
        return None

    def image(self):
        """ Returns the scene's current image pixmap as a QImage, or else None if no image exists.
        :rtype: QImage | None
        """
        if self.hasImage():
            return self._pixmapHandle.pixmap().toImage()
        return None

    def setImage(self, image):
        """ Set the scene's current image pixmap to the input QImage or QPixmap.
        Raises a RuntimeError if the input image has type other than QImage or QPixmap.
        :type image: QImage | QPixmap
        """
        if type(image) is QtGui.QPixmap:
            pixmap = image
        elif type(image) is QtGui.QImage:
            pixmap = QtGui.QPixmap.fromImage(image)
        else:
            raise RuntimeError("ImageViewer.setImage: Argument must be a QImage or QPixmap.")
        if self.hasImage():
            self._pixmapHandle.setPixmap(pixmap)
        else:
            self._pixmapHandle = self.scene.addPixmap(pixmap)
        self.setSceneRect(QtCore.QRectF(pixmap.rect()))  # Set scene size to image size.
        self.updateViewer()

            
            
            
def hasImage(self):
        """ Returns whether or not the scene contains an image pixmap.
        """
        return self._pixmapHandle is not None





if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MyApp()
    window.show()
#    window.graphicsView2.setHidden(True) 
    sys.exit(app.exec_())
    