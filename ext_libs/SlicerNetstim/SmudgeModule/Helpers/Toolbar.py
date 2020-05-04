import qt, vtk, slicer
from qt import QToolBar
import os
from slicer.util import VTKObservationMixin

import SmudgeModule
import ImportAtlas
import ImportSubject
import TransformsUtil
from . import WarpEffect

class reducedToolbar(QToolBar, VTKObservationMixin):

  def __init__(self):

    QToolBar.__init__(self)
    VTKObservationMixin.__init__(self)

    self.parameterNode = SmudgeModule.SmudgeModuleLogic().getParameterNode()
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateToolbarFromMRML)
    
    self.setWindowTitle(qt.QObject().tr("LeadDBS"))

    smw = slicer.util.mainWindow()
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    layoutManager = slicer.app.layoutManager()

    self.addObserver(interactionNode, interactionNode.InteractionModeChangedEvent, self.onInteractionModeChanged) 


    #
    # Layout
    #

    layoutFourUpAction = qt.QAction(smw)
    layoutFourUpAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','LayoutFourUpView.png'))))
    layoutFourUpAction.connect('triggered(bool)', lambda t: layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView))
    self.addAction(layoutFourUpAction)

    layoutTabbedAction = qt.QAction(smw)
    layoutTabbedAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','LayoutTabbedSliceView.png'))))
    layoutTabbedAction.connect('triggered(bool)', lambda t: layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutTabbedSliceView))

    self.addAction(layoutTabbedAction)

    sliceIntersectionAction = qt.QAction(smw)
    sliceIntersectionAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','SlicesCrosshair.png'))))
    sliceIntersectionAction.setCheckable(True)
    sliceIntersectionAction.connect('toggled(bool)', self.sliceIntersectionToggle)

    self.addAction(sliceIntersectionAction)
    
    #
    # Window Level
    #

    windowLevelModeActions = qt.QActionGroup(smw)
    windowLevelModeActions.setExclusive(True)

    windowLevelAdjustModeAction = qt.QAction(smw)
    windowLevelAdjustModeAction.setText('Adjust')
    windowLevelAdjustModeAction.setCheckable(True)
    windowLevelAdjustModeAction.setChecked(True)
    windowLevelAdjustModeAction.connect('triggered(bool)', lambda t: interactionNode.SetAttribute('AdjustWindowLevelMode', 'Adjust'))

    windowLevelRegionModeAction = qt.QAction(smw)
    windowLevelRegionModeAction.setText('Select Region')
    windowLevelRegionModeAction.setCheckable(True)
    windowLevelRegionModeAction.connect('triggered(bool)', lambda t: interactionNode.SetAttribute('AdjustWindowLevelMode', 'Rectangle'))

    windowLevelCenteredRegionModeAction = qt.QAction(smw)
    windowLevelCenteredRegionModeAction.setText('Select Region - centered')
    windowLevelCenteredRegionModeAction.setCheckable(True)
    windowLevelCenteredRegionModeAction.connect('triggered(bool)', lambda t: interactionNode.SetAttribute('AdjustWindowLevelMode', 'RectangleCentered'))

    windowLevelModeActions.addAction(windowLevelAdjustModeAction)
    windowLevelModeActions.addAction(windowLevelRegionModeAction)
    windowLevelModeActions.addAction(windowLevelCenteredRegionModeAction)

    windowLevelMenu = qt.QMenu(smw)
    windowLevelMenu.addActions(windowLevelModeActions.actions())
    
    self.windowLevelAction = qt.QAction(smw)
    self.windowLevelAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','MouseWindowLevelMode.png'))))
    self.windowLevelAction.setCheckable(True)
    self.windowLevelAction.setMenu(windowLevelMenu)
    self.windowLevelAction.connect('toggled(bool)', lambda t: interactionNode.SetCurrentInteractionMode(5) if t else interactionNode.SetCurrentInteractionMode(2))

    self.addAction(self.windowLevelAction)

    #
    # Warp visible in slice view
    #

    warpViewAction = qt.QAction(smw)
    warpViewAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','GlyphIcon.png'))))
    warpViewAction.setCheckable(True)
    warpViewAction.connect('toggled(bool)', self.onWarpViewAction)
    self.addAction(warpViewAction)


    #
    # Modality
    #
    self.modalityComboBox = qt.QComboBox()
    self.modalityComboBox.addItem('t1')
    self.modalityComboBox.view().pressed.connect(self.onModalityPressed)
    self.addWidget(self.modalityComboBox)

    #
    # B <-> F slider
    #

    templateSlider = qt.QSlider(1)
    templateSlider.singleStep = 10
    templateSlider.minimum = 0
    templateSlider.maximum = 100
    templateSlider.value = 0
    templateSlider.setFixedWidth(120)
    templateSlider.connect('valueChanged(int)', lambda value: slicer.util.setSliceViewerLayers(foregroundOpacity = value / 100.0))
    self.addWidget(templateSlider)

    #
    # Resolution
    #
    self.addWidget(qt.QLabel('Warp Resolution: '))
    self.resolutionComboBox = qt.QComboBox()
    avalibaleResolutions = [0.5, 1, 2, 5, 10]
    self.resolutionComboBox.addItems([str(r)+'mm' for r in avalibaleResolutions])
    self.resolutionComboBox.setCurrentIndex(avalibaleResolutions.index(float(self.parameterNode.GetParameter("resolution"))))
    self.resolutionComboBox.connect('currentIndexChanged(int)', self.onResolutionChanged)

    self.addWidget(self.resolutionComboBox)

    #
    # Space Separator
    #

    empty = qt.QWidget()
    empty.setSizePolicy(qt.QSizePolicy.Expanding,qt.QSizePolicy.Preferred)
    self.addWidget(empty)

    #
    # Subject
    #

    self.subjectNameLabel = qt.QLabel('Subject: ')    
    self.addWidget(self.subjectNameLabel)

    #
    # Save
    #
    self.saveButton = qt.QPushButton("Finish and Exit")
    self.saveButton.setFixedWidth(200)
    self.saveButton.setStyleSheet("background-color: green")
    self.addWidget(self.saveButton)
    ImportAtlas.ImportAtlasLogic().run(os.path.join(self.parameterNode.GetParameter("MNIAtlasPath"), 'DISTAL Minimal (Ewert 2017)'))
    self.saveButton.connect("clicked(bool)", self.onSaveButton)

    #
    # Update
    #


    self.updateModalities(self.parameterNode.GetParameter("subjectPath"))
    reducedToolbarLogic().loadSubjectTransforms()
    self.onModalityPressed([],self.modalityComboBox.currentText)
    self.updateToolbarFromMRML()


  def onWarpViewAction(self, t):
    warpID = self.parameterNode.GetParameter("warpID")
    if warpID != "":
      warpNode = slicer.util.getNode(warpID)
      warpNode.GetDisplayNode().SetVisibility(t)
      warpNode.GetDisplayNode().SetVisibility2D(t)

  def initializeTransforms(self, imageNode):
    glanatCompositeNode = slicer.util.getNode(self.parameterNode.GetParameter("glanatCompositeID"))
    # apply glanat to image
    imageNode.SetAndObserveTransformNodeID(glanatCompositeNode.GetID())
    # apply warp to glanat
    glanatCompositeNode.SetAndObserveTransformNodeID(self.parameterNode.GetParameter("warpID"))

  def onModalityPressed(self, item, modality=None):
    if modality is None:
      modality = self.modalityComboBox.itemText(item.row())
    self.parameterNode.SetParameter("modality",modality)
    # find old node and delete
    slicer.mrmlScene.RemoveNode(reducedToolbarLogic().getBackgroundNode())
    # initialize new image and init
    imageNode = ImportSubject.ImportSubjectLogic().importImage(self.parameterNode.GetParameter("subjectPath"), modality)
    self.initializeTransforms(imageNode)
    slicer.util.setSliceViewerLayers(background=imageNode.GetID())
    # load new temaplate image
    slicer.mrmlScene.RemoveNode(reducedToolbarLogic().getForegroundNode())
    # change to t1 in case modality not present
    modality = modality if modality in ['t1','t2','pca','pd'] else 't1'
    templateNode = slicer.util.loadVolume(os.path.join(self.parameterNode.GetParameter("MNIPath"), modality + ".nii"), properties={'show':False})
    templateNode.GetDisplayNode().AutoWindowLevelOff()
    templateNode.GetDisplayNode().SetWindow(100)
    templateNode.GetDisplayNode().SetLevel(70)
    slicer.util.setSliceViewerLayers(foreground=templateNode.GetID())


  def onInteractionModeChanged(self, caller, event):
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    self.windowLevelAction.setChecked(interactionNode.GetCurrentInteractionMode() == 5)

  def sliceIntersectionToggle(self, t):
    compositeCollection = slicer.mrmlScene.GetNodesByClass("vtkMRMLSliceCompositeNode")
    for i in range(compositeCollection.GetNumberOfItems()):
      compositeCollection.GetItemAsObject(i).SetSliceIntersectionVisibility(t)


  def updateToolbarFromMRML(self, caller=None, event=None):
    # subject text
    subjectN = int(self.parameterNode.GetParameter("subjectN"))
    subjectPaths = self.parameterNode.GetParameter("subjectPaths").split(self.parameterNode.GetParameter("separator"))
    self.subjectNameLabel.text = 'Subject: ' + os.path.split(os.path.abspath(self.parameterNode.GetParameter("subjectPath")))[-1]
    self.saveButton.text = 'Finish and Exit' if subjectN == len(subjectPaths)-1 else 'Finish and Next'
    # modality
    self.modalityComboBox.setCurrentText(self.parameterNode.GetParameter("modality"))
    # change resolution
    warpID = self.parameterNode.GetParameter("warpID")
    warpNode = slicer.util.getNode(warpID) if warpID != "" else None
    warpNumberOfComponents = TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode)
    self.resolutionComboBox.enabled = warpNumberOfComponents == 1


  def onSaveButton(self):
    WarpEffect.WarpEffectTool.empty()
    if reducedToolbarLogic().applyChanges():
      # remove nodes
      SmudgeModule.SmudgeModuleLogic().removeRedoNodes()
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(self.parameterNode.GetParameter("glanatCompositeID")))
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(self.parameterNode.GetParameter("warpID")))
      slicer.mrmlScene.RemoveNode(reducedToolbarLogic().getBackgroundNode())

      nextSubjectN = int(self.parameterNode.GetParameter("subjectN"))+1
      subjectPaths = self.parameterNode.GetParameter("subjectPaths").split(self.parameterNode.GetParameter("separator"))
    
      if nextSubjectN < len(subjectPaths):
        self.updateModalities(subjectPaths[nextSubjectN])
        self.parameterNode.SetParameter("subjectN", str(nextSubjectN))
        self.parameterNode.SetParameter("subjectPath", subjectPaths[nextSubjectN])
        reducedToolbarLogic().loadSubjectTransforms()
        self.onModalityPressed([],self.parameterNode.GetParameter("modality"))
        self.updateToolbarFromMRML()
        self.parameterNode.SetParameter("warpModified","0")
      else:
        slicer.util.exit()

    self.parameterNode.SetParameter("subjectChanged","1")

  def updateModalities(self, subjectPath):
    currentModality = self.modalityComboBox.currentText
    subjectModalities = ImportSubject.ImportSubjectLogic().getAvailableModalities(subjectPath)
    if currentModality not in subjectModalities:
      self.parameterNode.SetParameter("modality", subjectModalities[0])
    self.modalityComboBox.clear()
    self.modalityComboBox.addItems(subjectModalities)

  def onResolutionChanged(self, index):
    SmudgeModule.SmudgeModuleLogic().removeRedoNodes()
    newResolution = float(self.resolutionComboBox.itemText(index)[:-2]) # get resolution
    # apply to warp
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID")) # get warp node
    reducedToolbarLogic().resampleTransform(warpNode, newResolution)
    # apply to glanat comp
    glanatCompositeNode = slicer.util.getNode(self.parameterNode.GetParameter("glanatCompositeID")) # get warp node
    reducedToolbarLogic().resampleTransform(glanatCompositeNode, newResolution)
    # save
    self.parameterNode.SetParameter("resolution",str(newResolution))



#
# Logic
#

class reducedToolbarLogic(object):

  def __init__(self):
    self.parameterNode = SmudgeModule.SmudgeModuleLogic().getParameterNode()


  def loadSubjectTransforms(self):
    subjectPath = self.parameterNode.GetParameter("subjectPath")

    # update subject warp fields to new lead dbs specification
    if ImportSubject.ImportSubjectLogic().ish5Transform(subjectPath):
      ImportSubject.ImportSubjectLogic().updateTranform(subjectPath, self.parameterNode.GetParameter("antsApplyTransformsPath"))

    # load glanat composite
    glanatCompositeNode = ImportSubject.ImportSubjectLogic().importTransform(subjectPath, 'glanatComposite.nii.gz')
    self.parameterNode.SetParameter("glanatCompositeID", glanatCompositeNode.GetID())
    # resample
    self.resampleTransform(glanatCompositeNode, float(self.parameterNode.GetParameter("resolution")))

    # create warp
    size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(glanatCompositeNode)
    warpNode = TransformsUtil.TransformsUtilLogic().emptyGridTransform(size,origin,spacing)
    warpNode.CreateDefaultDisplayNodes()
    self.parameterNode.SetParameter("warpID", warpNode.GetID())

  def resampleTransform(self, transformNode, resolution):
    # check resolution
    size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(transformNode)
    if spacing[0] == resolution:
      return
    # aux nodes
    outNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
    size,origin,spacing = TransformsUtil.TransformsUtilLogic().getMNIGrid(resolution)
    referenceVolume = TransformsUtil.TransformsUtilLogic().createEmpyVolume(size,origin,spacing) # aux reference volume with specified resolution 
    # apply
    transformsLogic = slicer.modules.transforms.logic()
    transformsLogic.ConvertToGridTransform(transformNode, referenceVolume, outNode)
    transformNode.SetAndObserveTransformFromParent(outNode.GetTransformFromParent()) # set new transform to warp node
    # remove aux nodes
    slicer.mrmlScene.RemoveNode(outNode)
    slicer.mrmlScene.RemoveNode(referenceVolume)
  
  def applyChanges(self):

    parameterNode = self.parameterNode
    warpNode = slicer.util.getNode(parameterNode.GetParameter("warpID"))
    subjectPath = parameterNode.GetParameter("subjectPath")

    if not bool(int(parameterNode.GetParameter("warpModified"))):
      msgBox = qt.QMessageBox()
      msgBox.setText('No modifications in warp')
      msgBox.setInformativeText('Save subject as approved?')
      msgBox.setStandardButtons(qt.QMessageBox().Save | qt.QMessageBox().Discard | qt.QMessageBox().Cancel)
      ret = msgBox.exec_()
      if ret == qt.QMessageBox().Cancel:
        return False
      if ret == qt.QMessageBox().Save:
        FunctionsUtil.saveApprovedData(subjectPath)
      return True
    
    if TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode) == 1:
      msgBox = qt.QMessageBox()
      msgBox.setText('Only one layer')
      msgBox.setInformativeText('Save changes?')
      msgBox.setStandardButtons(qt.QMessageBox().Save | qt.QMessageBox().Discard | qt.QMessageBox().Cancel)
      ret = msgBox.exec_()
      if ret == qt.QMessageBox().Cancel:
        return False
      elif ret == qt.QMessageBox().Discard:
        return True
    
    # harden changes in glanat composite
    glanatCompositeNode = slicer.util.getNode(self.parameterNode.GetParameter("glanatCompositeID"))
    glanatCompositeNode.HardenTransform()
    TransformsUtil.TransformsUtilLogic().flattenTransform(glanatCompositeNode, True)

    # back to original resolution
    size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(glanatCompositeNode)
    if spacing[0] != 0.5:
      self.resampleTransform(glanatCompositeNode, 0.5)

    # save foreward
    slicer.util.saveNode(glanatCompositeNode, os.path.join(subjectPath,'glanatComposite.nii.gz'))

    # invert transform
    glanatCompositeNode.Inverse()
    # get image to set as reference 
    imageNode = self.getBackgroundNode()
    # get inverse
    outNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
    transformsLogic = slicer.modules.transforms.logic()
    transformsLogic.ConvertToGridTransform(glanatCompositeNode, imageNode, outNode)
    # save inverse
    slicer.util.saveNode(outNode, os.path.join(subjectPath,'glanatInverseComposite.nii.gz'))

    # delete aux node
    slicer.mrmlScene.RemoveNode(outNode)
    

    return True


  def getBackgroundNode(self):
    layoutManager = slicer.app.layoutManager()
    compositeNode = layoutManager.sliceWidget('Red').sliceLogic().GetSliceCompositeNode()
    if compositeNode.GetBackgroundVolumeID():
      return slicer.util.getNode(compositeNode.GetBackgroundVolumeID())
    else:
      return None

  def getForegroundNode(self):
    layoutManager = slicer.app.layoutManager()
    compositeNode = layoutManager.sliceWidget('Red').sliceLogic().GetSliceCompositeNode()
    if compositeNode.GetForegroundVolumeID():
      return slicer.util.getNode(compositeNode.GetForegroundVolumeID())
    else:
      return None