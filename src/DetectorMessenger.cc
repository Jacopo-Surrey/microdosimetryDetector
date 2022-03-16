//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Authors: Susanna Guatelli and Francesco Romano
// susanna@uow.edu.au, francesco.romano@ct.infn.it
//
// Created by Jacopo Magini: j.magini@surrey.ac.uk

#include "DetectorMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "G4UIcmdWithoutParameter.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"

DetectorMessenger::DetectorMessenger(AnalysisManager* analysis_manager)
{
    changeTheGeometryDir = new G4UIdirectory("/geometrySetup/");
    changeTheGeometryDir -> SetGuidance("Geometry setup");
    
    changeTheDetectorCmd = new G4UIcmdWithAString("/geometrySetup/selectDetector", this);
    changeTheDetectorCmd -> SetGuidance("Select the type of detector you wish to use");
    changeTheDetectorCmd -> SetParameterName("Material", false);
	changeTheDetectorCmd -> AvailableForStates(G4State_PreInit);
	
	changeDetectorPositionDir = new G4UIdirectory("/geometrySetup/detectorPosition/");
	changeDetectorPositionDir -> SetGuidance("Modify the placement of the detector");
	
	changeDetectorPositionDepthCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorPosition/setDepth", this);
	changeDetectorPositionDepthCmd -> SetGuidance("Set the detector depth inside the water phantom");
	changeDetectorPositionDepthCmd -> SetParameterName("Depth", false);
	changeDetectorPositionDepthCmd -> SetRange("Depth >= 10. && Depth <= 1000.");
	changeDetectorPositionDepthCmd -> SetUnitCategory("Length");
	changeDetectorPositionDepthCmd -> SetDefaultUnit("mm");
	changeDetectorPositionDepthCmd -> AvailableForStates(G4State_PreInit);
	
	changeDetectorDimensionDir = new G4UIdirectory("/geometrySetup/detectorDimension/");
	changeDetectorDimensionDir -> SetGuidance("Modify the dimensions of the detector");
	
	changeDetectorSizeWidthCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorDimension/setWidth", this);
	changeDetectorSizeWidthCmd -> SetGuidance("Set the width of the detector (DE stage in case of telescope detector)");
	changeDetectorSizeWidthCmd -> SetParameterName("Width", false);
	changeDetectorSizeWidthCmd -> SetRange("Width >= 0.1 && Width <= 1000.");
	changeDetectorSizeWidthCmd -> SetUnitCategory("Length");
	changeDetectorSizeWidthCmd -> SetDefaultUnit("um");
	changeDetectorSizeWidthCmd -> AvailableForStates(G4State_PreInit);
	
	changeDetectorSizeThicknessCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorDimension/setThickness", this);
	changeDetectorSizeThicknessCmd -> SetGuidance("Set the thickness of the detector");
	changeDetectorSizeThicknessCmd -> SetParameterName("Thickness", false);
	changeDetectorSizeThicknessCmd -> SetRange("Thickness >= 0.1 && Thickness <= 1000.");
	changeDetectorSizeThicknessCmd -> SetUnitCategory("Length");
	changeDetectorSizeThicknessCmd -> SetDefaultUnit("um");
	changeDetectorSizeThicknessCmd -> AvailableForStates(G4State_PreInit);
	
	changeDetectorSecondStageDir = new G4UIdirectory("/geometrySetup/detectorDimension/secondStage/");
	changeDetectorSecondStageDir -> SetGuidance("Modify the dimensions of the second stage (telescope only)");
	
	changeSecondStageSizeDimCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorDimension/secondStage/setWidth", this);
	changeSecondStageSizeDimCmd -> SetGuidance("Set the dimension (diameter) of the second-stage (telescope detector only)");
	changeSecondStageSizeDimCmd -> SetParameterName("Diameter", false);
	changeSecondStageSizeDimCmd -> SetRange("Diameter >= 0.1 && Diameter <= 1000.");
	changeSecondStageSizeDimCmd -> SetUnitCategory("Length");
	changeSecondStageSizeDimCmd -> SetDefaultUnit("um");
	changeSecondStageSizeDimCmd -> AvailableForStates(G4State_PreInit);

	changeSecondStageThicknessCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorDimension/secondStage/setThickness", this);
	changeSecondStageThicknessCmd -> SetGuidance("Set the thickness of the second-stage (telescope detector only)");
	changeSecondStageThicknessCmd -> SetParameterName("Thickness", false);
	changeSecondStageThicknessCmd -> SetRange("Thickness >= 10. && Thickness <= 1000.");
	changeSecondStageThicknessCmd -> SetUnitCategory("Length");
	changeSecondStageThicknessCmd -> SetDefaultUnit("um");
	changeSecondStageThicknessCmd -> AvailableForStates(G4State_PreInit);
	
	changeWaterPixelSettingsDir = new G4UIdirectory("/geometrySetup/setupWaterPixel/");
	changeWaterPixelSettingsDir -> SetGuidance("Change settings specific to the water pixel");
	
	selectPixelKinScoringTypeCmd = new G4UIcmdWithAnInteger("/geometrySetup/setupWaterPixel/selectKinScoring", this);
	selectPixelKinScoringTypeCmd -> SetGuidance("Set how the kinetic energy is scored:\n\
													0 : no scoring\n\
													1 : scored at the beginning (lowest depth) of the pixel\n\
													2 : scored in the middle (just before the SV) of the pixel\n\
													3 : UNIMPLEMENTED\n\
													4 : beginning + UNIMPLEMENTED\n\
													5 : middle + UNIMPLEMENTED");
	selectPixelKinScoringTypeCmd -> SetParameterName("Method", false);
	selectPixelKinScoringTypeCmd -> SetRange("Method >= 0 && Method <= 5");
	selectPixelKinScoringTypeCmd -> AvailableForStates(G4State_PreInit);
	
	enableWaterPhantomCmd = new G4UIcmdWithABool("/geometrySetup/enableWaterPhantom", this);
	enableWaterPhantomCmd -> SetGuidance("If true, the detector is placed inside a water phantom");
	enableWaterPhantomCmd -> SetParameterName("Phantom", false);
	enableWaterPhantomCmd -> AvailableForStates(G4State_PreInit);
	
	useMultipleSVCmd = new G4UIcmdWithABool("/geometrySetup/useMultipleSV", this);
	useMultipleSVCmd -> SetGuidance("Set to true for multiple SV, false for a single one");
	useMultipleSVCmd -> SetParameterName("MultipleSV", false);
	useMultipleSVCmd -> AvailableForStates(G4State_PreInit);
	
	changeMultiSVSetupDir = new G4UIdirectory("/geometrySetup/multipleSVsetup/");
	changeMultiSVSetupDir -> SetGuidance("Choose how multiple SV are placed");
	
	changeMaximumBreadthForMultiSVCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/multipleSVsetup/setBreadth", this);
	changeMaximumBreadthForMultiSVCmd -> SetGuidance("Set the maximum breadth transversal to the phantom along which multiple SV are placed");
	changeMaximumBreadthForMultiSVCmd -> SetParameterName("Breadth", false);
	changeMaximumBreadthForMultiSVCmd -> SetRange("Breadth >= 10 && Breadth <= 100000.");
	changeMaximumBreadthForMultiSVCmd -> SetUnitCategory("Length");
	changeMaximumBreadthForMultiSVCmd -> SetDefaultUnit("um");
	changeMaximumBreadthForMultiSVCmd -> AvailableForStates(G4State_PreInit);
	
	changeSpacingBetweenSV = new G4UIcmdWithADoubleAndUnit("/geometrySetup/multipleSVsetup/setSpacing", this);
	changeSpacingBetweenSV -> SetGuidance("Set the space between the edges of two SV");
	changeSpacingBetweenSV -> SetParameterName("Space", false);
	changeSpacingBetweenSV -> SetRange("Space >= 1 && Space <= 1000000.");
	changeSpacingBetweenSV -> SetUnitCategory("Length");
	changeSpacingBetweenSV -> SetDefaultUnit("um");
	changeSpacingBetweenSV -> AvailableForStates(G4State_PreInit);
	
	setCutsDir = new G4UIdirectory("/cuts/custom/");
	setCutsDir -> SetGuidance("Set cuts for this specific setup");
	
	changeCutForRegionCmd = new G4UIcmdWithADoubleAndUnit("/cuts/custom/setCutsAroundSV", this);
	changeCutForRegionCmd -> SetGuidance("Set the cuts in the region around the SV");
	changeCutForRegionCmd -> SetParameterName("Cut_SV", false);
	changeCutForRegionCmd -> SetRange("Cut_SV >= 0.05 && Cut_SV <= 50.");
	changeCutForRegionCmd -> SetUnitCategory("Length");
	changeCutForRegionCmd -> SetDefaultUnit("um");
	changeCutForRegionCmd -> AvailableForStates(G4State_PreInit);
	
	changeCutForWorldCmd = new G4UIcmdWithADoubleAndUnit("/cuts/custom/setCutsElsewhere", this);
	changeCutForWorldCmd -> SetGuidance("When using a water phantom, set the cuts farther away from the SV");
	changeCutForWorldCmd -> SetParameterName("Cut_ext", false);
	changeCutForWorldCmd -> SetRange("Cut_ext >= 0.1 && Cut_ext <= 1000.");
	changeCutForWorldCmd -> SetUnitCategory("Length");
	changeCutForWorldCmd -> SetDefaultUnit("um");
	changeCutForWorldCmd -> AvailableForStates(G4State_PreInit);
	
	changeRegionWidthCmd = new G4UIcmdWithADoubleAndUnit("/cuts/custom/setWidthAroundSV", this);
	changeRegionWidthCmd -> SetGuidance("When using a water phantom, set the width of the region around the SV where setCutsAroundSV is applied");
	changeRegionWidthCmd -> SetParameterName("Region_width", false);
	changeRegionWidthCmd -> SetRange("Region_width >= 0.5 && Region_width <= 5000.");
	changeRegionWidthCmd -> SetUnitCategory("Length");
	changeRegionWidthCmd -> SetDefaultUnit("um");
	changeRegionWidthCmd -> AvailableForStates(G4State_PreInit);
	
	applyChangesToGeometryCmd = new G4UIcmdWithoutParameter("/geometrySetup/applyChanges",this);
    applyChangesToGeometryCmd -> SetGuidance("Apply selected changes to the geometry");
	applyChangesToGeometryCmd -> AvailableForStates(G4State_PreInit);
	
	analysis = analysis_manager;
	
	// default values
	detectorType = "MicroDiamond";
	detectorDepth = 10*mm;
	detectorWidth = 100.*um;
	detectorThickness = 8.*um;
	secondStageDim = 500.*um;
	secondStageThickness = 500.*um;
	pixelKinScoring = 0;
	usingPhantom = false;	// used to crash the program if set to false under some conditions (why?)
	multiSV = false;
	multiSVbreadth = 5000*um;
	svSpacing = detectorThickness*2;
	
	cutForWorld = 100.*um;
	cutForRegion = 0.1*um;
	regionWidth = 1*mm;
	
	pendingChanges = false;
}

DetectorMessenger::~DetectorMessenger()
{
    delete changeTheDetectorCmd;
	delete changeDetectorPositionDepthCmd;
	delete changeDetectorSizeWidthCmd;
	delete changeSecondStageSizeDimCmd;
	delete changeSecondStageThicknessCmd;
	delete changeDetectorSizeThicknessCmd;
	delete selectPixelKinScoringTypeCmd;
	delete enableWaterPhantomCmd;
	delete useMultipleSVCmd;
	delete changeMaximumBreadthForMultiSVCmd;
	delete changeSpacingBetweenSV;
	delete changeCutForRegionCmd;
	delete changeCutForWorldCmd;
	delete changeRegionWidthCmd;
	
	delete applyChangesToGeometryCmd;
	
	delete changeDetectorPositionDir;
	delete changeDetectorDimensionDir;
	delete changeTheGeometryDir;
	delete changeDetectorSecondStageDir;
	delete changeWaterPixelSettingsDir;
	delete changeMultiSVSetupDir;
	delete setCutsDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String commandContent)
{

	if( command == changeTheDetectorCmd )
	{
		if( commandContent == "Diamond" || commandContent == "MicroDiamond" || commandContent == "WPMicroDiamond" || commandContent == "Telescope" || commandContent == "Silicon" || commandContent == "WaterPixel" )
		{
			detectorType = commandContent;
			
			G4cout << "Detector type changed to " << commandContent << G4endl;
			//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
			
			pendingChanges = true;
		}
		
		else
		{
			G4cout <<"Unknown detector type: " << commandContent << ". Geometry not changed" << G4endl;
		}
	}
	
	else if( command == changeDetectorPositionDepthCmd )
	{
		detectorDepth = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Detector depth in water changed to " << commandContent << G4endl;
		if( usingPhantom == false )	G4cout << "However the water phantom is not enabled. Enable it with '/geometrySetup/enableWaterPhantom true' or this value will be ignored" << G4endl;
		//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
		
		pendingChanges = true;
	}
	
	else if( command == changeDetectorSizeWidthCmd )
	{
		detectorWidth = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Detector width changed to " << commandContent << G4endl;
		//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
		
		pendingChanges = true;
	}

	else if( command == changeSecondStageSizeDimCmd )
	{
		secondStageDim = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Second-stage (telescope detector) diameter changed to " << commandContent << G4endl;
		//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
		if( detectorType != "Telescope" )
		{
			G4cout << "WARNING: the detector type is currently set to " << detectorType << G4endl;
			G4cout << "Unless this is changed to Telescope, the last command will be ignored" << G4endl;
		}
		
		pendingChanges = true;
	}

	else if( command == changeSecondStageThicknessCmd )
	{
		secondStageThickness = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Second-stage (telescope detector) thickness changed to " << commandContent << G4endl;
		//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
		if( detectorType != "Telescope" )
		{
			G4cout << "WARNING: the detector type is currently set to " << detectorType << G4endl;
			G4cout << "Unless this is changed to Telescope, the last command will be ignored" << G4endl;
		}
		
		pendingChanges = true;
	}
	
	else if( command == changeDetectorSizeThicknessCmd )
	{
		detectorThickness = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Detector thickness changed to " << commandContent << G4endl;
		//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
		
		pendingChanges = true;
	}
	
	else if( command == selectPixelKinScoringTypeCmd )
	{
		pixelKinScoring = G4UIcmdWithAnInteger::GetNewIntValue(commandContent);
		
		G4cout << "Kinetic energy scoring set to ";
		switch(pixelKinScoring)
		{
			case 0 :
				G4cout << "no scoring";
				break;
			case 1 :
				G4cout << "beginning";
				break;
			case 2 :
				G4cout << "middle";
				break;
			case 3 :
				G4cout << "UNIMPLEMENTED";
				break;
			case 4 :
				G4cout << "beginning + UNIMPLEMENTED";
				break;
			case 5 :
				G4cout << "middle + UNIMPLEMENTED";
				break;
		}
		G4cout << G4endl;
		
		if( detectorType != "WaterPixel" )
		{
			G4cout << "WARNING: the detector type is currently set to " << detectorType << G4endl;
			G4cout << "Unless this is changed to WaterPixel, the last command will be ignored" << G4endl;
		}
		
		pendingChanges = true;
	}
	
	else if( command == enableWaterPhantomCmd )
	{
		usingPhantom = G4UIcmdWithABool::GetNewBoolValue(commandContent);
		
		if( usingPhantom == true )	G4cout << "Water phantom enabled" << G4endl;
		else if( usingPhantom == false )	G4cout << "Water phantom disabled" << G4endl;
		else	G4cout << "Error: " << commandContent << "is not a bool" << G4endl;
		
		//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
		
		pendingChanges = true;
	}
	
	else if( command == useMultipleSVCmd )
	{
		multiSV = G4UIcmdWithABool::GetNewBoolValue(commandContent);
		
		if( multiSV == true )	G4cout << "Now using multiple SV" << G4endl;
		else if( multiSV == false )	G4cout << "Now using a single SV" << G4endl;
		else	G4cout << "Error: " << commandContent << "is not a bool" << G4endl;
		
		//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
		
		pendingChanges = true;
	}
	
	else if( command == changeMaximumBreadthForMultiSVCmd )
	{
		multiSVbreadth = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Breadth along which SV are placed changed to " << commandContent << G4endl;
		//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
		if( multiSV == false )
		{
			G4cout << "WARNING: currently only a single SV is enabled" << G4endl;
			G4cout << "Switch to multiple SV via /geometrySetup/useMultipleSV, otherwise the last command will be ignored" << G4endl;
		}
		
		pendingChanges = true;
	}
	
	else if( command == changeSpacingBetweenSV )
	{
		svSpacing = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Space between SV changed to " << commandContent << G4endl;
		//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
		if( multiSV == false )
		{
			G4cout << "WARNING: currently only a single SV is enabled" << G4endl;
			G4cout << "Switch to multiple SV via /geometrySetup/useMultipleSV, otherwise the last command will be ignored" << G4endl;
		}
		
		pendingChanges = true;
	}
	
	else if( command == changeCutForRegionCmd )
	{
		cutForRegion = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Cuts around the SV changed to " << commandContent << G4endl;
		//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
		
		pendingChanges = true;
	}
	
	else if( command == changeCutForWorldCmd )
	{
		cutForWorld = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Cuts away from the SV changed to " << commandContent << G4endl;
		//G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
		if( usingPhantom == false )
		{
			G4cout << "WARNING: this value only takes effect when using a water phantom" << G4endl;
			G4cout << "Enable it via /geometrySetup/enableWaterPhantom, otherwise the last command will be ignored" << G4endl;
		}
		
		pendingChanges = true;
	}
	
	else if( command == changeRegionWidthCmd )
	{
		regionWidth = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Width of the region around the SV with different cuts changed to " << commandContent << G4endl;
		
		if( usingPhantom == false )
		{
			G4cout << "WARNING: this value only takes effect when using a water phantom" << G4endl;
			G4cout << "Enable it via /geometrySetup/enableWaterPhantom, otherwise the last command will define the size of the world volume" << G4endl;
				// the last effect is not ideal and will be changed later
		}
		
		pendingChanges = true;
	}
	
	else if( command == applyChangesToGeometryCmd )
	{
		G4RunManager* runManager = G4RunManager::GetRunManager();
		
		if( pendingChanges == true )
		{
			// if geometryHasChanged ?
			DetectorConstruction* newDetector = new DetectorConstruction(analysis, this);
			runManager -> SetUserInitialization(newDetector);
			runManager -> GeometryHasBeenModified();
			
			// (nested?) if regionsHaveChanged ?
			G4VUserPhysicsList* newPhysics = new PhysicsList(this);
			runManager -> SetUserInitialization(newPhysics);
			runManager -> PhysicsHasBeenModified();
			
			G4cout << "Changes to geometry have been applied" << G4endl;
			
			pendingChanges = false;
		}
		
		else if( pendingChanges == false )
		{
			G4cout << "No pending changes. The last command has been ignored" << G4endl;
		}
	}
}
