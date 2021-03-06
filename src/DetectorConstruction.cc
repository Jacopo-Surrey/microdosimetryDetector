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

// Modified by Jacopo Magini: j.magini@surrey.ac.uk

#include "DetectorConstruction.hh"
#include "globals.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
//#include "G4SubtractionSolid.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "SensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4NistManager.hh"

DetectorConstruction::DetectorConstruction(AnalysisManager* analysis_manager, DetectorMessenger* detector_messenger)
{
	analysis = analysis_manager;
	messenger = detector_messenger;
	
	// the following are currently only implemented for MicroDiamond
	detectorType = messenger -> GetTheDetector();
	detectorPositionDepth = messenger -> GetDetectorPositionDepth();
	detectorSizeWidth = messenger -> GetDetectorSizeWidth();
	detectorSizeThickness = messenger -> GetDetectorSizeThickness();
	usingPhantom = messenger -> GetUsingPhantomBool();
	multiSV = messenger -> GetMultiSVBool();
		
	if( multiSV == true )	nOfSV = 4;
	else if( multiSV == false )	nOfSV = 1;
}

DetectorConstruction::~DetectorConstruction(){

}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	// build water phantom (or don't)
	// also sets physical_world and logical_motherVolumeForDetector private members for later
	if( usingPhantom == true ) ConstructWithWaterPhantom();
	
	else if( usingPhantom == false ) ConstructWithoutWaterPhantom();
	
	else
	{
		G4cout << "ERROR: usingPhantom not set to true/false";
		
		return 0;
	}
	
	// build a detector
	// and place it in whichever mother volume you've just build
	if( detectorType == "Diamond" ) ConstructDiamondDetector();
	
	else if( detectorType == "MicroDiamond" ) ConstructMicroDiamondDetector();
	
	else if( detectorType == "Telescope" ) ConstructMicroDiamondDetector();
	
	else if( detectorType == "Silicon" ) ConstructSiliconDetector();
	
	else
	{
		G4cout << "ERROR: " << detectorType << " is not an allowed detector type. ";
		G4cout << "Did you change some code in radioprotection.cc and/or DetectorMessenger.cc ?" << G4endl;
		
		return 0;
	}
	
	return physical_world;
}

void DetectorConstruction::ConstructWithWaterPhantom()	// make changes according to comments    &      check that Region is working       is the size big enough for every detector type??????
{
	//Define Hydrogen 
	G4double A = 1.01 * g/mole;
	G4double Z = 1;
	G4Element* elH = new G4Element ("Hydrogen", "H", Z, A);
	
	//Define Oxygen
	A = 16.0 * g/mole;
	Z = 8;
	G4Element* elO = new G4Element ("Oxygen", "O", Z, A);
	
	//Define Vacuum
	A = 1.01*g/mole;
	Z = 1.;
	G4double vacuumDensity = 1.e-25 *g/cm3;
	G4double pressure = 3.e-18*pascal;
	G4double temperature = 2.73*kelvin;
	G4Material* vacuum = new G4Material("Galactic", Z, A,
			         vacuumDensity,kStateGas,temperature,pressure);

	//define water
	G4Material* water = new G4Material("water", 1*g/cm3, 2);
	water -> AddElement(elH, 2);
	water -> AddElement(elO, 1);
	
	//Define volumes
	// World volume  has size 1cm
	G4double worldx = 0.5 * m;
	G4double worldy = 0.5 * m;
	G4double worldz = 0.5 * m;

	// World volume, containing all geometry
	G4Box* world = new G4Box("world_box", worldx, worldy, worldz);

	G4LogicalVolume* logical_world = new G4LogicalVolume(world, vacuum, "world_log", 0,0,0);

	//set the logical world volume invisible
	logical_world -> SetVisAttributes(G4VisAttributes::GetInvisible());

	physical_world = new G4PVPlacement(0,
								G4ThreeVector(),
								logical_world, 
								"world_phys",
								0, 
								false, 
								0);
	
	//water phantom
	G4double phantom_x = 300.*mm /2.;
	G4double phantom_y = 300.*mm /2.;
	G4double phantom_z = 300.*mm /2.;	// this value is also hardcoded in DetectorMessenger.cc (minus 10.*mm)
	
	G4Box* phantom_box = new G4Box("phantom_box", phantom_x, phantom_y, phantom_z);
	
	G4LogicalVolume* logical_phantom = new G4LogicalVolume(phantom_box, water, "phantom_log", 0,0,0);
	
	G4ThreeVector phantom_position = G4ThreeVector( 0., 0., phantom_z );	//the phantom starts at z=0
	
	 new G4PVPlacement(0, phantom_position, logical_phantom,"phantom_phys",
				logical_world, 
				false, 0, true);
	
	logical_phantom -> SetVisAttributes(G4VisAttributes(G4Colour(0., 0.2, 0.6)));
	
	//smaller volume where I can lower the cuts with G4Region
	G4double highPVol_x = 11.*mm /2.; 	// SET DEPENDING ON SIZE OF DETECTOR
	G4double highPVol_y = 11.*mm /2.;
	G4double highPVol_z = 6.*mm /2.;

	G4Box* highPVol_box = new G4Box("highPVol_box", highPVol_x, highPVol_y, highPVol_z);
 
	G4LogicalVolume* logical_highPVol = new G4LogicalVolume(highPVol_box, water, "highPVol_log",0,0,0);
	
	G4double detectorDepth = detectorPositionDepth;
	
	G4double detectorCentre_z = -phantom_z + detectorDepth;
	
	G4ThreeVector highP_position = G4ThreeVector( 0 , 0 , detectorCentre_z );
 
	new G4PVPlacement(0, highP_position, logical_highPVol,"highPVol_phys",
				logical_phantom, 
				false, 0, true);
 
	logical_highPVol -> SetVisAttributes(G4VisAttributes(G4Colour(0., 0., 1.)));
	
	// private member of DetectorConstruction,
	// needed as mother volume for Construct*Detector()
	logical_motherVolumeForDetector = logical_highPVol;
	
	// high precision region
	G4Region* highPRegion = new G4Region("highPRegion");
	highPRegion -> AddRootLogicalVolume( logical_highPVol );
}

void DetectorConstruction::ConstructWithoutWaterPhantom()
{
	//Define Vacuum
	G4double Z = 1.;
	G4double A = 1.01*g/mole;
	G4double vacuumDensity = 1.e-25 *g/cm3;
	G4double pressure = 3.e-18*pascal;
	G4double temperature = 2.73*kelvin;
	G4Material* vacuum = new G4Material("Galactic", Z, A,
						 vacuumDensity,kStateGas,temperature,pressure);
	
	G4double worldx = 0.5 * m;  //half length!!!!
	G4double worldy = 0.5 * m;
	G4double worldz = 0.5 * m;

	// World volume, containing all geometry
	G4Box* world = new G4Box("world_box", worldx, worldy, worldz);

	G4LogicalVolume* logical_world = new G4LogicalVolume(world, vacuum, "world_log", 0,0,0);

	//set the logical world volume invisible
	logical_world -> SetVisAttributes(G4VisAttributes::GetInvisible());

	physical_world = new G4PVPlacement(0,
								G4ThreeVector(),
								logical_world, 
								"world_phys",
								0, 
								false, 
								0);
	
	logical_motherVolumeForDetector = logical_world;
}

void DetectorConstruction::ConstructDiamondDetector() // change return value  --- COMMENTED OUT!!!
{/*

	//Define each individual element
	//Define Nitrogen
	G4double A = 14.01 * g/mole;
	G4double Z = 7;
	G4Element* elN = new G4Element ("Nitrogen", "N", Z, A);

	//Define Oxygen
	A = 16.0 * g/mole;
	Z = 8;
	G4Element* elO = new G4Element ("Oxygen", "O", Z, A);

	//Define Hydrogen 
	A = 1.01 * g/mole;
	Z = 1;
	G4Element* elH = new G4Element ("Hydrogen", "H", Z, A);

	//Define Boron
	A = 10.8 * g/mole;
	Z = 5;
	G4Element* elB = new G4Element ("Boron", "B", Z, A);

	//Define Carbon
	A = 12.01 * g/mole;
	Z = 6;
	G4Element* elC = new G4Element ("Carbon", "C", Z, A);

	//Define Air   
	G4Material* Air = new G4Material("Air", 1.29*mg/cm3, 2);
	Air -> AddElement(elN, 70*perCent);
	Air -> AddElement(elO, 30*perCent);

	//Define diamond
	A = 12.01 * g/mole;
	Z = 6;
	G4Material* diamond = new G4Material("diamond", Z, A, 3.515*g/cm3);
				
	//Define dopant (boron doped diamond)
	G4Material* dopant = new G4Material("dopant", 3.514*g/cm3, 2);
	dopant -> AddElement(elC, 99.9994*perCent);
	dopant -> AddElement(elB, 0.0006*perCent);

	//Define Aluminium contacts (AlContact)
	A = 26.981 * g/mole;
	Z = 13;
	G4Material* AlContact = new G4Material("AlContact", Z, A, 2.7 *g/cm3);

	//Define Gold contact (AuContact)
	A = 196.97 * g/mole;
	Z = 79;
	G4Material* AuContact = new G4Material("AuContact", Z, A, 19.3 *g/cm3);
	
	//Define PMMA (C502H8)
	// NIST reference 
	G4Material* PMMA = new G4Material("PMMA", 1.19*g/cm3, 3);
	PMMA -> AddElement(elC, 5);
	PMMA -> AddElement(elO, 2);
	PMMA -> AddElement(elH, 8);

	//define water
	G4Material* water = new G4Material("water", 1*g/cm3, 2);
	water -> AddElement(elH, 2);
	water -> AddElement(elO, 1);
		
	//Define Vacuum
	G4double vacuumDensity = 1.e-25 *g/cm3;
	G4double pressure = 3.e-18*pascal;
	G4double temperature = 2.73*kelvin;
	G4Material* vacuum = new G4Material("Galactic", Z=1., A=1.01*g/mole,
						vacuumDensity,kStateGas,temperature,pressure);

	//Define volumes
	// World volume  has size 1cm
	G4double worldx = 0.5 * m;  //half length!!!!
	G4double worldy = 0.5 * m;
	G4double worldz = 0.5 * m;

	// World volume, containing all geometry
	G4Box* world = new G4Box("world_box", worldx, worldy, worldz);

	G4LogicalVolume* logical_world = new G4LogicalVolume(world, vacuum, "world_log", 0,0,0);

	//set the logical world volume invisible
	logical_world -> SetVisAttributes(G4VisAttributes::GetInvisible());

	G4VPhysicalVolume* physical_world = new G4PVPlacement(0,
								G4ThreeVector(),
								logical_world, 
								"world_phys",
								0, 
								false, 
								0);

		
	// Define the geometry of the diamond microdosimeter
	// mother volume of the detector components
	G4double DiaVol_x = 300*micrometer;
	G4double DiaVol_y = 240*micrometer;
	G4double DiaVol_z = 150*micrometer; 

	G4Box* DiaVol_box = new G4Box("DiaVol_box",DiaVol_x,DiaVol_y,DiaVol_z);

	G4LogicalVolume* logical_DiaVol = new G4LogicalVolume(DiaVol_box, diamond, "DiaVol_log", 0,0,0);

	new G4PVPlacement(0, G4ThreeVector(0,0,0), logical_DiaVol,"DiaVol_phys",
				logical_world, 
				false, 0, true);

	//VacBlock for contact placement
	G4double vacblock_x = 300*um;
	G4double vacblock_y = 240*um;
	G4double vacblock_z = 0.25*um; 

	G4Box* vacblock_box = new G4Box("vacblock_box",vacblock_x,vacblock_y,vacblock_z);

	G4LogicalVolume* logical_vacblock = new G4LogicalVolume(vacblock_box, vacuum, "vacblock_log", 0,0,0);

	new G4PVPlacement(0, 
				G4ThreeVector(0,0,DiaVol_z - vacblock_z),
			logical_vacblock,
				"vacblock_phys",
				logical_DiaVol, 
				false, 
				0, true);
	//Bdl in DiaVol
	G4double Bdl_x = 300*micrometer;
	G4double Bdl_y = 240*micrometer;
	G4double Bdl_z = 0.69*micrometer; 
		
	G4Box* Bdl_box = new G4Box("Bdl_box",Bdl_x,Bdl_y,Bdl_z);

	G4LogicalVolume* logical_Bdl = new G4LogicalVolume(Bdl_box, dopant, "Bdl_log", 0,0,0);

	new G4PVPlacement(0, 
			G4ThreeVector(0,0,DiaVol_z - Bdl_z - vacblock_z- vacblock_z),
			logical_Bdl,
			"Bdl_phys",
				logical_DiaVol,   //mother volume 
					false, 
			0, true);

	//Diamond SV
	G4double SV_x = 75*um;
	G4double SV_y = 75*um;
	G4double SV_z = 0.69*um; 

	G4Box* SV_box = new G4Box("SV_box",SV_x,SV_y,SV_z);

	G4LogicalVolume* logical_SV = new G4LogicalVolume(SV_box, diamond, "SV_log", 0,0,0);

	new G4PVPlacement(0, G4ThreeVector(-45*um,105*um,0*um), logical_SV,"SV_phys1",
				logical_Bdl,false, 0, true);

	new G4PVPlacement(0, G4ThreeVector(165*um,105*um,0*um), logical_SV,"SV_phys2",
			logical_Bdl, false, 0, true);

	new G4PVPlacement(0, G4ThreeVector(-45*um,-105*um,0*um),logical_SV,"SV_phys3", 
			logical_Bdl, false, 0, true);

	new G4PVPlacement(0, G4ThreeVector(165*um,-105*um,0*um),logical_SV,"SV_phys4",
			logical_Bdl, false, 0, true);

	//Al strips
	//10 nm thickness 
	G4double AlStrip_x = 240*um;
	G4double AlStrip_y = 240*um;
	G4double AlStrip_z = vacblock_z; 

	G4Box* AlStrip = new G4Box("AlStrip",AlStrip_x,AlStrip_y,AlStrip_z);

	G4LogicalVolume* logical_AlStrip = new G4LogicalVolume(AlStrip, AlContact, "AlStrip_log", 0,0,0);

	new G4PVPlacement(0, G4ThreeVector(60*um,0,0), logical_AlStrip, "AlStrip_phys",
					logical_vacblock, false, 0, true);

	//gold cylinder in vacblock
	G4double innerRadiusOfTheTube1 = 0.*um;
	G4double outerRadiusOfTheTube1 = 45.*um;
	G4double heightOfTheTube1 = 10*nm;
	G4double startAngleOfTheTube1 = 0.*deg;
	G4double spanningAngleOfTheTube1 = 360.*deg;

	G4Tubs* GoldCylinder1 = new G4Tubs("GoldCylinder1", innerRadiusOfTheTube1, 
								outerRadiusOfTheTube1,
								heightOfTheTube1,
								startAngleOfTheTube1, 
								spanningAngleOfTheTube1);
		
	G4LogicalVolume* logical_GoldCylinder1 = new G4LogicalVolume(GoldCylinder1, AuContact, "GoldCylinder1_log", 0,0,0);

	new G4PVPlacement(0,G4ThreeVector(-245*um,0,-vacblock_z + heightOfTheTube1),
						logical_GoldCylinder1,
					"GoldCylinder1_phys",
					logical_vacblock, false, 0, true);

	//gold contacts
	G4double innerRadiusOfTheTube2 = 0.*um;
	G4double outerRadiusOfTheTube2 = 45.*um;
	G4double heightOfTheTube2 = Bdl_z;
	G4double startAngleOfTheTube2 = 0.*deg;
	G4double spanningAngleOfTheTube2 = 360.*deg;

	G4Tubs* GoldCylinder2 = new G4Tubs("GoldCylinder2",
								innerRadiusOfTheTube2, 
								outerRadiusOfTheTube2,
								heightOfTheTube2,
								startAngleOfTheTube2, 
								spanningAngleOfTheTube2);
		
	G4LogicalVolume* logical_GoldCylinder2 = new G4LogicalVolume(GoldCylinder2, AuContact, "GoldCylinder2_log", 0,0,0);

	new G4PVPlacement(0, G4ThreeVector(-245*um,0,0), logical_GoldCylinder2, "GoldCylinder2_phys",
			logical_Bdl, false, 0, true);

	//gold cylinder in DiaVol
	G4double innerRadiusOfTheTube3 = 0.*um;
	G4double outerRadiusOfTheTube3 = 45.*um;
	G4double heightOfTheTube3 = 75.*um -heightOfTheTube2 - heightOfTheTube1 ;
	G4double startAngleOfTheTube3 = 0.*deg;
	G4double spanningAngleOfTheTube3 = 360.*deg;

	G4Tubs* GoldCylinder3 = new G4Tubs("GoldCylinder3",
								innerRadiusOfTheTube3, 
								outerRadiusOfTheTube3,
								heightOfTheTube3,
								startAngleOfTheTube3, 
								spanningAngleOfTheTube3);
		
	G4LogicalVolume* logical_GoldCylinder3 = new G4LogicalVolume(GoldCylinder3, AuContact, "GoldCylinder3_log", 0,0,0);

	new G4PVPlacement(0, G4ThreeVector(-245*um,0,DiaVol_z - heightOfTheTube3 - Bdl_z - Bdl_z - vacblock_z- vacblock_z),
					logical_GoldCylinder3,
					"GoldCylinder3_phys",
					logical_DiaVol, 
					false, 
					0, true);

// Visualisation attributes

        logical_DiaVol -> SetVisAttributes(G4VisAttributes(G4Colour(255,255,255))); //white
	logical_Bdl -> SetVisAttributes(G4VisAttributes(G4Colour(0,255,0)));        //green
	
	G4VisAttributes vis_SV(G4Colour(198, 226, 255));
	vis_SV.SetForceSolid(true);
	logical_SV -> SetVisAttributes(vis_SV);
        logical_vacblock -> SetVisAttributes(G4VisAttributes::GetInvisible());		
	logical_AlStrip -> SetVisAttributes(G4VisAttributes(G4Colour(0, 255, 255)));//cyan
	
	G4VisAttributes vis_GoldCylinder1(G4Colour(255, 255, 0));                    
	vis_GoldCylinder1.SetForceAuxEdgeVisible(true);
	logical_GoldCylinder1 -> SetVisAttributes(vis_GoldCylinder1);
	
	G4VisAttributes vis_GoldCylinder2(G4Colour(255, 255, 0));                    
	vis_GoldCylinder2.SetForceAuxEdgeVisible(true);
	logical_GoldCylinder2 -> SetVisAttributes(vis_GoldCylinder2); 
	
	G4VisAttributes vis_GoldCylinder3(G4Colour(255, 255, 0));                    
	vis_GoldCylinder3.SetForceAuxEdgeVisible(true);
	logical_GoldCylinder3 -> SetVisAttributes(vis_GoldCylinder3); 
        
return physical_world; 
*/
}

void DetectorConstruction::ConstructMicroDiamondDetector()
{
	//Define each individual element
	//Define Nitrogen
	//G4double A = 14.01 * g/mole;
	//G4double Z = 7;
	//G4Element* elN = new G4Element ("Nitrogen", "N", Z, A);

	//Define Boron
	G4double A = 10.8 * g/mole;
	G4double Z = 5;
	G4Element* elB = new G4Element ("Boron", "B", Z, A);

	//Define Carbon
	A = 12.01 * g/mole;
	Z = 6;
	G4Element* elC = new G4Element ("Carbon", "C", Z, A);

	//Define Air   
	//G4Material* Air = new G4Material("Air", 1.29*mg/cm3, 2);
	//Air -> AddElement(elN, 70*perCent);
	//Air -> AddElement(elO, 30*perCent);

	//Define diamond
	A = 12.01 * g/mole;
	Z = 6;
	G4Material* diamond = new G4Material("diamond", Z, A, 3.515*g/cm3);
			
	//Define p-type diamond (boron doped diamond)
	G4Material* p_diamond = new G4Material("p_diamond", 3.514*g/cm3, 2);
	// Boron concentration used is 1e20 cm-3, considering the diamond density and a Boron atomic weight of 10.811u
	p_diamond -> AddElement(elC, 99.94887*perCent);
	p_diamond -> AddElement(elB, 0.05113*perCent);
	
	//Define chromium contact
	G4Material* chromium = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cr");
	
	// sentive volume
	G4double SVside = detectorSizeWidth /2.;
	G4double SVthickness = detectorSizeThickness /2.;
	G4double SVspacing = 200.*um; //edge-edge distance
	
	G4Box* SV_box = new G4Box("SV_box", SVside, SVside, SVthickness);

	G4LogicalVolume* logical_SV = new G4LogicalVolume(SV_box, diamond, "SV_log", 0,0,0);
	
	G4VisAttributes SVcolour(G4Colour(0.5, 0.5, 0.5));
	SVcolour.SetForceSolid(true);
	logical_SV -> SetVisAttributes(SVcolour);

	// chromium front-electrode
	G4double feThickness = 50.*nm /2.; // front-electrode thickness

	G4Box* fe_box = new G4Box("frontElec_box", SVside, SVside, feThickness);

	G4LogicalVolume* logical_fe = new G4LogicalVolume(fe_box, chromium, "frontElec_log", 0,0,0);
	
	G4VisAttributes fe_colour(G4Colour::Brown());
	fe_colour.SetForceSolid(false);
	logical_fe -> SetVisAttributes(fe_colour);
	
	// p-type diamond
	G4double pDthickness = 1.*um /2.; // p-type diamond back-electrode thickness
	
	G4Box* pD_box = new G4Box("pDiam_box", SVside, SVside, pDthickness);
	
	G4LogicalVolume* logical_pD = new G4LogicalVolume(pD_box, p_diamond, "pDiam_box", 0,0,0);
	
	G4VisAttributes pDcolour(G4Colour::Blue());
	pDcolour.SetForceSolid(false);
	
	logical_pD -> SetVisAttributes(pDcolour);

	// put them in place
	G4ThreeVector SVposition = {0., 0., SVthickness};	//position of the first edge (left)
	G4ThreeVector fePosition = {0., 0., -feThickness};
	G4ThreeVector pDposition = {0., 0., 2.*SVthickness + pDthickness};
	
	// if there are multiple sensitive volumes, each needs its own x-coordinate
	G4double* SVposition_x = new G4double[ nOfSV ];
	
	if( nOfSV == 4 )
	{
		//SVposition_x[0] = { -3.*SVside -1.5*SVspacing, -SVside -0.5*SVspacing, +SVside +0.5*SVspacing, +3.*SVside +1.5*SVspacing };
		SVposition_x[0] = -3.*SVside -1.5*SVspacing;
		SVposition_x[1] = -SVside -0.5*SVspacing;
		SVposition_x[2] = +SVside +0.5*SVspacing;
		SVposition_x[3] = +3.*SVside +1.5*SVspacing;
	}
	
	else if( nOfSV == 1 )	SVposition_x[0] = 0.;

	std::ostringstream PVName;

	for( int i=0; i<nOfSV; i++)
	{	
		// sensitive volume
		SVposition[0] = SVposition_x[i];
		PVName << "SV_phys_" << (i+1) ;
		new G4PVPlacement(0, SVposition, logical_SV, PVName.str(),
					logical_motherVolumeForDetector,
					false, 0, true);
		PVName.str("");	//reset the string
		
		// chromium front-electrode
		PVName << "frontElec_phys_" << (i+1);
		fePosition[0] = SVposition[0];
		new G4PVPlacement(0, fePosition, logical_fe, PVName.str(),
					logical_motherVolumeForDetector,
					false, 0, true);
		PVName.str("");
		
		// p-type diamond back-electrode
		PVName << "pD_phys_" << (i+1);
		pDposition[0] = SVposition[0];
		new G4PVPlacement(0, pDposition, logical_pD, PVName.str(),
					logical_motherVolumeForDetector,
					false, 0, true);
		PVName.str("");		
	}

	// HPHT diamond substrate (only one big substrate for simplicity)	// CAMBIA LE MIE DIMENSIONI IN BASE AI SV !!!!
	G4double subs_x = 2.*mm /2.;
	G4double subs_y = 0.5*mm /2.; 
	G4double sub_z = 300.*micrometer /2.; 
	
	G4Box* sub_box = new G4Box("sub_box", subs_x, subs_y, sub_z);
	
	G4LogicalVolume* logical_sub = new G4LogicalVolume(sub_box, diamond, "sub_log", 0,0,0);
	
	G4ThreeVector subPosition = {0,0, 2.*SVthickness + 2.*pDthickness +sub_z};
	
	new G4PVPlacement(0, subPosition, logical_sub, "sub_phys",
				logical_motherVolumeForDetector,
				false, 0, true);

	G4VisAttributes subColour(G4Colour(0.5, 0.5, 0.5));
	subColour.SetForceSolid(false);
	logical_sub -> SetVisAttributes(subColour);
	
	delete SVposition_x;
}

void DetectorConstruction::ConstructTelescopeDetector()		// stub, will be written by Gabriele (?)
{
	// READ ME BEFORE EDITING:
	// use logical_motherVolumeForDetector as if it were the world volume
	// don't return anything at the end calling function takes care of that
	// call the logical volume where you take the microdosimetric spectrum SV_log
	// call the other stage whatever you want
}

void DetectorConstruction::ConstructSiliconDetector()	// change return value   --- COMMENTED OUT!!!
{/*
	//load NIST database
	G4NistManager* nist = G4NistManager::Instance();
	nist->SetVerbose(1);

	//Define each individual element
	//Define Nitrogen
	G4double A = 14.01 * g/mole;
	G4double Z = 7;
	G4Element* elN = new G4Element ("Nitrogen", "N", Z, A);

	//Define Oxygen
	A = 16.0 * g/mole;
	Z = 8;
	G4Element* elO = new G4Element ("Oxygen", "O", Z, A);

	//Define Hydrogen 
	A = 1.01 * g/mole;
	Z = 1;
	G4Element* elH = new G4Element ("Hydrogen", "H", Z, A);

	//Define Carbon
	A = 12.01 * g/mole;
	Z = 6;
	G4Element* elC = new G4Element ("Carbon", "C", Z, A);
	
	//Define Air   
	G4Material* Air = new G4Material("Air", 1.29*mg/cm3, 2);
	Air -> AddElement(elN, 70*perCent);
	Air -> AddElement(elO, 30*perCent);
 
	//Define PMMA (C502H8)
	// NIST reference 
	G4Material* PMMA = new G4Material("PMMA", 1.19*g/cm3, 3);
	PMMA -> AddElement(elC, 5);
	PMMA -> AddElement(elO, 2);
	PMMA -> AddElement(elH, 8);

	//define materials
	G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");
	G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	
	//Define Vacuum
	 G4double vacuumDensity = 1.e-25 *g/cm3;
	 G4double pressure = 3.e-18*pascal;
	 G4double temperature = 2.73*kelvin;
	 G4Material* vacuum = new G4Material("Galactic", Z=1., A=1.01*g/mole,
						 vacuumDensity,kStateGas,temperature,pressure);

	 //Define volumes
	 // World volume  has size 1cm
	 G4double worldx = 0.5 * m;  //half length!!!!
	 G4double worldy = 0.5 * m;
	 G4double worldz = 0.5 * m;

	 // World volume, containing all geometry
	 G4Box* world = new G4Box("world_box", worldx, worldy, worldz);

	 G4LogicalVolume* logical_world = new G4LogicalVolume(world, vacuum, "world_log", 0,0,0);

	 //set the logical world volume invisible
	 logical_world -> SetVisAttributes(G4VisAttributes::GetInvisible());

	 G4VPhysicalVolume* physical_world = new G4PVPlacement(0,
								   G4ThreeVector(),
								   logical_world, 
								"world_phys",
								0, 
								false, 
								0);
		
	// I need to set the size of the SV now, because some other parameters depend on it
	G4double SV_thick = 25.*um /2.; 
	
	// PMMA
	G4double PMMA_x = 200.*um /2.;
	G4double PMMA_y = 200.*um /2.;
	G4double PMMA_z = SV_thick;
	
	G4Box* PMMA_box = new G4Box("PMMA_box", PMMA_x, PMMA_y, PMMA_z);
	
	G4LogicalVolume* logical_PMMA = new G4LogicalVolume(PMMA_box, PMMA, "PMMA_log", 0,0,0);
	
	new G4PVPlacement(0, G4ThreeVector(), logical_PMMA, "PMMA_phys",
					logical_world,
					false, 0, true);
	
	logical_PMMA -> SetVisAttributes(G4VisAttributes(G4Colour(0., 1., 0.)));
	
	// sensitive volumes
	G4double SV_radius = SV_thick *2;	// full length
	
	G4Tubs* SV_cyl = new G4Tubs("SV_cyl", 0., SV_radius, SV_thick, 0.*deg, 360.*deg);
	//G4RotationMatrix* cylRot = new G4RotationMatrix;
	//cylRot->rotateY(M_PI/2.*rad);
		
	G4LogicalVolume* logical_SV = new G4LogicalVolume(SV_cyl, silicon, "SV_log", 0,0,0);
	
	G4VisAttributes SVcolour(G4Colour(0.5, 0.5, 0.5));
	SVcolour.SetForceSolid(true);
	logical_SV -> SetVisAttributes(SVcolour);
	
	G4ThreeVector SVposition;	//if(volumeName != "SV_phys1")
		
	G4double SVspacing = 10.*um;	// distance between the edges of two SV
	
	SVposition = { +SVspacing/2. +SV_radius, +SVspacing/2. +SV_radius, 0. };
	new G4PVPlacement(0, SVposition, logical_SV, "SV_phys1",
						logical_PMMA,
						false, 0, true);
	//new G4PVPlacement(cylRot, SVposition, logical_SV, "SV_phys1",
	//				logical_PMMA,
	//				false, 0, true);

	SVposition = { -SVspacing/2. -SV_radius, +SVspacing/2. +SV_radius, 0. };
	new G4PVPlacement(0, SVposition, logical_SV, "SV_phys2",
						logical_PMMA,
						false, 0, true);
	
	SVposition = { -SVspacing/2. -SV_radius, -SVspacing/2. -SV_radius, 0. };
	new G4PVPlacement(0, SVposition, logical_SV, "SV_phys3",
						logical_PMMA,
						false, 0, true);
						
	SVposition = { +SVspacing/2. +SV_radius, -SVspacing/2. -SV_radius, 0. };
	new G4PVPlacement(0, SVposition, logical_SV, "SV_phys4",
						logical_PMMA,
						false, 0, true);
	
	// Si02 layer
	G4double oxyde_x = PMMA_x;
	G4double oxyde_y = PMMA_y;
	G4double oxyde_z = 1.*um /2.;
	
	G4Box* oxyde_box = new G4Box("oxyde_box", oxyde_x, oxyde_y, oxyde_z);
	
	G4LogicalVolume* logical_oxyde = new G4LogicalVolume(oxyde_box, SiO2, "oxyde_log", 0,0,0);
	
	G4ThreeVector oxyde_position = G4ThreeVector( 0, 0, PMMA_z + oxyde_z );
	new G4PVPlacement(0, oxyde_position, logical_oxyde, "oxyde_phys",
					logical_world,
					false, 0, true);
	
	logical_oxyde -> SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 0.6)));

	// !!! ADD THE G4Region AGAIN, and uncomment the lines in PhysicsList.cc !!!
	
	return physical_world; 
*/
}

void DetectorConstruction::ConstructSDandField()
{
   SensitiveDetector* SD = new SensitiveDetector("SD", "DetectorHitsCollection", analysis);
   G4SDManager::GetSDMpointer()->AddNewDetector(SD);
   SetSensitiveDetector("SV_log", SD);


}
