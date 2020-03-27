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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//

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

DetectorConstruction::DetectorConstruction(AnalysisManager* analysis_manager)
{
analysis = analysis_manager;
}

DetectorConstruction::~DetectorConstruction(){

}

#ifndef USING_SILICON
//#ifndef USING_SILICON	// by default, build the diamond detector
G4VPhysicalVolume* DetectorConstruction::Construct()
{

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
 //A = 26.981 * g/mole;
 //Z = 13;
 //G4Material* AlContact = new G4Material("AlContact", Z, A, 2.7 *g/cm3);

 //Define Gold contact (AuContact)
 //A = 196.97 * g/mole;
 //Z = 79;
 //G4Material* AuContact = new G4Material("AuContact", Z, A, 19.3 *g/cm3);
 
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
 G4double worldx = 1*m /2.;  //half length!!!!
 G4double worldy = 1*m /2.;
 G4double worldz = 1*m /2.;

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

 // I ignore the Aluminium walls. The inside of the chamber is filled with air
 G4double chamVol_x = 250*millimeter /2.;
 G4double chamVol_y = 150*millimeter /2.;
 G4double chamVol_z = 116*millimeter /2.;

 G4Box* chamVol_box = new G4Box("chamVol_box", chamVol_x, chamVol_y, chamVol_z);
 
 G4LogicalVolume* logical_chamVol = new G4LogicalVolume(chamVol_box, Air, "chamVol_log",0,0,0);
 
 new G4PVPlacement(0, G4ThreeVector(0,0,0), logical_chamVol,"chamVol_phys",
			  logical_world, 
			  false, 0, true);
 
 //logical_chamVol -> SetVisAttributes(G4VisAttributes(G4Colour(255,255,255))); //white
 logical_chamVol -> SetVisAttributes(G4VisAttributes::GetInvisible());
 
  // 4 later: smaller air volume where I can lower the cuts with G4Region
 G4double highPVol_x = 8*millimeter /2.; 
 G4double highPVol_y = 8*millimeter /2.;
 G4double highPVol_z = 6*millimeter /2.;

 G4Box* highPVol_box = new G4Box("highPVol_box", highPVol_x, highPVol_y, highPVol_z);
 
 G4LogicalVolume* logical_highPVol = new G4LogicalVolume(highPVol_box, Air, "highPVol_log",0,0,0);
 
 G4double sourceDetectorAxis_x = chamVol_x - 70*millimeter;
 
 G4ThreeVector highP_position = G4ThreeVector( sourceDetectorAxis_x , 0 , - chamVol_z + 20*millimeter );
 
 new G4PVPlacement(0, highP_position, logical_highPVol,"highPVol_phys",
			  logical_chamVol, 
			  false, 0, true);
 
 logical_highPVol -> SetVisAttributes(G4VisAttributes(G4Colour(255,255,255))); //white
 
 // Sensitive volumes
 G4double SVspacing_cc = 500.*micrometer; //centre-centre distance
 G4double SVthickness = 8.*micrometer /2.;
 //G4double SVthickness = 1.*micrometer /2.;
 
 // 4 sensitive volumes
 G4double SVside[4] = { 50.*micrometer /2., 300.*micrometer /2., 100.*micrometer /2., 200.*micrometer /2.};
 
 G4Box* SV_box[4];
 G4LogicalVolume* logical_SV[4];
 G4ThreeVector SVposition;
 
 std::ostringstream name;
 
 G4VisAttributes SVcolour(G4Colour(198, 226, 255));
 SVcolour.SetForceSolid(true);
 
 for( int i=0; i<4; i++)
 {
		name << "SV_box_" << i;
		SV_box[i] = new G4Box(name.str(), SVside[i], SVside[i], SVthickness);
		name.str(""); //clears the string 
		
		name << "SV_log_" << i;
		logical_SV[i] = new G4LogicalVolume(SV_box[i], diamond, name.str(), 0,0,0);
		name.str("");
		
		name << "SV_phys_" << i;
		SVposition = { -1.5 * SVspacing_cc + i * SVspacing_cc, 0, SVthickness };
		new G4PVPlacement(0, SVposition, logical_SV[i], name.str(),
					logical_highPVol,
					false, 0, true);
		name.str("");
		
		logical_SV[i] -> SetVisAttributes(SVcolour);
 }
 
 //substrate
 //SCRIVI LE DIMENSIONI DEL SUBSTRATO
        
return physical_world; 

}
#else	// if the flag is on, build the reference silicon instead
G4VPhysicalVolume* DetectorConstruction::Construct()	//UNIMPLEMENTED, will give an error
{
	return 0;
}
#endif

void DetectorConstruction::ConstructSDandField()
{
	G4int sensitiveVolumeToOutput = 1;
	
	std::ostringstream outName; outName << "SV_log_" << sensitiveVolumeToOutput;

	SensitiveDetector* SD = new SensitiveDetector("SD", "DetectorHitsCollection", analysis);
	G4SDManager::GetSDMpointer()->AddNewDetector(SD);
	SetSensitiveDetector(outName.str(), SD);
}
