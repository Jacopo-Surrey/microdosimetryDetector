#ifndef PHYSICS_LIST_HH
#define PHYSICS_LIST_HH

#include <G4VModularPhysicsList.hh>
#include "DetectorMessenger.hh"

class PhysicsList : public G4VModularPhysicsList
{
public:
	PhysicsList(DetectorMessenger* detectorMessenger);
	~PhysicsList(){;};

	//Set cuts per region!
	void 	SetCuts ();

private:
	G4double fCutForGamma;
	G4double fCutForElectron;
	G4double fCutForPositron;
	G4double fCutForProton;
	
	G4bool usingRegion;
};

#endif
