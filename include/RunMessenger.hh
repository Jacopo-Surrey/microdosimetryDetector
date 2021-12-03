#ifndef RunMessenger_h
#define RunMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4String.hh"

#include "RunAction.hh"

class G4UIdirectory;
class G4UIcmdWithAnInteger;

class RunAction;

class RunMessenger: public G4UImessenger
{
public:
  RunMessenger(RunAction* action);
  ~RunMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
  
  // functions to read parameters
  G4int GetRequiredCounts() { return requiredCounts; }
  
private:
  G4UIdirectory *customizeRunDir;      // UI directory
  
  G4UIcmdWithAnInteger *changeRequiredCounts;
	// Set how many counts to collect before stopping the simulation
  
  RunAction* runAction;
  
  // parameters to store
  G4int requiredCounts;
};
#endif

