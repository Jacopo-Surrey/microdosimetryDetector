#include "RunMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

RunMessenger::RunMessenger(RunAction* action)
{
	customizeRunDir = new G4UIdirectory("/run/custom/");
	customizeRunDir -> SetGuidance("Custom run parameters");

	changeRequiredCounts = new G4UIcmdWithAnInteger("/run/custom/setRequiredCounts", this);
	changeRequiredCounts -> SetGuidance("Set the number of events on any SV after which the simulation is stopped");
	changeRequiredCounts -> SetParameterName("Counts", false);
	changeRequiredCounts -> SetRange("Counts >= 100");
	changeRequiredCounts -> AvailableForStates(G4State_PreInit, G4State_Idle);
	
	runAction = action;
}

RunMessenger::~RunMessenger()
{
	delete changeRequiredCounts;

	delete customizeRunDir;
}

void RunMessenger::SetNewValue(G4UIcommand* command, G4String commandContent)
{

	if( command == changeRequiredCounts )
	{
		requiredCounts = G4UIcmdWithAnInteger::GetNewIntValue(commandContent);
		
		runAction -> setHitsRequired(requiredCounts);
		
		G4cout << "Required hits in the SV changed to " << commandContent << G4endl;
	}
}
