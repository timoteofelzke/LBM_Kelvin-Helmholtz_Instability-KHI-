#include <iostream>
#include <string>
#include <iomanip>

#include "boost/program_options.hpp"
#include "options.h"

using namespace std;
using namespace boost::program_options;

stOptions readCommandLine(int argc, char** argv)
{
	stOptions opt;

	opt.backupInterval = -1;
	opt.filename = "data.ini";
	using namespace boost::program_options;

	options_description desc("Allowed options");

	desc.add_options()
	("help,h", "produce help message")
	("backup,b", value<int>(&opt.backupInterval), "set the interval for saving the backup file sim.bak")
	("endbackup,e","save the backup file at the of the simulation")
	("continue,c", "continue the simulation from using the backup file (sim.bak) ")
	("input,i", value<string>(&opt.filename), "set the INI file for the simulation (default: data.ini)");

	variables_map vm;
	store( parse_command_line(argc, argv, desc), vm );
	notify(vm);

	if ( vm.count("help")  )
	{
		cout << endl << desc << endl;
		exit(0);
	}

	opt.continueFlag 	= vm.count("continue");
	opt.backupFlag	= vm.count("backup");
	opt.endBackupFlag = vm.count("endbackup");

	return opt;
}

void reportOptions(stOptions opt)
{
    cout << left << setw(40) << "Paremeters filename (INI):";
    cout << right << setw(30) << opt.filename;
}

void printInfo()
{
    cout << endl;
    cout << setw(4) << setfill(' ') << " " << setw(62) << setfill('-') << "-" << endl;
    cout << setw(4) << setfill(' ') << " " << " Lattice Boltzmann Algorithm of the Grad Research Laboratory " << endl;
    cout << setw(4) << setfill(' ') << " " << " Two Relaxation Times Model " << endl;
    cout << setw(4) << setfill(' ') << " " << " Version: 0.1 " << endl;
    cout << setw(4) << setfill(' ') << " " << setw(62) << setfill('-') << "-" << endl;
    cout << setfill(' ');
};
