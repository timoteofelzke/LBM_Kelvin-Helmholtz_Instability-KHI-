/**
 * @file options.h
 * @author Diogo Nardellli Siebert
 * @date 20/07/15
 * @brief Header file for functions to read and print (report) options from the commandline
 */
 
#ifndef OPTIONS_H_INCLUDED__
#define OPTIONS_H_INCLUDED__

#include <string>

/**
 * @struct stOptions
 * @brief Struct to store options read from command line 
 * 
 * This struct is used to store all the parameters that can be read
 * from command line. These parameters are limited to input,output and backup
 * and not related to physical parameters of the simulation (these kind of 
 * parameters will be read from the ini file) b
 */

struct stOptions
{
    std::string filename;  /*!< The Name of INI file  */
    bool continueFlag;     /*!< A bool to indicate if the current simulation will continue from previous backup */
    bool backupFlag;	   /*!< A bool to indicate if a backup file is to be generated  */
    bool endBackupFlag;    /*!< A bool to indicate if a backup must be generated and end of the simulation */
    int  backupInterval;   /*!< The interval of step in which the backup must be written to the disk. */
};

/**
 * @brief Parse command line arguments.
 * @param argc Number of arguments passed trought the command line (argc main function argument)
 * @param argv array containings the chars passed throught the commadn line (argv main function argument)
 * 
 * Function to parse command line arguments (flags) and store them into a stOptions struct
 **/
 
stOptions readCommandLine(int argc, char** argv);

/**
 * @brief Report Parsed Command Lines Arguments
 * @param opt
 * 
 * Print in to the the screen (cout) a formated report with the arguments stored in a stOptions struct. 
 * This function will be executed in the begin of the program in order for the user to check if the options
 * were correcly parsed by the algorithm
 */
 
void reportOptions(stOptions opt);

/**
 * @brief Print the program name and version in the screen
 * 
 * A simple function to print (cout) in the screen the name of the program, the version and the name
 * of the research group.
 */
 
void printInfo();

#endif
