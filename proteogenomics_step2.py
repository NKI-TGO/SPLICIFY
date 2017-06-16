#########################################################################
#                                                                       #
# Proteogenomic pipeline for differential splice variant identification #
#                                                                       #
#                              STEP 2                                   #
#                            27/01/2017                                 #
#                                                                       #
#                             MA Komor                                  #
#                                                                       #
#########################################################################

# python 2.*

__author__ = "Malgorzata Anna Komor"
__email__ = "g.komor@nki.nl"

# import 

import ConfigParser
import sys
import subprocess
import shlex
import os.path
import argparse


# read config file
def readInput(arg='src/config.ini'):
    print "#--------------------------#\n|Reading configuration file|\n#--------------------------#"
    config = ConfigParser.ConfigParser()
    config.optionxform = str

    config.read(arg)

    # check if all sections present
    for section in ["MaxQuant output", "Extract"]:
        if section not in config.sections():
            raise Exception("Section %s missing from config file" % section)

    # check if splice variant database file exist
    sv_database = config.get("matsToFasta", "outputPepFasta")

    if not os.path.exists(sv_database):
        raise Exception("Splice variant database not found. \n File %s doesn't exist" % sv_database)

    # check if canonical protein database file exist
    can_database = config.get("Extract", "canonical")

    if not os.path.exists(sv_database):
        raise Exception("Canonical database not found. \n File %s doesn't exist" % can_database)

    # check if evidence.txt and peptides.txt files exist
    evidence = config.get("MaxQuant output", "evidence")

    if not os.path.exists(evidence.strip()):
        raise Exception("MaxQuant output evidence.txt file not found. \n File %s doesn't exist" % evidence)

    peptides = config.get("MaxQuant output", "peptides")

    if not os.path.exists(evidence.strip()):
        raise Exception("MaxQuant output peptides.txt file not found. \n File %s doesn't exist" % peptides)

    # check if output directory exists

    output_dir = config.get("Extract", "output_prefix")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("Created a directory %s" % output_dir)

    return config



# Obtain splice variant specific peptides
def extractVariantPeptides(config):
    arg_maxquant = dict(config.items("MaxQuant output"))
    arg = dict(config.items("Extract"))

    print "#--------------------------------------------------#\n" \
          "|Extracting variant peptides from the evidence file|\n" \
          "#--------------------------------------------------#"

    cmd = "Rscript src/extractVariantPeptides.R " + arg_maxquant["evidence"] + " " + \
          config.get("matsToFasta", "outputPepFasta")  + " " + arg["output_prefix"] + " " + \
          arg["threads"] + " " + arg["canonical"]

    print cmd
    try:
        output = subprocess.check_output(shlex.split(cmd))
        print output

    except:
        raise Exception(" extractVariantPeptides.R failed. Please check the R error. ")

# Add splice variant information on RNA level to the peptide list
def getRNAInformation(config):
    arg = dict(config.items("matsToFasta"))

    print "#----------------------------------------------------------#\n" \
          "|Adding RNA splice variant information to the peptide table|\n" \
          "#----------------------------------------------------------#"

    cmd = "Rscript src/getRNAInformation.R " + config.get("RMATS", "o") + " " + arg["input"] + " " + \
          arg["comparisonName"] + " " + config.get("Extract", "output_prefix") + " " + arg["ifMXE"]

    print cmd
    try:
        output = subprocess.check_output(shlex.split(cmd))
        print output

    except:
        raise Exception(" getRNAInformation.R failed. Please check the R error. ")
        
    print "Isoform-specific peptides were identified."
    print ("You can find the output in the directory %s ." %config.get("Extract", "output_prefix"))



def quantitativeAnalysis(config):
    arg = dict(config.items("Quantitative Analysis"))

    print "#---------------------------------------#\n" \
          "|Running differential peptide expression|\n" \
          "#---------------------------------------#"

    cmd = "Rscript src/quantitativeAnalysis.R " + config.get("MaxQuant output", "peptides") + " " + \
          arg["group1_column_nr"].replace(" ", "") + " " + arg["group2_column_nr"].replace(" ", "") + " " + config.get("Extract", "output_prefix") + \
          " " + arg["imputation"]

    print cmd
    try:
        output = subprocess.check_output(shlex.split(cmd))
        print output

    except:
        raise Exception(" quantitativeAnalysis.R failed. Please check the R error. ")
     
    print "Quantitative analysis of the isoform-specific peptides succeeded."
    print ("You can find the output in the directory %s ." %config.get("Extract", "output_prefix"))
    print "This is the end of the step2 of the proteogenomic pipeline SPLICIFY"   




def main():

    parser = argparse.ArgumentParser(description='Run step 2 of the proteogenomic pipeline. '
                                                 'By default the full analysis will be performed. '
                                                 'If you want to skip some steps use False in the commands: '
                                                 '-e, --extract, -g, --getrna, -q, --quantify, '
                                                 'but make sure the config file has paths to intermediate outputs '
                                                 'and the file names are correct. In case of no replicates, '
                                                 'it is recommended to skip the quantitative analysis',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-e", "--extract", help=" run ExtractVariantPeptides, True/False", default=True, nargs='?')
    parser.add_argument("-g", "--getrna", help=" run getRNAInformation, True/False", default=True, nargs='?')
    parser.add_argument("-q", "--quantify", help=" run quantitativeAnalysis , True/False", default=True, nargs='?')

    requiredNamed = parser.add_argument_group('required argument')
    requiredNamed.add_argument("-c", "--config", help=" configuration file, required (e.g.: src/configPE.ini)",
                               required=True, nargs=1)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        parameters = readInput(args.config)
        if args.extract != "False":
            extractVariantPeptides(parameters)
        if args.getrna != "False":
            getRNAInformation(parameters)
        if args.quantify != "False":
            quantitativeAnalysis(parameters)

if __name__ == "__main__":
    main()
