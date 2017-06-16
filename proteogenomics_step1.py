#########################################################################
#                                                                       #
# Proteogenomic pipeline for differential splice variant identification #
#                                                                       #
#                               27/01/2017                              #
#                                                                       #
#                               MA Komor                                #
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
    for section in ["INPUT FILES", "TRIMMOMATIC", "STAR", "RMATS", "matsToFasta"]:
        if section not in config.sections():
            raise Exception("Section %s missing from config file" % section)

    # check if paired/single-end files exist
    seq = config.get("INPUT FILES", "seq")
    group1 = config.get("INPUT FILES", "Group1").split(",")
    group2 = config.get("INPUT FILES", "Group2").split(",")

    if seq == "paired":

        for sample in group1 + group2:
            pairedend = sample.split(":")
            for files in pairedend:
                if not os.path.exists(files.strip()):
                    raise Exception("File %s doesn't exist" % files)
            if len(pairedend) != 2:
                raise Exception("Matching file  missing for sample %s" % sample)

    elif seq == "single":

        for sample in group1 + group2:
            singleend = sample.split(":")
            for files in singleend:
                if not os.path.exists(files.strip()):
                    raise Exception("File %s doesn't exist" % files)
            if len(singleend) != 1:
                raise Exception("Sample %s single end data should not be separated by :" % sample)
    else:
        raise Exception("%s is a not valid seq argument, choose from single/paired" % seq)

    return config

# Trim the reads
def runTrimmomatic(config):
    print "#-------------------#\n|Running Trimmomatic|\n#-------------------#"
    arg = dict(config.items("TRIMMOMATIC"))

    group1 = config.get("INPUT FILES", "Group1").split(",")
    group2 = config.get("INPUT FILES", "Group2").split(",")

    seq = config.get("INPUT FILES", "seq")

    for group in [group1, group2]:
        for sample_nr in range(1, len(group) + 1):
            sample_r1 = group[sample_nr - 1].split(":")[0]
            if seq == "paired":
                sample_r2 = group[sample_nr - 1].split(":")[1]

            sample_nr = str(sample_nr)

            if group == group1:
                group_name = 'group1'
            else:
                group_name = 'group2'
            if seq == "paired":

                cmd = 'java -jar ' + arg['trimmomatic_jar'] + ' PE -threads ' + arg['threads'] + \
                      ' -phred' + arg['phred'] + ' ' + sample_r1 + ' ' + sample_r2 + ' ' + \
                      os.path.join(arg['output'], 'paired_' + group_name + '_sample' + sample_nr + '_R1.fq.gz ') + \
                      os.path.join(arg['output'], 'unpaired_' + group_name + '_sample' + sample_nr + '_R1.fq.gz ') + \
                      os.path.join(arg['output'], 'paired_' + group_name + '_sample' + sample_nr + '_R2.fq.gz ') + \
                      os.path.join(arg['output'], 'unpaired_' + group_name + '_sample' + sample_nr + '_R2.fq.gz') + \
                      ' ILLUMINACLIP:' + arg['illuminaclip'] + ' LEADING:' + arg['leading'] + ' TRAILING:' + \
                      arg['trailing'] + ' CROP:' + arg['crop'] + ' AVGQUAL:' + arg['avgqual'] + ' SLIDINGWINDOW:' + \
                      arg['slidingwindow'] + ' MINLEN:' + arg['minlen']

            elif seq == "single":

                cmd = 'java -jar ' + arg['trimmomatic_jar'] + ' SE -threads ' + arg['threads'] + ' -phred' + \
                      arg['phred'] + ' ' + sample_r1 + ' ' + \
                      os.path.join(arg['output'], group_name + '_sample' + sample_nr + '.fq.gz ') + ' ILLUMINACLIP:' + \
                      arg['illuminaclip'] + ' LEADING:' + arg['leading'] + ' TRAILING:' + arg['trailing'] + ' CROP:' + \
                      arg['crop'] + ' AVGQUAL:' + arg['avgqual'] + ' SLIDINGWINDOW:' + arg['slidingwindow'] + ' MINLEN:' \
                      + arg['minlen']
            else:

                raise Exception("Wrong sequencing parameter %s. Choose from single/paired" % seq)

            print "\nTrimmomatic command:\n"
            print cmd, '\n'
            try:
                output = subprocess.check_output(shlex.split(cmd))
            except:
                raise Exception(" Trimmomatic failed for sample %s " % group[int(sample_nr - 1)])


            # Align the reads to the genome

# Map the reads
def runSTAR(config):
    print "#------------#\n|Running STAR|\n#------------#"

    arg = dict(config.items("STAR"))
    trim_output = config.get("TRIMMOMATIC", "output")

    group1 = config.get("INPUT FILES", "Group1").split(",")
    group2 = config.get("INPUT FILES", "Group2").split(",")

    for group in [group1, group2]:

        if group == group1:
            group_name = 'group1'
        else:
            group_name = 'group2'

        for sample_nr in range(1, len(group) + 1):
            star_input = ''

            if config.get("INPUT FILES", "seq") == "paired":

                file1 = os.path.join(trim_output, 'paired_' + group_name + '_sample' + str(sample_nr) +
                                     '_R1.fq.gz ')
                file2 = os.path.join(trim_output, 'paired_' + group_name + '_sample' + str(sample_nr) +
                                     '_R2.fq.gz ')

                if not os.path.exists(file1.strip()):
                    raise Exception("Trimmomatic ouput %s does not exist. Check Trimmomatic output in the "
                                    "configuration file." % file1)

                if not os.path.exists(file2.strip()):
                    raise Exception("Trimmomatic ouput %s does not exist. Check Trimmomatic output in the "
                                    "configuration file." % file2)
                star_input = file1 + ' ' + file2

            elif config.get("INPUT FILES", "seq") == "single":

                star_input = os.path.join(trim_output, group_name + '_sample' + str(sample_nr) + '.fq.gz ')

                if not os.path.exists(star_input.strip()):
                    raise Exception("Trimmomatic ouput %s does not exist. Check Trimmomatic output "
                                    "in the configuration file." % star_input)

            cmd = arg["pathToStar"]

            for key, value in arg.items():
                if key == "outFileNamePrefix":
                    cmd = cmd + " --" + key + " " + os.path.join(value, group_name + "_sample" + str(sample_nr))
                elif key != "pathToStar":
                    cmd = cmd + " --" + key + " " + value

            cmd = cmd + " --readFilesIn " + star_input

            print "\nSTAR command\n"
            print cmd
            try:
                output = subprocess.check_output(shlex.split(cmd))
            except:
                raise Exception(" STAR failed for sample %s. Please check the STAR error. " % group[int(sample_nr)])

# Run rMATS
def runrMATS(config):
    print "#-------------#\n|Running rMATS|\n#-------------#"

    arg = dict(config.items("RMATS"))

    group1 = config.get("INPUT FILES", "Group1").split(",")
    group2 = config.get("INPUT FILES", "Group2").split(",")

    bams_dir = config.get("STAR", "outFileNamePrefix")

    b1 = []
    b2 = []

    for sample_nr in range(1, len(group1) + 1):
        b1.append(os.path.join(bams_dir, 'group1_sample' + str(sample_nr) + 'Aligned.sortedByCoord.out.bam'))

    for sample_nr in range(1, len(group2) + 1):
        b2.append(os.path.join(bams_dir, 'group2_sample' + str(sample_nr) + 'Aligned.sortedByCoord.out.bam'))

    if arg['analysis'] != 'P' and arg['analysis'] != 'U':
        print 'WARNING: Wrong analysis argument - %s, choose from U/P. ' \
              'Running unpaired rMATS analysis (U). ' % arg['analysis']
        arg['analysis'] = 'U'

    if not os.path.exists(arg['gtf']):
        raise Exception("No gtf file. File %s doesn not exist" % arg['gtf'])

    cmd = 'python ' + os.path.join(arg['pathTorMATS'], 'RNASeq-MATS.py') + ' -b1 ' + ','.join(b1) + ' -b2 ' + \
          ','.join(b2) + ' -gtf ' + arg['gtf'] + ' -t ' + config.get("INPUT FILES", "seq") + ' -len ' + \
          config.get("TRIMMOMATIC", "minlen") + ' -o ' + arg['o'] + ' -analysis ' + arg['analysis']

    print cmd
    try:
        output = subprocess.check_output(shlex.split(cmd))
    except:
        raise Exception(" rMATS failed. Please check the rMATS error. ")

# Obtain splice variant database
def matsToFasta(config):
    arg = dict(config.items("matsToFasta"))
    print "#-----------------------------------------------#\n|Obtaining fasta sequences from the rMATS output|\n" \
          "#-----------------------------------------------#"
    print "#-------------------------------#\n|Change rMATS output to bed file|\n#-------------------------------#"
    cmd = "Rscript src/matsToBed.R " + os.path.join(config.get("RMATS", "o"), "MATS_output/") + " " + \
          arg["comparisonName"] + " " + \
          arg["outputBed"] + " " + arg["ifMXE"] + " " + arg["input"] + " " + arg["fdr"]

    print cmd
    try:
        output = subprocess.check_output(shlex.split(cmd))
    except:
        raise Exception(" matsToBed.R failed. Please check the R error. ")
    #print output

    print "#----------------------#\n|Bed file to Fasta file|\n#----------------------#"

    # check if genoma fasta sequence available
    if not os.path.exists(arg["Genome_Fasta"]):
        raise Exception ("Fasta file %s does not exist. \nProvide genome fasta file to obtain fasta sequences from "
                         "the splice variants." %arg["Genome_Fasta"])

    cmd = "bedtools getfasta -s -name -split -fi " + arg["Genome_Fasta"] + " -bed " + arg["outputBed"] + \
          " -fo " + arg["outputFasta"]

    print cmd
    try:
        output = subprocess.check_output(shlex.split(cmd))
    except:
        raise Exception(" Bed to Fasta failed. Please check the error. ")
    print output

    print "#-----------------------------------#\n|Nucleotide sequences to amino acids|\n" \
          "#-----------------------------------#"

    if arg["pathToTranseq"] == "FALSE":
        cmd = "transeq " + arg["outputFasta"] + " " + arg["outputPepFasta"] + " -frame=F"
    else:
        cmd = arg["pathToTranseq"] + " " + arg["outputFasta"] + " " + arg["outputPepFasta"] + " -frame=F"

    print cmd

    try:
        output = subprocess.check_output(shlex.split(cmd))
    except (IOError, OSError) as e:
        if e.errno == os.errno.ENOENT:
            print "Transeq not installed."
        print "Running alternative script for nucleotide to amino acid translation." \
              " Check src/threeFrame.R for details."

        cmd = "Rscript src/threeFrame.R " + arg["outputFasta"] + " " + arg["outputPepFasta"]
        try:
            output = subprocess.check_output(shlex.split(cmd))
        except:
            raise Exception ("Something went wrong with the src/threeFrame.R script")

    print output
    print "The genomics part of the proteogenomic pipeline is finished."
    print ("You can find the splice variant database here %s ." %arg["outputPepFasta"])
    print ("Please proceed to the proteomics analysis with the new database.")

def main():
    parser = argparse.ArgumentParser(description='Run step 1 of the proteogenomic pipeline. '
                                                 'By default the full analysis will be performed. '
                                                 'If you want to skip some steps use False in the commands: '
                                                 '-t, --trim, -s, --star, -r, --rmats, -m, --matstofasta'
                                                 ', but make sure the config file has paths to intermediate outputs '
                                                 'and the file names are correct.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-t", "--trim", help=" run trimmomatic, True/False", default=True, nargs='?')
    parser.add_argument("-s", "--star", help=" run STAR, True/False", default=True, nargs='?')
    parser.add_argument("-r", "--rmats", help=" run rMATS, True/False", default=True, nargs='?')
    parser.add_argument("-m", "--matstofasta", help=" run matsToFasta, True/False", default=True, nargs='?')

    requiredNamed = parser.add_argument_group('required argument')
    requiredNamed.add_argument("-c", "--config", help=" configuration file, required (e.g.: src/configPE.ini)",
                               required=True, nargs=1)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        parameters = readInput(args.config)
        if args.trim != "False":
            runTrimmomatic(parameters)
        if args.star != "False":
            runSTAR(parameters)
        if args.rmats != "False":
            runrMATS(parameters)
        if args.matstofasta != "False":
            matsToFasta(parameters)

if __name__ == "__main__":
    main()
