__author__ = "Inodb, Alneberg"
__license__ = "MIT"


configfile: "config.json"

import os
import glob


# Add all reads to the bowtie2_rules
config["bowtie2_rules"]["samples"] = {os.path.basename(p).replace("_R1.cut.sync.fasta", ""): [os.path.basename(p).replace("_R1.cut.sync.fasta", "")] for p in glob.glob("{0}/*_R1.cut.sync.fasta".format( config["reads_dir"]))}
config["bowtie2_rules"]["units"] = {os.path.basename(p).replace("_R1.cut.sync.fasta", ""): [p, p.replace("R1", "R2")] for p in glob.glob("{0}/*_R1.cut.sync.fasta".format( config["reads_dir"] )) }

# Add the assembly to bowtie2_rules
config["bowtie2_rules"]["references"] = {os.path.basename(p).replace(".fna", ""): p for p in config["assemblies"]}

# add newbler merged assemblies to concoct assemblies
config["concoct_rules"]["assemblies"] = {k:v for k,v in config["bowtie2_rules"]["references"].items()}

# Show all bowtie2 logs and markduplicate percent duplications in the mapping report
config["mapping_report_rules"]["bowtie2_logs"] = sorted(expand("mapping/bowtie2/{mapping_params}/{reference}/units/{unit}.log",
                                                        mapping_params=config["bowtie2_rules"]["mapping_params"],
                                                        reference=config["bowtie2_rules"]["references"],
                                                        unit=config["bowtie2_rules"]["units"]))

config["mapping_report_rules"]["markduplicates_metrics"] = sorted(expand("mapping/bowtie2/{mapping_params}/{reference}/units/{unit}.sorted.removeduplicates.metrics",
                                                                  mapping_params=config["bowtie2_rules"]["mapping_params"],
                                                                  reference=config["bowtie2_rules"]["references"],
                                                                  unit=config["bowtie2_rules"]["units"]))

#SM_WORKFLOW_LOC="https://raw.githubusercontent.com/alneberg/snakemake-workflows/fe913ac3a40387dbe26558ac35c8a807236e466a/"
SM_WORKFLOW_LOC="/proj/b2014214/repos/snakemake-workflows/"
#SM_WORKFLOW_LOC="https://raw.githubusercontent.com/inodb/snakemake-workflows/fe913ac3a40387dbe26558ac35c8a807236e466a/"
#SM_WORKFLOW_LOC = "/glob/inod/github/snakemake-workflows/"
include: SM_WORKFLOW_LOC + "common/rules/track_dir.rules"
include: SM_WORKFLOW_LOC + "bio/ngs/rules/assembly/report.rules"
include: SM_WORKFLOW_LOC + "bio/ngs/rules/mapping/bowtie2.rules"
include: SM_WORKFLOW_LOC + "bio/ngs/rules/mapping/samtools.rules"
include: SM_WORKFLOW_LOC + "bio/ngs/rules/mapping/report.rules"
include: SM_WORKFLOW_LOC + "bio/ngs/rules/binning/concoct.rules"
include: SM_WORKFLOW_LOC + "bio/ngs/rules/annotation/prodigal.rules"
include: SM_WORKFLOW_LOC + "bio/ngs/rules/blast/rpsblast.rules"
include: SM_WORKFLOW_LOC + "bio/ngs/rules/annotation/hmmer.rules"
localrules: track_changes


rule report:
    input:
        "report/mapping/index.html",
        "report/concoct/index.html",
        "report/notebooks_output/bin_overview.html"
    output:
        "report/index.html"
    shell:
        """
        (
            echo '<html><head><style>body {{ text-align: center }}</style></head><body>'
            echo "<a href='fastqc/index.html'>FastQC Results</a><br />"
            echo "<a href='assemblies/index.html'>Assembly Results</a><br />"
            echo "<a href='mapping/index.html'>Mapping Results</a><br />"
            echo "<a href='concoct/index.html'>CONCOCT Results</a><br />"
            echo "<a href='notebooks_output/bin_overview.html'>Binning Overview</a><br />"
            echo "<a href='http://nbviewer.ipython.org/urls/github.com/inodb/2014-05-mdopson-viral/tree/master/notebooks'>Notebooks</a><br />"
            echo '</body></html>'
        ) > {output}
        """

rule track_changes:
    input:
        "results_track.txt"
