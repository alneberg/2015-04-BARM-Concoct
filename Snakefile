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

SM_WORKFLOW_LOC="https://raw.githubusercontent.com/alneberg/snakemake-workflows/a2a6cab285adad51f1997a03638cdf53b07b529b/"
#SM_WORKFLOW_LOC="/proj/b2014214/repos/snakemake-workflows/"
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

# Reads are fasta
ruleorder: bowtie2_map_fasta > bowtie2_map > sam_to_bam

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

#  add regular assemblies for prodigal to predict genes for
for a_name, a in config["concoct_rules"]["assemblies"].items():
    config["prodigal_rules"]["assemblies"][a_name] = a 

#  add prodigal predicted genes as query for rpsblast
config["rpsblast_rules"]["query_aas"] = {a: "annotation/prodigal/default-meta/{a}/proteins/proteins.faa".format(a=a) for a in config["prodigal_rules"]["assemblies"]}

#  add prodigal predicted genes as query for hmmer
config["hmmer_rules"]["query_aas"] = config["rpsblast_rules"]["query_aas"]

rule merge_concoct_results:
    input:
        "concoct/{assembly}/output/{concoct_params}/clustering.csv"
    output:
        "concoct/{assembly}/output/{concoct_params}/clustering_merged.csv"
    shell:
        """
            {config[concoct_rules][load_env]}
            python scripts/majority_merge_cutup_clustering.py {input} > {output}
        """

rule concoct_eval_cog_table_merged:
    """
    Generate COG table from rpsblast output and concoct binning results
    """
    input:
        clust="concoct/{assembly}/output/{concoct_params}/clustering_merged.csv",
        rpsblast="blast/rpsblast/default-concoct/cog/{assembly}/rpsblast.out"
    output:
        "concoct/{assembly}/evaluation/scg/{concoct_params}/clustering_scg_merged.tsv"
    shell:
        """
        {config[concoct_rules][load_env]}
        python {config[concoct_rules][scripts_dir]}/COG_table.py \
            -b {input.rpsblast} \
            -m {config[concoct_rules][scripts_dir]}/../scgs/scg_cogs_min0.97_max1.03_unique_genera.txt \
            -c {input.clust} \
            --cdd_cog_file {config[concoct_rules][scripts_dir]}/../scgs/cdd_to_cog.tsv \
            > {output}
        """

rule concoct_eval_cog_plot_merged:
    """
    Plot COGs using COG table
    """
    input:
        "concoct/{assembly}/evaluation/scg/{concoct_params}/clustering_scg_merged.tsv"
    output:
        "concoct/{assembly}/evaluation/scg/{concoct_params}/clustering_scg_merged.pdf"
    shell:
        """
        {config[concoct_rules][load_env]}
        Rscript {config[concoct_rules][scripts_dir]}/COGPlot.R \
            -s {input} \
            -o {output}
        """

rule concoct_eval_all:
    """
    Plot COGs using COG table for both merged and cutup
    """
    input:
        expand("concoct/{assembly}/evaluation/scg/{concoct_params}/clustering_scg_merged.pdf",
                assembly=config["concoct_rules"]["assemblies"],
                concoct_params = config["concoct_rules"]["concoct_params"]),
        expand("concoct/{assembly}/evaluation/scg/{concoct_params}/clustering_scg.pdf",
                assembly=config["concoct_rules"]["assemblies"],
                concoct_params = config["concoct_rules"]["concoct_params"])

rule concoct_extract_approved_scg_bins_merged:
    input:
        scg_tsvs=expand("concoct/{assembly}/evaluation/scg/{concoct_params}/clustering_scg_merged.tsv",
            assembly=sorted(config["concoct_rules"]["assemblies"]),
            concoct_params=sorted(config["concoct_rules"]["concoct_params"])),
        asms=config["assemblies"]
    output:
        dynamic("concoct/approved_merged_scg_bins/{cluster_name}.fa")
    params:
        names=expand("{assembly}_{concoct_params}",
            assembly=sorted(config["concoct_rules"]["assemblies"]),
            concoct_params=sorted(config["concoct_rules"]["concoct_params"])),
        groups=expand("{assembly}",
            assembly=sorted(config["concoct_rules"]["assemblies"]),
            concoct_params=sorted(config["concoct_rules"]["concoct_params"]))
    shell:
        """
            {config[concoct_rules][load_env]}
            python {config[concoct_rules][scripts_dir]}/extract_scg_bins.py \
                --output_folder concoct/approved_merged_scg_bins \
                --scg_tsvs {input.scg_tsvs} \
                --fasta_files {input.asms} \
                --names {params.names} \
                --groups {params.groups} \
                --max_missing_scg 5 \
                --max_multicopy_scg 2
         """

rule concoct_extract_approved_scg_bins_all_merged:
    input:
        dynamic("concoct/approved_merged_scg_bins/{cluster_name}.fa")


rule concoct_dnadiff_dist_matrix_merged:
    """Get distance matrix from approved SCG bins"""
    input:
        clusters=dynamic("concoct/approved_merged_scg_bins/{cluster_name}.fa")
    output:
        "concoct/dnadiff_dist_matrix_merged/dist_matrix.tsv",
        "concoct/dnadiff_dist_matrix_merged/hclust_heatmap.pdf",
        "concoct/dnadiff_dist_matrix_merged/hclust_dendrogram.pdf"
    run:
        mags="concoct/mags/*.fa"
        sorted_input = sorted(input.clusters)
        shell("""
        {config[concoct_rules][load_env]}
        python {config[concoct_rules][scripts_dir]}/dnadiff_dist_matrix.py \
            concoct/dnadiff_dist_matrix_merged {sorted_input} {mags}
        """)

rule concoct_dnadiff_dist_matrix_mags:
    """Get distance matrix from approved SCG bins"""
    input:
        clusters=dynamic("concoct/approved_scg_bins/{cluster_name}.fa")
    output:
        "concoct/dnadiff_dist_matrix_mags/dist_matrix.tsv",
        "concoct/dnadiff_dist_matrix_mags/hclust_heatmap.pdf",
        "concoct/dnadiff_dist_matrix_mags/hclust_dendrogram.pdf"
    run:
        mags="concoct/mags/*.fa"
        sorted_input = sorted(input.clusters)
        shell("""
        {config[concoct_rules][load_env]}
        python {config[concoct_rules][scripts_dir]}/dnadiff_dist_matrix.py \
            concoct/dnadiff_dist_matrix_mags {sorted_input} {mags}
        """)


rule track_changes:
    input:
        "results_track.txt"
