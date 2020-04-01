#!/usr/bin/python
# Copyright (C) 2019 Ion Torrent Systems, Inc. All Rights Reserved
#
# Plugin: AssemblerTrinity
# This plugin is developed for Ion Ampliseq Coronavirus panel sequencing data assembly
#
# Author: Lucius Zheng
# Last modified: 2020/03/07
#

import json
import os
from django.utils.functional import cached_property
from ion.plugin import *
import subprocess
from subprocess import check_output
import shutil

from Bio import SeqIO
from Bio.Seq import Seq

from django.conf import settings
from django.template.loader import render_to_string

def createReport(reportName,reportTemplate,reportData):
    with open(reportName,'w') as bcsum:
        bcsum.write( render_to_string(reportTemplate,reportData) )

class AssemblerTrinity(IonPlugin):
    # The version number for this plugin
    version = "2.0.0.0"

    # this plugin can run on fullchip runs, thumbnail runs, and composite (merged via project page) runs
    # note that when the plugin is manually launched, only the 'launch' method will be called
    runtypes = [RunType.FULLCHIP, RunType.THUMB, RunType.COMPOSITE]

    # specify when the plugin is called.  For log parsing, stay simple and just get called when the run completes.
    # but can also be called before the run starts, at the block level, or after all other default plugins run
    runlevels = [RunLevel.DEFAULT]
    
    # a simple cached version of the start plugin property
    @cached_property
    def startplugin_json(self):
        return self.startplugin

    @cached_property
    def barcodes_json(self):
        with open('barcodes.json', 'r') as barcodes_handle:
            return json.load(barcodes_handle)
    
    def launch(self, data=None):
        """This is the primary launch method for the plugin."""
    
        # configure django to use the templates folder        
        #settings.configure(TEMPLATE_DIRS=(self.startplugin["runinfo"]["plugin_dir"] + '/templates'),)
        
        if not settings.configured:
            settings.configure( DEBUG=False, TEMPLATE_DEBUG=False,
                                INSTALLED_APPS=('django.contrib.humanize',),
                                TEMPLATE_DIRS=(os.path.join(self.startplugin["runinfo"]["plugin_dir"],'templates'),) 
                            )
        
        
        # define Trinity Singularity image
        trinity_image = self.startplugin_json['runinfo']['plugin_dir'] + "/bin/trinityrnaseq.v2.9.1.simg"
        
        # start to analyze bam files
        
        Trinity_report = []
        for barcode_name, barcode_values in self.barcodes_json.iteritems():
            # do you work per barcode here!    
            
            # first check to see if the barcode was excluded using the frame work barcodes configuration table        
            selected = True
            barcodeData = self.startplugin_json['pluginconfig'].get('barcodetable',None)
            if barcodeData:
                #print(barcodeData)
                for bc in barcodeData:
                    if  bc.get('barcode_name',"") == barcode_name:
                        selected = bc.get('selected',True)
                        break
            
            if not selected:
                continue
            
            print("Barcode Name: " + barcode_name)
            print("Bam Filepath: " + barcode_values['bam_filepath'])
            print("Read count: " + str(barcode_values['read_count']))
            
            # if no BAM file or file size is 0, then skip the sample
            if not os.path.exists(barcode_values['bam_filepath']):
                print "BAM file does not exist. We will skip the sample in the followed analysis.\n"
                continue
            if os.path.getsize(barcode_values['bam_filepath']) == 0:
                print "BAM file size is zero. We will skip the sample in the followed analysis.\n"
                continue
                        
            if barcode_values['read_count'] < int(self.startplugin_json['pluginconfig']['min_reads']):
                print "BAM file reads number is less than minimum required. We will skip the sample in the followed analysis.\n"
                continue
            
            # define variables
            results_dir = self.startplugin_json['runinfo']['results_dir']
            ram = self.startplugin_json['pluginconfig']['RAM']
            cpu = self.startplugin_json['pluginconfig']['CPU']
            min_kmer_cov = self.startplugin_json['pluginconfig']['min_kmer_cov']
            
            if self.startplugin_json['pluginconfig']['genome_guided'] == "Yes": 
            
                # copy bam file to current dir, ensuring it's accessible in the singularity container
                cp_results = check_output(["cp",barcode_values['bam_filepath'], barcode_values['bam_file']], cwd=results_dir)
            
                try:
                    check_output(["singularity", "exec", "-e", trinity_image, "Trinity", 
                            "--genome_guided_bam", barcode_values['bam_file'], 
                            "--genome_guided_max_intron", "1000",
                            "--max_memory", ram,
                            "--CPU", cpu,
                            "--min_kmer_cov", min_kmer_cov,
                            "--output", barcode_name + "_trinity_out"],
                            cwd=results_dir)
                except subprocess.CalledProcessError as e:
                    print e.output
            
                # then delete it when finished to save space
                rm_results = check_output(["rm",barcode_values['bam_file']], cwd=results_dir)
            
                Trinity_output = os.path.join(results_dir, barcode_name + "_trinity_out", "Trinity-GG.fasta")
            
            else:
                # convert sam/bam to fastq
                fastq_filename = barcode_name + '.fastq'
                cmd = "samtools bam2fq " + barcode_values['bam_filepath'] + ">" + fastq_filename
                fastq_results = check_output(cmd, shell=True, cwd=results_dir)

                try:
                    check_output(["singularity", "exec", "-e", trinity_image, "Trinity", 
                            "--seqType", "fq", 
                            "--single", fastq_filename,
                            "--max_memory", ram,
                            "--CPU", cpu,
                            "--min_kmer_cov", min_kmer_cov,
                            "--output", barcode_name + "_trinity_out"],
                            cwd=results_dir)
                except subprocess.CalledProcessError as e:
                    print e.output
                
                Trinity_output = os.path.join(results_dir, barcode_name + "_trinity_out", "Trinity.fasta")
            
            max_len = 1
            longest_seq = ""
            
            data_entry = {}
            with open(Trinity_output, "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if len(record.seq) > max_len:
                        max_len = len(record.seq)
                        longest_seq = record
            
            longest_output = os.path.join(results_dir, barcode_name + "_trinity_out", barcode_name + "_Trinity_longest.fasta")
            SeqIO.write(longest_seq, longest_output, "fasta")
            
            ###########################################################################
            # We further do 2 rounds of align-call-refine cycles to polish the contig
            ###########################################################################
            
            trinity_outdir = cwd=os.path.join(results_dir,barcode_name + "_trinity_out")
            
            # select the parameters used for TVC
            if self.startplugin_json['pluginconfig']['chip'] == "Proton P1":
                parameters = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/ampliseq_germline_cov_p1_parameters.json"
                
            if self.startplugin_json['pluginconfig']['chip'] == "PGM":
                parameters = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/ampliseq_germline_cov_pgm_parameters.json"
            
            if self.startplugin_json['pluginconfig']['chip'] == "510" or self.startplugin_json['pluginconfig']['chip'] == "520" or self.startplugin_json['pluginconfig']['chip'] == "530":
                parameters = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/ampliseq_germline_cov_530_parameters.json"
            
            if self.startplugin_json['pluginconfig']['chip'] == "540":
                parameters = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/ampliseq_germline_cov_540_parameters.json"
            
            if self.startplugin_json['pluginconfig']['chip'] == "550":
                parameters = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/ampliseq_germline_cov_550_parameters.json"
            
            # 1st align-call cycle: TMAP aligns the data to assembled longest contig with human transcriptome, then TVC calls variants
            
            hg19_trxo = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/hg19_human_trxo.fasta"
            check_output(["cat",hg19_trxo,">>",longest_output],shell=True, cwd=trinity_outdir)
            
            build_cmd = "build_genome_index.pl -f " + longest_output + " -s trinity_longest -l  trinity_longest -v 1.0" 
            build_results = check_output(build_cmd, shell=True, cwd=trinity_outdir)
            reference = "trinity_longest/trinity_longest.fasta"
            
            target_bed = os.path.join(results_dir, barcode_name + "_trinity_out","target.bed")
            with open(target_bed,"w") as file:
                file.write(longest_seq.id + "\t\0\t" + str(max_len) + "\n")
            
            ## tmap
            aligned_filename = barcode_name + '.1.bam'
            align_cmd = "tmap mapall -n " + cpu + " -f " + reference + \
                        " -r " + barcode_values['bam_filepath'] + \
                        " -v -Y -u --prefix-exclude 5 -o 1 -J 25 --end-repair 15 --do-repeat-clip --context stage1 map4" + \
                        " | samtools sort -T align.sorted -@ " + cpu + " -o " + aligned_filename
            
            align_results = check_output(align_cmd, shell=True, cwd=trinity_outdir)
            
            ## tvc
            call_cmd = "/results/plugins/variantCaller/bin/variant_caller_pipeline.py --input-bam " + aligned_filename + \
                        " --reference-fasta " + reference + \
                        " --output-dir variantCall1 " + \
                        " --parameters-file " + parameters + \
                        " --region-bed " + target_bed + \
                        " --bin-dir /results/plugins/variantCaller/bin " + \
                        " --num-threads " + cpu
            
            call_results = check_output(call_cmd, shell=True, cwd=trinity_outdir)
            
            ## refine the reference using major alleles
            
            refine_cmd = "cat " + reference + " | vcf-consensus variantCall1/TSVC_variants.vcf.gz > variantCall1/TSVC_variants.consensus.fasta"
            refine_results = check_output(refine_cmd, shell=True, cwd=trinity_outdir)
            
            # 2nd align-call cycle to minimize reference bias in the assembly
            ## first build genome index on the new reference
            build_cmd = "build_genome_index.pl -f variantCall1/TSVC_variants.consensus.fasta -s trinity_polished1 -l trinity_polished1 -v 1.0" 
            build_results = check_output(build_cmd, shell=True, cwd=trinity_outdir)
            
            trinity_polished1 = os.path.join(trinity_outdir, "variantCall1/TSVC_variants.consensus.fasta")
            with open(trinity_polished1, "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id == longest_seq.id:
                        max_len = len(record.seq)
            
            with open(target_bed,"w") as file:
                file.write(longest_seq.id + "\t\0\t" + str(max_len) + "\n")
            
            ## 2nd tmap
            aligned_filename = barcode_name + '.2.bam'
            align_cmd = "tmap mapall -n " + cpu + " -f " + "trinity_polished1/trinity_polished1.fasta" + \
                        " -r " + barcode_values['bam_filepath'] + \
                        " -v -Y -u --prefix-exclude 5 -o 1 -J 25 --end-repair 15 --do-repeat-clip --context stage1 map4" + \
                        " | samtools sort -T align.sorted -@ " + cpu + " -o " + aligned_filename
            
            align_results = check_output(align_cmd, shell=True, cwd=trinity_outdir)
            
            ## 2nd tvc
            
            call_cmd = "/results/plugins/variantCaller/bin/variant_caller_pipeline.py --input-bam " + aligned_filename + \
                        " --reference-fasta " + "trinity_polished1/trinity_polished1.fasta" + \
                        " --output-dir variantCall2 " + \
                        " --parameters-file " + parameters + \
                        " --region-bed " + target_bed + \
                        " --bin-dir /results/plugins/variantCaller/bin " + \
                        " --num-threads " + cpu
            
            call_results = check_output(call_cmd, shell=True, cwd=trinity_outdir)
            
            ## refine the reference again using major alleles to get a consensus genome
            
            # extract 2019-nCoV consensus sequence in the 1st round
            consensus1 = trinity_outdir + "/variantCall1/TSVC_variants.consensus.2019-nCoV.fasta"
            with open(trinity_outdir + "/variantCall1/TSVC_variants.consensus.fasta", "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id == longest_seq.id:
                        SeqIO.write(record, consensus1, "fasta")
                        break
            
            refine_cmd = "cat " + consensus1 + " | vcf-consensus variantCall2/TSVC_variants.vcf.gz > variantCall2/TSVC_variants.consensus.fasta"
            refine_results = check_output(refine_cmd, shell=True, cwd=trinity_outdir)
            
            consen_seq = SeqIO.read(trinity_outdir + "/variantCall2/TSVC_variants.consensus.fasta","fasta")
            consen_seq.id = barcode_name + "-Trinity-Polished"
            consensus_output = os.path.join(trinity_outdir, barcode_name + ".trinity.polished.fasta")
            SeqIO.write(consen_seq, consensus_output, "fasta")
            
            # delete intermediate files
            os.remove(os.path.join(trinity_outdir, barcode_name + '.1.bam'))
            os.remove(os.path.join(trinity_outdir, barcode_name + '.1.bam.bai'))
            shutil.rmtree(os.path.join(trinity_outdir,"trinity_longest"))
            shutil.rmtree(os.path.join(trinity_outdir,"trinity_polished1"))
            
            data_entry['barcode'] = barcode_name
            data_entry['sample'] = barcode_values['sample']
            data_entry['longest'] = os.path.join(barcode_name + "_trinity_out", barcode_name + ".trinity.polished.fasta")
            if self.startplugin_json['pluginconfig']['genome_guided'] == "Yes": 
                data_entry['all'] =  os.path.join(barcode_name + "_trinity_out", "Trinity-GG.fasta")
            else:
                data_entry['all'] =  os.path.join(barcode_name + "_trinity_out", "Trinity.fasta")
            
            Trinity_report.append(data_entry)
            
        ###################################################################################################################
        ## output in HTML
        ###################################################################################################################
        
        render_context = {
            "autorefresh" : False,
            "TrinityData" : Trinity_report
        }
        
        createReport(os.path.join(results_dir,'Trinity_report.html'), 'barcode_summary_all.html', render_context )
        
        return True
    
    # Return list of columns you want the plugin table UI to show.
    # Columns will be displayed in the order listed.
    def barcodetable_columns(self):
        return [
            { "field": "selected", "editable": True },
            { "field": "barcode_name", "editable": False },
            { "field": "sample", "editable": False } ]

# Devel use - running directly
if __name__ == "__main__":
    PluginCLI()
