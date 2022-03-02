#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from hl_functions import *
hl.init(log='./log.log')
#hl.init(master='spark://c23n02.ruddle.hpc.yale.internal:7077')


# In[ ]:


import csv
import pandas as pd
pd.set_option("display.max_columns", None)
pd.set_option("display.max_colwidth", None)
from pcgc_hail.hail_scripts.utils import *
#from new_names import *
import gnomad
from gnomad_methods import *
from gnomad_methods.gnomad.sample_qc import *

def hl_to_txt(hl_df, name, delim='\t'):
    """Convert table to pandas dataframe and output to file"""
    df = hl_df.to_pandas()
    df.to_csv(name, sep=delim)

#from matplotlib import pyplot as plt
    
import bokeh.io
from bokeh.io import * 
from bokeh.layouts import *
from bokeh.models import *
hl.plot.output_notebook()


# # User Updates 

# In[ ]:


# Update with desired paths 

# Make sure these directories have a slash at the end
hl_outdir = hl_outdir = "/gpfs/gibbs/pi/brueckner/ken_rotation/output/consortium_new/"
#dnv_outdir = 

######################################################################################################

#Yale generated Joint VCF based on Yale gvcfs (build 37)
#vcf = '/gpfs/gibbs/pi/brueckner/Phil_Pipeline/yale_callset/yale_pipeline.vcf.gz'

#Yale generated Joint VCF based on ACC gvcfs (build 38)
#vcf = '/gpfs/gibbs/pi/brueckner/Phil_Pipeline/callset/phil_pipeline.vcf.gz'

#ACC generated Joint VCF based on ACC gvcfs (build 38)
#vcf = '/gpfs/gibbs/pi/brueckner/Phil_Joint_Calling_Pipeline/pilotcalls.vqsr.vcf.gz'

#unsplit_mt = hl_outdir + 'unsplit.mt' 
split_mt = hl_outdir + 'split.mt' 
vep_mt = hl_outdir + 'vep.mt' 
vep_r1_mt = hl_outdir + 'vep_r1.mt' 
vep_r1c1_mt = hl_outdir + 'vep_r1c1.mt'


# # VCF -> Matrix Table (unsplit) -> Matrix Table (split)

# In[ ]:


# VCF to hail matrix table 

#vcf = ###

################################################################################
print('Reads in:' + vcf)
print('Writes out:' + mt_output)

# if using genome 37
#mt = hl.import_vcf(vcf, force_bgz=True, reference_genome='GRCh37')

# if using genome 38
# for whether to use recode see https://hail.is/docs/0.2/methods/impex.html#hail.methods.import_vcf
recode = {f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
mt = hl.import_vcf(vcf, force_bgz=True, reference_genome='GRCh38', contig_recoding=recode)

# count before splitting 
print(str(mt.count()) + '\n')

# split multi-allelic sites  
mt = generate_split_alleles(mt)

# count after splitting 
print(str(mt.count()) + '\n')

mt.write(mt_output, overwrite = True)


# In[ ]:


#Restrict to the exome intervals
mt_input = "/gpfs/gibbs/pi/brueckner/ken_rotation/output/consortium/split2.mt"
mt_output = "/gpfs/gibbs/pi/brueckner/ken_rotation/output/consortium_new/split_intervals_acc.mt"
intervals = "/gpfs/gibbs/pi/brueckner/ken_rotation/intervals_3.bed"   #Intersected bed file
exome_intervals = hl.import_locus_intervals(intervals, reference_genome='GRCh38')
mt = hl.read_matrix_table(mt_input)
#mt = mt.annotate_rows(old_locus=mt.locus) #If lifting over

# Liftover from hg19 to hg38
#rg37 = hl.get_reference('GRCh37')  
#rg38 = hl.get_reference('GRCh38')  
#rg37.add_liftover("/gpfs/gibbs/pi/brueckner/ken_rotation/references_grch37_to_grch38.over.chain.gz", rg38) 
#mt = mt.annotate_rows(new_locus=hl.liftover(mt.locus, 'GRCh38'))  
#mt = mt.filter_rows(hl.is_defined(mt.new_locus))  
#mt = mt.key_rows_by(locus=mt.new_locus, alleles=mt.alleles)  

mt = mt.filter_rows(hl.is_defined(exome_intervals[mt.locus]), keep=True)
mt.write(mt_output, overwrite = True)

# Generate sample QC metrics
mt_input = "/gpfs/gibbs/pi/brueckner/ken_rotation/output/consortium_new/split_intervals_acc.mt"
mt_output = "/gpfs/gibbs/pi/brueckner/ken_rotation/output/consortium_new/split_sampleqc_acc.mt"
mt = hl.read_matrix_table(mt_input)
mt = hl.sample_qc(mt, name='sample_qc')
mt.write(mt_output, overwrite = True)


# In[ ]:


get_ipython().run_cell_magic('capture', '', '\nvep_json = \'/gpfs/gibbs/pi/brueckner/14_v0/vep85-loftee-ruddle-b38.json\'\n\n#mt_input ="/gpfs/gibbs/pi/brueckner/ken_rotation/output/consortium_new/split_intervals_acc.mt"\n#mt_output ="/gpfs/gibbs/pi/brueckner/ken_rotation/output/consortium_new/vep_acc.mt"\n\n#print(\'Reads in the following matrix table: \\n\' + mt_input + \'\\n\')\nmt = hl.read_matrix_table(mt_input)\n\nprint(\'Writes out to the following matrix table: \\n\' + mt_output + \'\\n\')\n\nmt_vep = hl.vep(mt, vep_json)\nmt_vep.write(mt_output, overwrite = True)')


# In[ ]:


### Filter out PCGC4 Samples ###
rf = open("/gpfs/gibbs/pi/brueckner/ken_rotation/pilot_info_for_ken.tsv", 'r+')   #Open the file for reading and writing
rf_line = []   #Initialize a list
for line in rf:   #Interate through all lines
    line_s = line.split('\t')   #Split with a tab delimeter
    rf_line.append(line_s)   #Append each line to the initialized list
rf_line = rf_line[1:]
new_lis = []  #LIST OF NON PCGC4

for i in rf_line:
    batch.append(i[3])
    if i[3]!="NPCGC\n":
        new_lis.append(i[0])


# In[ ]:


mt_acc = "/gpfs/gibbs/pi/brueckner/ken_rotation/output/consortium_new/vep_acc.mt"
mt_a = hl.read_matrix_table(mt_acc)
mt_yale = "/gpfs/gibbs/pi/brueckner/ken_rotation/output/consortium_new/vep.mt"
mt_y = hl.read_matrix_table(mt_yale)


# In[ ]:


#Remove PCGC4 samples
mt_a = mt_a.filter_cols(hl.set(new_lis).contains(mt_a.s))
mt_y = mt_y.filter_cols(hl.set(new_lis).contains(mt_y.s))


# In[ ]:


mt_a = mt_a.filter_rows(hl.agg.any(mt_a.GT.is_non_ref()))
mt_y = mt_y.filter_rows(hl.agg.any(mt_y.GT.is_non_ref()))


# In[ ]:


mt_a = hl.sample_qc(mt_a, name='sample_qc')
mt_y = hl.sample_qc(mt_y, name='sample_qc')


# In[ ]:


mt_a=mt_a.annotate_entries(AB = (mt_a.AD[1]/(mt_a.AD[0]+mt_a.AD[1])))
#mt_a = mt_a.filter_entries(mt_a.GT.is_het())
mt_y=mt_y.annotate_entries(AB = (mt_y.AD[1]/(mt_y.AD[0]+mt_y.AD[1])))
#mt_y = mt_y.filter_entries(mt_y.GT.is_het())


# In[ ]:


#This is for variant concordance, not call concordance
# Compute variants unique in each dataset
ds1 = mt_y.rows()   #yale
ds2 = mt_a.rows()  #acc

# unique in dataset1 -- filter out any variants in ds2
ds1_unique = ds1.anti_join(ds2)

# unique in dataset2 -- filter out any variants in ds1
ds2_unique = ds2.anti_join(ds1)

ds_shared = ds1.join(ds2)

print("There are "+ str(ds1_unique.count()) + " unique variants present in the Yale dataset.")
print("There are "+ str(ds2_unique.count()) + " unique variants present in the ACC dataset.")
print("There are "+ str(ds_shared.count()) + " shared variants.")


# In[ ]:


y_unique = mt_y.filter_rows(hl.is_defined(ds1_unique[mt_y.row_key]))
a_unique = mt_a.filter_rows(hl.is_defined(ds2_unique[mt_a.row_key]))
shared = mt_y.filter_rows(hl.is_defined(ds_shared[mt_y.row_key]))
shared_a = mt_a.filter_rows(hl.is_defined(ds_shared[mt_a.row_key]))


# In[ ]:


y_unique = hl.sample_qc(y_unique, name='sample_qc')
a_unique = hl.sample_qc(a_unique, name='sample_qc')
shared = hl.sample_qc(shared, name='sample_qc')
shared_a = hl.sample_qc(shared_a, name='sample_qc')


# In[ ]:


y_unique = hl.variant_qc(y_unique, name='variant_qc')
a_unique = hl.variant_qc(a_unique, name='variant_qc')
shared = hl.variant_qc(shared, name='variant_qc')
shared_a = hl.variant_qc(shared_a, name='variant_qc')


# In[ ]:


y_unique=y_unique.annotate_entries(AB = (y_unique.AD[1]/(y_unique.AD[0]+y_unique.AD[1])))
y_unique = y_unique.filter_entries(y_unique.GT.is_het())

a_unique=a_unique.annotate_entries(AB = (a_unique.AD[1]/(a_unique.AD[0]+a_unique.AD[1])))
a_unique = a_unique.filter_entries(a_unique.GT.is_het())

shared=shared.annotate_entries(AB = (shared.AD[1]/(shared.AD[0]+shared.AD[1])))
shared = shared.filter_entries(shared.GT.is_het())

shared_a=shared_a.annotate_entries(AB = (shared_a.AD[1]/(shared_a.AD[0]+shared_a.AD[1])))
shared_a = shared_a.filter_entries(shared_a.GT.is_het())


# In[ ]:


y_unique = y_unique.annotate_rows(mean_AB = hl.agg.mean(y_unique.AB))
a_unique = a_unique.annotate_rows(mean_AB = hl.agg.mean(a_unique.AB))
shared = shared.annotate_rows(mean_AB = hl.agg.mean(shared.AB))
shared_a = shared_a.annotate_rows(mean_AB = hl.agg.mean(shared_a.AB))


# In[ ]:


y_unique = y_unique.annotate_rows(mean_DP = hl.agg.mean(y_unique.DP))
a_unique = a_unique.annotate_rows(mean_DP = hl.agg.mean(a_unique.DP))
shared = shared.annotate_rows(mean_DP = hl.agg.mean(shared.DP))
shared_a = shared_a.annotate_rows(mean_DP = hl.agg.mean(shared_a.DP))


# In[ ]:


#Fitler PASS variants

y_unique_p = y_unique.filter_rows(hl.len(y_unique.filters) == 0)
a_unique_p = a_unique.filter_rows(hl.len(a_unique.filters) == 0)
shared_p = shared.filter_rows(hl.len(shared.filters) == 0)

y_unique_np = y_unique.filter_rows(hl.len(y_unique.filters) != 0)
a_unique_np = a_unique.filter_rows(hl.len(a_unique.filters) != 0)
shared_np = shared.filter_rows(hl.len(shared.filters) != 0)


# # Plot Metrics

# In[ ]:


#DP PER GENOTYPE
dp_hist = y_unique_p.aggregate_entries(hl.expr.aggregators.hist(y_unique_p.AB, 0, 1, 30))
p = hl.plot.histogram(dp_hist, title='Distribution of AB (Yale-Unique (PASS))')

fig = p
fig.xaxis.axis_label = 'Allele Balance'
fig.yaxis.axis_label = 'Frequency'

fig.title.align = "center"
fig.title.text_color = "black"
fig.title.text_font_size = "20px"

fig.x_range.start = 0
fig.x_range.end = 1

fig.background_fill_color = "white"

fig.xgrid.grid_line_color = None
fig.ygrid.grid_line_color = None

fig.xaxis.axis_label_text_font_size = "17px"
fig.xaxis.axis_label_text_color = "black"
fig.xaxis.axis_label_text_font_style = "normal"

fig.yaxis.axis_label_text_font_size = "17px"
fig.yaxis.axis_label_text_color = "black"
fig.yaxis.axis_label_text_font_style = "normal"

fig.outline_line_width = 0
fig.outline_line_color = "white"

fig.legend.border_line_width = 0
fig.legend.label_width= 0
fig.legend.label_height= 0
#fig.legend.location = "top_left"
fig.legend.visible = False
show(p)


# In[ ]:


#DP PER VARIANT
p = hl.plot.histogram(shared_np.mean_AB, title='Distribution of AB per variant (Shared (NO PASS))', bins = 30)

fig = p
fig.xaxis.axis_label = 'AB'
fig.yaxis.axis_label = 'Frequency'

fig.title.align = "center"
fig.title.text_color = "black"
fig.title.text_font_size = "20px"

fig.x_range.end = 1
#fig.x_range.start = 0.7
#fig.y_range.end = 5000

fig.background_fill_color = "white"

fig.xgrid.grid_line_color = None
fig.ygrid.grid_line_color = None

fig.xaxis.axis_label_text_font_size = "17px"
fig.xaxis.axis_label_text_color = "black"
fig.xaxis.axis_label_text_font_style = "normal"

fig.yaxis.axis_label_text_font_size = "17px"
fig.yaxis.axis_label_text_color = "black"
fig.yaxis.axis_label_text_font_style = "normal"

fig.outline_line_width = 0
fig.outline_line_color = "white"

fig.legend.border_line_width = 0
fig.legend.label_width= 0
fig.legend.label_height= 0
#fig.legend.location = "top_left"
fig.legend.visible = False
show(p)


# In[ ]:


var_plot = hl.plot.scatter(shared_a.sample_qc_new.r_het_hom_var, shared_a.sample_qc_new.n_snp, 
                           #colors=color_mapper,
                           title ='n_snps vs Het/Hom (Shared)')

#bokeh.models.mappers.ColorMapper or Mapping[str, bokeh.models.mappers.ColorMapper]

fig = var_plot
fig.xaxis.axis_label = 'Het/Hom Ratio'
fig.yaxis.axis_label = 'Number of variants'

#fig.y_range.start = 2.5
#fig.y_range.end = 3000
#fig.x_range.start = 25000
fig.x_range.end = 2.2

fig.title.align = "center"
fig.title.text_color = "black"
fig.title.text_font_size = "20px"

fig.background_fill_color = "white"

fig.xgrid.grid_line_color = None
fig.ygrid.grid_line_color = None

fig.xaxis.axis_label_text_font_size = "17px"
fig.xaxis.axis_label_text_color = "black"
fig.xaxis.axis_label_text_font_style = "normal"

fig.yaxis.axis_label_text_font_size = "17px"
fig.yaxis.axis_label_text_color = "black"
fig.yaxis.axis_label_text_font_style = "normal"

fig.outline_line_width = 0
fig.outline_line_color = "white"

'''fig.legend.border_line_width = 0
fig.legend.label_width= 0
fig.legend.label_height= 0
#fig.legend.location = "top_left"
#fig.legend.visible = False'''

hl.plot.show(fig)


# # Concordance between call sets

# In[ ]:


summary, samples, variants = hl.concordance(shared, shared_a)


# In[ ]:


left_homref_right_homvar = summary[2][4]
left_het_right_missing = summary[3][1]
left_het_right_something_else = sum(summary[3][:]) - summary[3][3]
total_concordant = summary[2][2] + summary[3][3] + summary[4][4]
total_discordant = sum([sum(s[2:]) for s in summary[2:]]) - total_concordant

samples.export("samples.tsv", delimiter = "\t")
variants.export("variants.tsv", delimiter = '\t')


# # Functional Annotation

# In[ ]:


##############################################################################
# Sort VEP transcript consequences array 
# annotate with desired VEP annotations as new row fields
# add gnomad 
# drop parent VEP struct (optional)
# write new matrix table 
##############################################################################

ref_data = '/gpfs/gibbs/pi/brueckner/hail_resources/combined_reference_data_grch38.ht'

################################
# specify inputs and outputs 
################################
#mt_input = "/gpfs/gibbs/pi/brueckner/ken_rotation/output/consortium_new/vep.mt" # specify the matrix table you wish to annotate. At minimum should already be annotated with VEP. 
#mt_output = "/gpfs/gibbs/pi/brueckner/ken_rotation/output/consortium_new/vep_r1.mt"

################################
#mt = hl.read_matrix_table(mt_input)
#print('Reads in the following matrix table: \n' + mt_input + '\n')
#print('count of the input dataframe:' + str(mt.count()) + '\n')
mt = y_unique_p

ref_data_ht = hl.read_table(ref_data)
print('Reads in the following hail table for annotation: \n' + ref_data + '\n')

################################
# add annotations (from vep)
################################

# Sort VEP array annotations 
##  sort transcript consequences 
mt = mt.annotate_rows(
    sortedTranscriptConsequences=get_expr_for_vep_sorted_transcript_consequences_array(vep_root=mt.vep)
)

# Annotate rows with index zero element of array sorted transcript annotations 
# Annotate rows with additional fields from VEP parent structure before dropping
# Annotate rows with gnomad AF (Assign 0.0 to variants without a gnomAD frequency)
mt = mt.annotate_rows(
    gene_symbol=hl.if_else(mt.sortedTranscriptConsequences.size() > 0, mt.sortedTranscriptConsequences[0].gene_symbol, hl.missing(hl.tstr)), 
    major_consequence=hl.if_else(mt.sortedTranscriptConsequences.size() > 0, mt.sortedTranscriptConsequences[0].major_consequence, hl.missing(hl.tstr)), 
    hgvs=hl.if_else(mt.sortedTranscriptConsequences.size() > 0, mt.sortedTranscriptConsequences[0].hgvs, hl.missing(hl.tstr)), 
    category=hl.if_else(mt.sortedTranscriptConsequences.size() > 0, mt.sortedTranscriptConsequences[0].category, hl.missing(hl.tstr)), 
    canonical=hl.if_else(mt.sortedTranscriptConsequences.size() > 0, mt.sortedTranscriptConsequences[0].canonical, -1),
    polyphen_pred=hl.if_else(mt.sortedTranscriptConsequences.size() > 0, mt.sortedTranscriptConsequences[0].polyphen_prediction, hl.missing(hl.tstr)),
    sift_pred=hl.if_else(mt.sortedTranscriptConsequences.size() > 0, mt.sortedTranscriptConsequences[0].sift_prediction, hl.missing(hl.tstr)),
    most_severe_consequence = mt.vep.most_severe_consequence,
    variant_class = mt.vep.variant_class,
    gene_id_array = mt.vep.transcript_consequences.gene_id,
    gene_symbol_array = mt.vep.transcript_consequences.gene_symbol,
    #gnomad_af = hl.if_else(hl.is_defined(ref_data_ht[mt.row_key]), ref_data_ht[mt.row_key].gnomad_exomes.AF, 0.0)
    gnomad_af =  hl.if_else(hl.is_defined(ref_data_ht[mt.row_key]),(hl.if_else(hl.is_defined(ref_data_ht[mt.row_key].gnomad_exomes.AF),ref_data_ht[mt.row_key].gnomad_exomes.AF, 0.0)), 0.0),
    meta_svm =  hl.if_else(hl.is_defined(ref_data_ht[mt.row_key]),(hl.if_else(hl.is_defined(ref_data_ht[mt.row_key].dbnsfp.MetaSVM_pred),ref_data_ht[mt.row_key].dbnsfp.MetaSVM_pred, 'na')), 'na'),
    cadd =  hl.if_else(hl.is_defined(ref_data_ht[mt.row_key]),(hl.if_else(hl.is_defined(ref_data_ht[mt.row_key].cadd.PHRED),ref_data_ht[mt.row_key].cadd.PHRED, 0.0)), 0.0)
)


# drop parent vep structure 
mt = mt.drop(mt.vep)

################################
# write new matrix table
################################

print('Writes out to the following matrix table: \n' + mt_output + '\n')
mt.write(mt_output, overwrite = True)


# In[ ]:


#mt_input = vep_r1c1_mt
#mt = hl.read_matrix_table(mt_input)
#Could not save the vep_r1. Hence use the variable from above

#print('Reads in the following matrix table: \n' + mt_input)
print('count of the input dataframe:' + str(mt.count()) + '\n')

# mt.describe()

# How many of each variant class? 
vep_variant_class = mt.variant_class
variant_class_dict = mt.aggregate_rows(hl.agg.counter(vep_variant_class))
print('Number of each variant class:')
pprint.pprint(variant_class_dict) 
print('\n')

# How many of each consequence type?
vep_most_severe_consequence = mt.most_severe_consequence
most_severe_term_dict = mt.aggregate_rows(hl.agg.counter(vep_most_severe_consequence))
print('Number of each consequence type:')
pprint.pprint(most_severe_term_dict)
print('\n')


# In[ ]:


#Filter for LOF

y_unique_lof_p = y_unique_p.filter_rows(hl.set(['frameshift_variant', 'stop_gained', "splice_acceptor_variant", "splice_donor_variant"]).contains(y_unique_p.most_severe_consequence))
a_unique_lof_p = a_unique_p.filter_rows(hl.set(['frameshift_variant', 'stop_gained', "splice_acceptor_variant", "splice_donor_variant"]).contains(a_unique_p.most_severe_consequence))


# In[ ]:


#Split the SNPs and Indels for analysis
y_unique_snp = y_unique.filter_rows(hl.set(["SNV"]).contains(y_unique.variant_class))
y_unique_indel = y_unique.filter_rows(hl.set(["insertion", "deletion"]).contains(y_unique.variant_class))

a_unique_snp = a_unique.filter_rows(hl.set(["SNV"]).contains(a_unique.variant_class))
a_unique_indel = a_unique.filter_rows(hl.set(["insertion", "deletion"]).contains(a_unique.variant_class))

shared_snp = shared.filter_rows(hl.set(["SNV"]).contains(shared.variant_class))
shared_indel = shared.filter_rows(hl.set(["insertion", "deletion"]).contains(shared.variant_class))

