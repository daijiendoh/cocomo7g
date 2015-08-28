require 'fileutils'

p File.expand_path("./")
comb_fasta_name=File.basename(File.expand_path("./")).sub("cocomo7_","")
short_tgt=comb_fasta_name

p comb_fasta_name
p short_tgt