require 'bio'
require 'fileutils'

def gbfile_to_hash(gbfile)
	sequences = Bio::FlatFile.auto(gbfile)
	core_name=gbfile.split("/")[2].sub(".gb","")
	h_vprot=Hash.new
	sequences.each do |seq|
	  # do stuff with the sequence entry
	  seq.entry_id
	  naseq=seq.naseq
	  # puts seq.definition
	  vprot=""
	  seq.features.each{|feature|
	  	position= feature.position
	  	p position
	  	cds_naseq= seq.naseq.splicing(position)
	  	if feature.feature == 'CDS'
	  		f_hash= feature.to_hash
	  		# p f_hash
	  		
	  		if f_hash["translation"] && f_hash["protein_id"] then
	  			protein_id= f_hash["protein_id"][0]
	  			h_vprot[protein_id]=[seq.organism, f_hash["translation"][0], cds_naseq] 
	  		end
	  	end
	  }

	end
	return h_vprot
end

def hash_to_prfasta(h_vprot, prot_dir, comb_fasta_name)
	File.open("./#{prot_dir}/#{comb_fasta_name}.fa","w"){|file|
		h_vprot.each{|k, v|
			file.puts(">#{k}|#{v[0]}")
			file.puts("#{v[1]}")
			# file.puts("#{v[2]}")
		}
	}	
end


##########################
prot_dir="prot_fasta"
repr_na_dir="repr_na"
div_na_dir="div_na"
alignment_repr_dir="alignment_repr"

p File.expand_path("./")
comb_fasta_name=File.basename(File.expand_path("./")).sub("cocomo7_","")
short_tgt=comb_fasta_name


## hash -> fasta


# homol analysis
limit_primer=22
limit_min_homol=98


mcl_abc_file="./#{prot_dir}/mcl_#{comb_fasta_name}.abc"
mcl_cluster_file="./#{prot_dir}/mcl_#{comb_fasta_name}.clt"

repr_file="./#{repr_na_dir}/repr_#{comb_fasta_name}"

min_na_length_for_design=150
uniq_fasta_name="./#{repr_na_dir}/#{comb_fasta_name}_uniq.fasta"
uniq_db_name="./#{repr_na_dir}/#{comb_fasta_name}_uniq.db"
uniq_homol_name="./#{repr_na_dir}/#{comb_fasta_name}_homol.b6"
uniq_fasta_div_file="./#{repr_na_dir}/#{comb_fasta_name}_uniq_div.fa"
uniq_div_db="./#{repr_na_dir}/#{comb_fasta_name}_uniq_div.db"
uniq_div_homol_name="./#{repr_na_dir}/#{comb_fasta_name}_uniq_div_homol.b6"
sel_div_core="./#{div_na_dir}/#{comb_fasta_name}_div"

mcl_div_abc_file="./#{repr_na_dir}/mcl_#{comb_fasta_name}_div.abc"
mcl_div_cluster_file="./#{repr_na_dir}/mcl_#{comb_fasta_name}_div.clt"
# run program

gbfiles=Dir.glob("./virus_gbfile/*.gb")
FileUtils.rm_r Dir.glob("./#{prot_dir}/*")
FileUtils.rm_r Dir.glob("./#{repr_na_dir}/*")
FileUtils.rm_r Dir.glob("./#{alignment_repr_dir}/*")

h_vprot=Hash.new
gbfiles.each_with_index{|gbfile,idx|
	h_vprot=gbfile_to_hash(gbfile)
}
hash_to_prfasta(h_vprot, prot_dir, comb_fasta_name)
# p h_vprot
