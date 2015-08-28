require 'bio'

# read the sequence entry by entry through the files listed in ARGV.

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
	  	if feature.feature == 'CDS'
	  		f_hash= feature.to_hash
	  		# p f_hash["translation"]
	  		vprot=vprot.concat(f_hash["translation"][0]) 
	  	end
	  }
	  if vprot.length > 10 then 
	  	unless h_vprot[seq.entry_id] then
	  		h_vprot[seq.entry_id]=[]
	  	end
		  h_vprot[seq.entry_id][0] = seq.organism
		  h_vprot[seq.entry_id][1]=vprot
		  h_vprot[seq.entry_id][2]=naseq
		end
	end
	return h_vprot
end

def hash_to_prfasta(h_vprot, prot_dir, comb_fasta_name)
	File.open("./#{prot_dir}/#{comb_fasta_name}.fa","w"){|file|
		h_vprot.each{|k, v|
			file.puts(">#{k}|#{v[0]}")
			file.puts("#{v[1]}")
		}
	}	
end

def hash_to_nafasta(h_vprot, repr_na_dir, comb_fasta_name)
	File.open("./#{repr_na_dir}/#{comb_fasta_name}.fa","w"){|file|
		h_vprot.each{|k, v|
			file.puts(">#{k}|#{v[0]}")
			file.puts("#{v[2]}")
		}
	}	
end
######
gbfiles=Dir.glob("./virus_gbfile/*.gb")
prot_dir="prot_fasta"
repr_na_dir="repr_na"
comb_fasta_name="bunyaviridae"
h_vprot=Hash.new
gbfiles.each_with_index{|gbfile,idx|
	p idx
	p gbfile
	if idx==0 then
		h_vprot=gbfile_to_hash(gbfile)
	else
		h_vprot=h_vprot.merge(gbfile_to_hash(gbfile))
	end
	p h_vprot.length
}
hash_to_prfasta(h_vprot, prot_dir, comb_fasta_name)
hash_to_nafasta(h_vprot, repr_na_dir, comb_fasta_name)