require 'bio'
require 'csv'
require 'oj'
require 'fileutils'
require './CMyRubyHelper.so'

def import_primers_pos(p3_result_file)
	seq_id=""
	primer_no=0
	primer_direction=""
	primer_seq=""
	h_pseq_f={}
	h_pseq_r={}
	File.open(p3_result_file){|file|
		while line=file.gets
			#p line
			if line =~ /\ASEQUENCE_ID/ then
				seq_id=line.strip.sub("SEQUENCE_ID=","")
				p seq_id
			elsif line =~ /\APRIMER\_LEFT\_\d+\_SEQUENCE\=/ then
				primer_core=line.strip.gsub(/\APRIMER\_|\_SEQUENCE/,"")
				primer_data= primer_core.split(/\=/)
				primer_set= primer_data[0].split("_")
				pset_id=primer_set[1].to_i
				#p pset_id
				h_pseq_f[pset_id]=[primer_data[1]]
				#p h_pseq_f[pset_id]
			elsif line =~ /\APRIMER\_RIGHT\_\d+\_SEQUENCE\=/ then
				primer_core=line.strip.gsub(/\APRIMER\_|\_SEQUENCE/,"")
				primer_data= primer_core.split(/\=/)
				primer_set= primer_data[0].split("_")
				pset_id=primer_set[1].to_i
				h_pseq_r[pset_id]=[primer_data[1]]
			elsif line =~ /\APRIMER\_LEFT\_\d+\=/ then
				pset_id=line.strip.gsub(/\APRIMER\_LEFT\_|\=\d+\,\d+/,"").to_i
				pos_len=line.strip.gsub(/\APRIMER\_LEFT\_\d+\=/,"")
				pos_data=pos_len.strip.split(",")
				if h_pseq_f[pset_id] then
					h_pseq_f[pset_id] << pos_data[0].to_i
					h_pseq_f[pset_id] << pos_data[1].to_i
				else
					p "error in import"
				end
			elsif line =~ /\APRIMER\_RIGHT\_\d+\=/ then
				pset_id=line.strip.gsub(/\APRIMER\_RIGHT\_|\=\d+\,\d+/,"").to_i
				pos_len=line.strip.gsub(/\APRIMER\_RIGHT\_\d+\=/,"")
				pos_data=pos_len.strip.split(",")
				if h_pseq_r[pset_id] then
					h_pseq_r[pset_id] << pos_data[0].to_i
					h_pseq_r[pset_id] << pos_data[1].to_i
				else
					p "error in import"
				end
			else
				#p "No selection muched"
			end
		end
	}

	h_pset=Hash.new
	h_pset[seq_id]=[h_pseq_f,h_pseq_r]
	return h_pset
end
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
	  	position= feature.position
	  	# p position
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

def hash_to_na_uniqfasta(repr_gbid, repr_gbsp, repr_na_dir, uniq_fasta_name)
	repr_gbinf=Hash.new
	File.open(uniq_fasta_name,"w"){|file|
		i=1
		repr_gbid.each{|gid, naseq|
			file.puts(">#{i}")
			file.puts("#{naseq}")
			repr_gbinf[i]=[repr_gbsp[gid], naseq]
			i+=1			
		}
	}	
	return repr_gbinf
end

def max_group(h_grno)
	max_grno=Hash.new
	max_gr=h_grno.max{|x,y| x[1].length <=> y[1].length}

	return max_gr
end

def select_repr_gbfile(gbfiles, prot_dir, comb_fasta_name, repr_na_dir, uniq_fasta_name, limit_min_homol)
	h_vprot=Hash.new
	gbfiles.each_with_index{|gbfile,idx|
		h_vprot=gbfile_to_hash(gbfile)
	}
	hash_to_prfasta(h_vprot, prot_dir, comb_fasta_name)

	system("./usearch8 -makeudb_usearch ./#{prot_dir}/#{comb_fasta_name}.fa -output ./#{prot_dir}/#{comb_fasta_name}.udb")
	system("./usearch8 -usearch_global ./#{prot_dir}/#{comb_fasta_name}.fa -db ./#{prot_dir}/#{comb_fasta_name}.udb -id 0.98 -blast6out ./#{prot_dir}/prot_#{comb_fasta_name}_homol.b6 -maxaccepts 8 -maxrejects 256")

	bhomols= Dir.glob("#{prot_dir}/*.b6")
	acc_homol=Hash.new
	bhomols.each{|file|
		set_arr=Array.new
		if File.size(file) > 0 then
			pair_acc= file.sub("#{prot_dir}/","").sub(".b6","")
			CSV.foreach(file) do |row|
				set_arr << row[0].split("\t")
			end
		end
		acc_homol[pair_acc]=set_arr
	}
	# p acc_homol
	node_homol=Array.new
	acc_homol.each{|accs, homol|
		# 
		homol.each{|harr|
			len1=harr[7].to_i-harr[6].to_i+1
			len2=harr[9].to_i-harr[8].to_i+1
			node1=harr[0].split("|")[0]
			node2=harr[1].split("|")[0]
			if harr[2].to_f > limit_min_homol && len1 == len2 && node1 != node2 then			
				node_homol << [node1, node2]
			end
		}
	}
	# p node_homol
	mcl_abc_file="./#{prot_dir}/mcl_#{comb_fasta_name}.abc"
	mcl_cluster_file="./#{prot_dir}/mcl_#{comb_fasta_name}.clt"

	File.open(mcl_abc_file, "w"){|file|
		node_homol.each{|h|
			file.puts("#{h[0]} #{h[1]} 1")
		}
	}

	system("mcl #{mcl_abc_file} --abc -o #{mcl_cluster_file}")

	##  h_vprot[seq.entry_id][0] = seq.organism, h_vprot[seq.entry_id][1]=vprot, h_vprot[seq.entry_id][2]=naseq
	## repr_gbsp:  representation of species
	repr_gbid=Hash.new # representation of gbid
	repr_gbsp=Hash.new # representation of species

	File.open(mcl_cluster_file){|file|
		while line=file.gets
			line_gid=[]
			line_sp=[]
			aline= line.chomp.split(/\t/)
			reprid=aline[0]
			unless repr_gbsp[reprid] then
				repr_gbsp[reprid]=Array.new
			end
			aline.each{|gid|
				line_gid << gid
				line_sp << h_vprot[gid][0]
			}
			
			repr_gbid[reprid] = h_vprot[aline[0]][2]
			repr_gbsp[reprid] = [line_gid.sort, line_sp.uniq.sort]

		end
	}

	## repr_file="./#{repr_na_dir}/repr_#{comb_fasta_name}"
	# p repr_gbid
	# p repr_gbsp
	hash_to_nafasta(h_vprot, repr_na_dir, comb_fasta_name)
	repr_gbinf = hash_to_na_uniqfasta(repr_gbid, repr_gbsp, repr_na_dir, uniq_fasta_name)
	return repr_gbinf
end

def first_round_na_homology(repr_gbinf, repr_na_dir, comb_fasta_name, uniq_fasta_name, uniq_db_name, uniq_homol_name, min_na_length_for_design)

	system("makeblastdb -in #{uniq_fasta_name} -dbtype nucl -out #{uniq_db_name}")
	system("blastn -db #{uniq_db_name} -query #{uniq_fasta_name} -outfmt 6 -out #{uniq_homol_name}")

	# create homology cluster
	homol_sets=Hash.new
	cum_homol=Hash.new
	CSV.foreach(uniq_homol_name) do |row|
		aline= row[0].split("\t")
		unless cum_homol[[aline[0],aline[1]]] then
			cum_homol[[aline[0],aline[1]]] =[]
		end
		cum_homol[[aline[0],aline[1]]] << [aline[2],aline[3], aline[6], aline[7], aline[8], aline[9]]
	end
	sum_homol=Hash.new
	cum_homol.each{|set, homol|
		sum_h=0
		homol.each{|data|
			sum_h+=data[1].to_i
		}
		sum_homol[set]=sum_h
	}
	p sum_homol
	cum_homol.each{|set, data|
		len_query=repr_gbinf[set[0]][1].length
		len_tgt=repr_gbinf[set[1]][1].length
		# print len_query,",",len_tgt;puts
		if set[0].to_i != set[1].to_i && len_query/len_tgt.to_f > 0.7 && len_query/len_tgt.to_f < 1.3 && len_query > 500 && sum_homol[set] > len_query.to_f*0.8 then
			unless homol_sets[set] then
				homol_sets[set] = []
			end
			homol_sets[set] << sum_homol[set].to_i
		end
	}
		# length of homologous entries

	# p homol_sets
	mcl_abc_file="./#{repr_na_dir}/mcl_#{comb_fasta_name}.abc"
	mcl_cluster_file="./#{repr_na_dir}/mcl_#{comb_fasta_name}.clt"
	File.open(mcl_abc_file, "w"){|file|
		homol_sets.each{|set, h|
			file.puts("#{set[0]} #{set[1]} 1")
		}
	}

	system("mcl #{mcl_abc_file} --abc -o #{mcl_cluster_file}")

	max_gbno=repr_gbinf.length
	ori_gbnos=(1..max_gbno).to_a
	gr_gbfa=Hash.new
	gr_gblen=Hash.new
	File.open(mcl_cluster_file){|file|
		i=1
		while line=file.gets
			aline= line.chomp.split(/\t/)
			gr_gbfa[i]=aline
			gr_gblen[i]=aline.map{|m| [m,repr_gbinf[m.to_s][1].length]}.sort_by{|x| x[1]}
			i+=1
		end
	}

	return gr_gbfa
end

def import_primers_pos(p3_result_file)
	seq_id=""
	primer_no=0
	primer_direction=""
	primer_seq=""
	h_pseq_f={}
	h_pseq_r={}
	File.open(p3_result_file){|file|
		while line=file.gets
			#p line
			if line =~ /\ASEQUENCE_ID/ then
				seq_id=line.strip.sub("SEQUENCE_ID=","")
				p seq_id
			elsif line =~ /\APRIMER\_LEFT\_\d+\_SEQUENCE\=/ then
				primer_core=line.strip.gsub(/\APRIMER\_|\_SEQUENCE/,"")
				primer_data= primer_core.split(/\=/)
				primer_set= primer_data[0].split("_")
				pset_id=primer_set[1].to_i
				#p pset_id
				h_pseq_f[pset_id]=[primer_data[1]]
				#p h_pseq_f[pset_id]
			elsif line =~ /\APRIMER\_RIGHT\_\d+\_SEQUENCE\=/ then
				primer_core=line.strip.gsub(/\APRIMER\_|\_SEQUENCE/,"")
				primer_data= primer_core.split(/\=/)
				primer_set= primer_data[0].split("_")
				pset_id=primer_set[1].to_i
				h_pseq_r[pset_id]=[primer_data[1]]
			elsif line =~ /\APRIMER\_LEFT\_\d+\=/ then
				pset_id=line.strip.gsub(/\APRIMER\_LEFT\_|\=\d+\,\d+/,"").to_i
				pos_len=line.strip.gsub(/\APRIMER\_LEFT\_\d+\=/,"")
				pos_data=pos_len.strip.split(",")
				if h_pseq_f[pset_id] then
					h_pseq_f[pset_id] << pos_data[0].to_i
					h_pseq_f[pset_id] << pos_data[1].to_i
				else
					p "error in import"
				end
			elsif line =~ /\APRIMER\_RIGHT\_\d+\=/ then
				pset_id=line.strip.gsub(/\APRIMER\_RIGHT\_|\=\d+\,\d+/,"").to_i
				pos_len=line.strip.gsub(/\APRIMER\_RIGHT\_\d+\=/,"")
				pos_data=pos_len.strip.split(",")
				if h_pseq_r[pset_id] then
					h_pseq_r[pset_id] << pos_data[0].to_i
					h_pseq_r[pset_id] << pos_data[1].to_i
				else
					p "error in import"
				end
			else
				#p "No selection muched"
			end
		end
	}

	h_pset=Hash.new
	h_pset[seq_id]=[h_pseq_f,h_pseq_r]
	return h_pset
end

def import_primers_set(p3_result_file)
	seq_id=""
	primer_no=0
	primer_direction=""
	primer_seq=""
	h_pseq={}
	
	File.open(p3_result_file){|file|
		while line=file.gets
			#p line
			if line =~ /\ASEQUENCE_ID/ then
				seq_id=line.strip.sub("SEQUENCE_ID=","")
				# p seq_id
			elsif line =~ /\APRIMER\_LEFT\_\d+\_SEQUENCE\=/ then
				primer_core=line.strip.gsub(/\APRIMER\_|\_SEQUENCE/,"")
				primer_data= primer_core.split(/\=/)
				primer_set= primer_data[0].split("_")
				pset_id=primer_set[1].to_i
				#p pset_id
				h_pseq[pset_id]=[[primer_data[1]]]
				#p h_pseq_f[pset_id]
			elsif line =~ /\APRIMER\_RIGHT\_\d+\_SEQUENCE\=/ then
				primer_core=line.strip.gsub(/\APRIMER\_|\_SEQUENCE/,"")
				primer_data= primer_core.split(/\=/)
				primer_set= primer_data[0].split("_")
				pset_id=primer_set[1].to_i
				h_pseq[pset_id] << [primer_data[1]]
			elsif line =~ /\APRIMER\_LEFT\_\d+\=/ then
				pset_id=line.strip.gsub(/\APRIMER\_LEFT\_|\=\d+\,\d+/,"").to_i
				pos_len=line.strip.gsub(/\APRIMER\_LEFT\_\d+\=/,"")
				pos_data=pos_len.strip.split(",")
				if h_pseq[pset_id] then
					h_pseq[pset_id][0] << pos_data[0].to_i
					h_pseq[pset_id][0] << pos_data[1].to_i
				else
					p "error in import"
				end
			elsif line =~ /\APRIMER\_RIGHT\_\d+\=/ then
				pset_id=line.strip.gsub(/\APRIMER\_RIGHT\_|\=\d+\,\d+/,"").to_i
				pos_len=line.strip.gsub(/\APRIMER\_RIGHT\_\d+\=/,"")
				pos_data=pos_len.strip.split(",")
				if h_pseq[pset_id] then
					h_pseq[pset_id][1] << pos_data[0].to_i
					h_pseq[pset_id][1] << pos_data[1].to_i
				else
					p "error in import"
				end
			else
				#p "No selection muched"
			end
		end
	}

	h_pset=[seq_id, h_pseq]
	return h_pset
end

def na_to_int(naseq)
	a_naseq=naseq.downcase.split("")
	base_int = {
		'a' => 1,	'c' => 2,	'g' => 4,	't'=>8,
	        'w' => 9, 's' => 6, 'm' => 3, 'k' => 12, 'r' => 5, 'y' => 10, 
	        'b' => 14,
	        'd' => 13,
	        'h' => 11,
	        'v' => 7,
	        'n' => 15
	}
	# p base_int
	# p a_naseq
	a_naint=a_naseq.map{|x| base_int[x]}
	return a_naint
end

def na_to_deg(naseq)
	# p naseq
	a_naseq=naseq.downcase.split("")
	base_deg = {
		'a' => 0,	'c' => 0,	'g' => 0,	't'=>0,
	        'w' => 1, 's' => 1, 'm' => 1, 'k' => 1, 'r' => 1, 'y' => 1, 
	        'b' => 2, 'd' => 2, 'h' => 2, 'v' => 2, 
	        'n' => 3
	}
	# p base_deg
	# p a_naseq
	a_naint=a_naseq.map{|x| base_deg[x]}
	return a_naint
end

def index_deg_seth(ai_oligo)
	stindex=ai_oligo[-3]+ai_oligo[-2]
	return stindex
end

def index_sum_deg(ai_oligo)
	# p ai_oligo
	ai_oligo.inject {|sum, n| sum + n }
end

def index_number_deg(ai_oligo)
	ai_oligo.count {|item| item > 0 }
end

###### Main program

# Definitions
## gb-> hash

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
# Gbfile -> hash
## arguments: gbfile, prot_dir, comb_fasta_name, repr_na_dir

repr_gbinf=select_repr_gbfile(gbfiles, prot_dir, comb_fasta_name, repr_na_dir, uniq_fasta_name, limit_min_homol)
## Maximum of gbno
Oj.to_file("./repr_gbinf.oj", repr_gbinf, :mode => :compat)

# restart from here
repr_gbinf=Oj.load_file("./repr_gbinf.oj", :mode => :compat)
max_gbno=repr_gbinf.length
p max_gbno

# p repr_gbinf
## first round na homol
##method - first round homology: arguments: uniq_fasta_name, uniq_db_name, uniq_homol_name, 
gr_gbfa=first_round_na_homology(repr_gbinf, repr_na_dir, comb_fasta_name, uniq_fasta_name, uniq_db_name, uniq_homol_name, min_na_length_for_design)

p gr_gbfa
Oj.to_file("./gr_gbfa.oj", gr_gbfa, :mode => :compat)

gr_gbfa=Oj.load_file("./gr_gbfa.oj", :mode => :compat)
# check length
# gr_gbfa.each{|grno, faset|
# 	faset_len=faset.map{|x| repr_gbinf[x][1].length}
# 	p faset_len
# }

# Make Greedy max_group
gr_gbspc=Hash.new
gr_gbfa.each{|grno, grmem|
	gr_gbspc[grno]=grmem.map{|x| repr_gbinf[grno.to_s][0][1]}.uniq.sort.flatten
}

residual_grno=gr_gbspc.dup
greedy_repr_grno=[]
while max_group(residual_grno)[1].length > 0 do
	max_grmember= max_group(residual_grno)[1]
	# p max_group(residual_grno)
	greedy_repr_grno << max_group(residual_grno)[0]
	residual_grno.each{|grno, gr|
		residual_grno[grno]=gr-max_grmember
	}
end
# p greedy_repr_grno
greedy_gb_fa=Hash.new
greedy_repr_grno.each{|grno|
	greedy_gb_fa[grno]=gr_gbfa[grno]
}

p greedy_gb_fa
Oj.to_file("./gr_gbspc.oj", gr_gbspc, :mode => :compat)

greedy_gb_fa.each{|grno, grmem|
	File.open("./#{alignment_repr_dir}/#{grno}.fa","w"){|file|
		grmem.each{|gr|
			file.puts(">#{gr}\n")
			# p gr
			file.puts("#{repr_gbinf[gr][1].upcase}\n")
		}
	}
	system("mafft  ./alignment_repr/#{grno}.fa > ./alignment_repr/align_mafft_#{grno}.txt")
}

deg_base = {
	'a' => ['a'], 'c' => ['c'], 'g' => ['g'], 't'=>['t'],
  'w' => ['a', 't'], 's' => ['c', 'g'], 'm' => ['a', 'c'], 'k' => ['g', 't'], 'r' => ['a', 'g'],
        'y' => ['c', 't'],
        'b' => ['c', 'g', 't'],
        'd' => ['a', 'g', 't'],
        'h' => ['a', 'c', 't'],
        'v' => ['a', 'c', 'g'],
        'n' => ['a', 'c', 'g', 't']
}
base_deg=deg_base.invert
p base_deg

align_dir="./alignment_repr"
mafft_alines=Dir.glob("#{align_dir}/align_mafft_*.txt")
nseq_deg=Hash.new
mafft_alines.each{|mfile|
	p mfile
	s =File::stat(mfile)
	if  s.size > 0 then
		grn=mfile.sub("#{align_dir}/align_mafft_","").sub(".txt","")
		ff=Bio::FlatFile.auto(mfile)
		h_nseq=Hash.new
		ff.each do |entry|
			h_nseq[entry.definition]= entry.naseq
		end
		comb_seq=Hash.new
		h_nseq.each{|entry, seq|
			aseq=seq.split("")
			aseq.each_with_index{|n,idx|
				unless comb_seq[idx] then
					comb_seq[idx]=[]
				end
				unless n=="-" then
					comb_seq[idx] << n
				end
			}
		}
		seq_freq=Hash.new
		comb_seq.each{|i, n|
			if n.length >0 then
				seq_freq[i]=[n.length, n.inject(Hash.new(0)){|r, i| r[i] += 1; r}.sort_by{|key,val| -val}]
			end
		}
		len_degseq=seq_freq.length
		seq_range=(0..len_degseq-1)
		dg_seq=Array.new
		seq_range.each{|i|
			
			sumf=0
			ndata=seq_freq[i]
			bases=[]
			j=0
			while  sumf/ndata[0].to_f < 0.9 do
				bases << ndata[1][j][0]
				sumf += ndata[1][j][1]
				j+=1
			end
			if ndata[1].length > 2 then
				p i
				p seq_freq[i]
				p bases.sort!
				p base_deg[bases.sort]
			end
			if base_deg[bases.sort] then
				dg_seq<< base_deg[bases.sort]
			else
				p bases
				p "no_nucleotide"
			end

			# p dg_seq
		}
		# p dg_seq
		nseq_deg[grn]=dg_seq.join("")
		# p nseq_deg[grn].upcase
		# primer="aaagggaccagga".upcase
		# p primer
		# p CMyRubyHelper.clibFindPattern(nseq_deg[grn].upcase, primer)
	end
}
# p nseq_deg
 Oj.to_file("./nseq_deg.oj", nseq_deg, :mode => :compat)
nseq_deg=Oj.load_file("./nseq_deg.oj", :mode => :compat)
# p nseq_deg["27"]
align_dir="./alignment_repr"
primer_folder="./primers"
primer3_settings_file="primer3_settings.txt"
mafft_alines=Dir.glob("#{align_dir}/align_mafft_*.txt")

prpair0_inf=Hash.new
prpair_inf=Hash.new
mafft_alines.each{|mfile|
	grn=mfile.sub("#{align_dir}/align_mafft_","").sub(".txt","")
	s =File::stat(mfile)
	offset=0
	if  s.size > 0 then
		ff=Bio::FlatFile.auto(mfile)
		tgt_seq=""
		# tgt_seq_1=""
		# select last sequence
		i=0
		ff.each do |f|
			if i==0 then
				tgt_seq=f.naseq
			end
			i += 1
		end
		# if /\A\-+/.match(tgt_seq_0) then
		# 	p grn
		# 	offset= /\A\-+/.match(tgt_seq_0).to_s.length
		# 	tgt_seq_0.gsub!("-","")
		# else
		# 	offset=0
		# end
		# ngseq_len=tgt_seq_0.length
		# st_target=offset+1
		# end_target=offset+ngseq_len
		# p [st_target, end_target]
		# tgt_seq=tgt_seq_1[st_target-1..end_target-1]
		# p tgt_seq.upcase
		if tgt_seq.length>200 then
			primer3_settings_file="primer3_settings.txt"
			
			p3_seq_file="#{primer_folder}/target_seq"
			p3_result_file="#{primer_folder}/p_primers.out"
			h_pset=Hash.new
			h_pset_f={} 
			h_pset_r={}
			# Contract primers and import results as primer set array 

			File.open(p3_seq_file,"w"){|file|
				file.puts("SEQUENCE_ID=#{grn}")
				file.puts("SEQUENCE_TEMPLATE=#{tgt_seq.gsub("-","n").upcase}")
				file.puts("=")
			}
			# File.open(p3_seq_file, "r") do |file|
			# 	file.each do |line|
			# 		p line
			# 	end
			# end 

			system("primer3_core -p3_settings_file=./#{primer3_settings_file} -output=./#{p3_result_file}  < #{p3_seq_file}")
			# h_pset=import_primers_pos(p3_result_file)
			# p h_pset
			h_pset=import_primers_set(p3_result_file)
			grno=h_pset[0].to_s
			unless prpair_inf[grno] then
				prpair_inf[grno]=Hash.new
			end
			unless prpair0_inf[grno] then
				prpair0_inf[grno]=Hash.new
			end
			# p h_pset[1]
			h_pset[1].each{|prno, prdata|
				# print "grno",grno;puts
				# print "prno= " ,prno;puts
				unless prpair_inf[grno][prno] then
					prpair_inf[grno][prno]=[]
				end
				unless prpair0_inf[grno][prno] then
					prpair0_inf[grno][prno]=[]
				end
				prpair0_inf[grno][prno]=prdata
				# prdata [["GGGTAAGGAAGCAGTCAATCAC", 1202, 22], ["ACACCCCTAATCTACCACCCTT", 1552, 22]]
				# print "prdata ", prdata; puts
				prdata.each{|primer|					
					# print "primer = ", primer;puts

					if nseq_deg[grno] then
						dseq= nseq_deg[grno].upcase
						# p dseq
						primer_seq=primer[0].upcase
						# p primer_seq
						pr_homol_pos=CMyRubyHelper.clibFindPattern(dseq, primer_seq)
						if pr_homol_pos.length==1 then
							prpair_inf[grno][prno] << pr_homol_pos[0]
						elsif pr_homol_pos.length > 1 then
							prpair_inf[grno][prno] << ["too_frequent"]
						else
							prpair_inf[grno][prno] << ["no_data"]
						end
					end
				}
			}
		end
	end
}
Oj.to_file("./prpair_inf.oj", prpair_inf, :mode => :compat)
Oj.to_file("./prpair0_inf.oj", prpair0_inf, :mode => :compat)
prpair_inf=Oj.load_file("./prpair_inf.oj", :mode => :compat)
prpair2_inf=Hash.new
prpair_inf.each{|grno, pset|
	unless prpair2_inf[grno] then
		prpair2_inf[grno]=Hash.new
	end
	pset.each{|prno, prpair|
		if prpair[0].length == prpair[1].length && prpair[0].length==3 then
			# p prpair
			p prno
			if na_to_deg(prpair[0][2])[-1]==0 && na_to_deg(prpair[1][2])[-1]==0 then
				unless prpair2_inf[grno][prno] then
					prpair2_inf[grno][prno]=Hash.new
				end
				# p prno
				prpair2_inf[grno][prno]= prpair
				print na_to_deg(prpair[0][2]), na_to_deg(prpair[1][2]);puts
			end
		end
	}
}
Oj.to_file("./prpair2_inf.oj", prpair2_inf, :mode => :compat)
p prpair2_inf

prpair2_inf=Oj.load_file("./prpair2_inf.oj", :mode => :compat)
prpair0_inf=Oj.load_file("./prpair0_inf.oj", :mode => :compat)

primer2_index=Hash.new
prpair2_inf.each{|grno, pset|
	unless primer2_index[grno] then
		primer2_index[grno]=Hash.new
	end
	pset.each{|prno, prpair|
		# print "prno = ", prno;puts
		unless primer2_index[grno][prno] then
			primer2_index[grno][prno]=Hash.new
		end
		pidx=0
		ideg_f= na_to_deg(prpair[0][2])
		ideg_r=na_to_deg(prpair[1][2])
		st_idx= (index_deg_seth(ideg_f)+ index_deg_seth(ideg_r))*100
		nd_idx=(index_number_deg(ideg_f)+index_number_deg(ideg_r))*5
		sumd_idx=(index_sum_deg(ideg_f)+index_sum_deg(ideg_r))
		# p st_idx
		# p nd_idx
		# p sumd_idx
		primer_index=st_idx+nd_idx+sumd_idx
		primer2_index[grno][prno]=primer_index
	}
}
repr3_primers=Hash.new
primer2_index.each{|grno, pridx|
	p grno
	unless repr3_primers[grno] then
		repr3_primers[grno]=Hash.new
	end
	repr_primer_pair=pridx.sort{|a,b| a[1]<=>b[1]}[0,3]
	repr_primer_pair.each{|pr|
		print pr[0],",",pr[1];puts
		# unless repr3_primers[grno][pr[0]] then
		# 	repr3_primers[grno][pr[0]]=Hash.new
		# end
		p prpair2_inf[grno][pr[0].to_s]
		repr3_primers[grno][pr[0]] = prpair2_inf[grno][pr[0].to_s]
	}
}

prpair_inf=Oj.load_file("./prpair_inf.oj", :mode => :compat)

uniq_seq=Hash.new
repr3_primers.each{|grno, prdata|
	p grno
	prdata.each{|prno, prset|
		p prno
		p prset
		unless uniq_seq[prset[0][2]] then
			uniq_seq[prset[0][2]]=Array.new
		end
		unless uniq_seq[prset[1][2]] then
			uniq_seq[prset[1][2]]=Array.new
		end
		prset0= prpair0_inf[grno][prno]
		uniq_seq[prset[0][2]] << [grno, prno, "f", prset0[0][0]]
		uniq_seq[prset[1][2]] << [grno, prno, "r", prset0[1][0]]
	}
}
Oj.to_file("./uniq_seq.oj", uniq_seq, :mode => :compat)
# Order_file
File.open("Order_#{short_tgt}.csv","w"){|file|
	uniq_seq.each{|prseq, pdata|
		file.puts("#{short_tgt}_#{pdata[0][0]}_#{pdata[0][1]}_#{pdata[0][2]}: #{prseq}")
	}
}
# p uniq_seq
pcr_prorder=Hash.new
repr3_primers.each{|grno, prdata|
	p grno
	unless pcr_prorder[grno] then
		pcr_prorder[grno]=Hash.new
	end
	prdata.each{|prno, prset|
		p prno
		p prset
		prset0= prpair0_inf[grno][prno]
		order_prf=uniq_seq[prset[0][2]][0][0,3].join("_")
		order_prr=uniq_seq[prset[1][2]][0][0,3].join("_")
		pcr_prorder[grno][prno] =[grno, prno, order_prf, order_prr,prset0[1][1].to_i-prset0[0][1].to_i, prset0[0][1].to_i,prset0[1][1].to_i]
	}
}
Oj.to_file("./pcr_prorder.oj", pcr_prorder, :mode => :compat)

File.open("PCR_pair_#{short_tgt}.csv","w"){|file|
	file.puts("group_no, prset_no, forword_order, reverse_order, pcr_size, position_f, position_r")
	pcr_prorder.each{|grno, prdata|
		prdata.each{|prno, ppair|
			file.puts("#{ppair.join(",")}")
		}
		
	}
}

gr_gbspc=Oj.load_file("./gr_gbspc.oj", :mode => :compat)

File.open("Group_species_#{short_tgt}.csv","w"){|file|
	gr_gbspc.each{|grno, spcdata|
		file.puts("#{grno},#{spcdata.join(",")}")
	}
}