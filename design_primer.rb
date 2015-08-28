require 'bio'
require 'oj'
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
	p base_int
	p a_naseq
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
################################################
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