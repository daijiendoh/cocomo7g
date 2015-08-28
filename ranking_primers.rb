require 'oj'

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

## Method test
# ai_oligo=[0,0,0,1,2,0,0]
# p index_sum_deg(ai_oligo)
# p index_number_deg(ai_oligo)


# Main
comb_fasta_name="bunyaviridae"
short_tgt="bunya"

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