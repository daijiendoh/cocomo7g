require 'bio'
require 'oj'
require './CMyRubyHelper.so'

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
