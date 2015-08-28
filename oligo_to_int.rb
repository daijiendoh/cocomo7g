require './CMyRubyHelper.so'

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
	a_naseq=naseq.downcase.split("")
	base_deg = {
		'a' => 0,	'c' => 0,	'g' => 0,	't'=>0,
	        'w' => 1, 's' => 1, 'm' => 1, 'k' => 1, 'r' => 1, 'y' => 1, 
	        'b' => 2, 'd' => 2, 'h' => 2, 'v' => 2, 
	        'n' => 3
	}
	p base_deg
	p a_naseq
	a_naint=a_naseq.map{|x| base_deg[x]}
	return a_naint
end

oligo="CYCARGCAAAATGGMGGTTAAC"
p na_to_int(oligo)
p na_to_deg(oligo)