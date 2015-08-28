require 'bio'
require 'oj'
require './CMyRubyHelper.so'

gr_gbfa=Oj.load_file("../gr_gbfa.oj", :mode => :compat)
p gr_gbfa

repr_gbinf=Oj.load_file("../repr_gbinf.oj", :mode => :compat)
p repr_gbinf["337"]

ch_primer=[[0, 1886, "GYCARAATGGWGATTGGGACTT"], [1, 2484, "GGGAAKGCWKWTTCRWRKAGTG"]]
primer=ch_primer[1][2]
p CMyRubyHelper.clibFindPattern(repr_gbinf["337"][1].upcase, primer)