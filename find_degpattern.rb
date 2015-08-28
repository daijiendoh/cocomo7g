require './CMyRubyHelper.so'

seq="KKKKKKKSACRRARWWWWWWWWSSSSSS"
seq2="ACAAATACRRRWWWWWWWWSSSSSS"
retarr=CMyRubyHelper.clibFindPattern(seq, "ACAAAAT")
p retarr
retarr2=CMyRubyHelper.clibFindPattern(seq, "TTTT")
p retarr2