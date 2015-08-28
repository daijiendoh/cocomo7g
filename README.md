# cocomo7g
Degenerate Primer Desgin Program from Virus Genome

## Description

Degenerate Primer Desgin Program from Virus Genome

## Requirement

OS: Ubuntu > 14.04, 64bit
ruby > 2.2
g++
mcl
mafft
graphviz
gems: bio, oj, fileutils

## Install

setup ruby environment according to https://gorails.com/setup/ubuntu/15.04.
cd ~
sudo apt-get install g++ mcl mafft graphviz
gem install bio
gem install oj
gem install fileutils
git clone 

## Usage

1 Search viral genome on GenBank with a search word "(virus name) complete genome"
2 Download results with GenBank-format into the folder virus_gbfile.
3 cd ~/cocomo7g
4 ruby cocomo7.rb
5 get results in the created several txt files

## Licence

[MIT](https://github.com/tcnksm/tool/blob/master/LICENCE)

## Author

Daiji Endoh, Rakuno Gakuen University
dendoh@rakuno.ac.jp
