max_iterations = 0.0
max_iterations = ARGV[3].to_i if ARGV.size > 3
if ARGV.size < 3
  puts "Usage: provide 3 parameters"
  puts " - FASTA file"
  puts " - k-mer length (55 by default)"
  puts "- output file path"
  # puts "(additional parameter) max number of iterations"
  exit
end
# start of processing
require 'graphviz'
require './graph_utils'
#now reference_string contains reference genome
# puts reference_string

k = ARGV[1].to_i
k = 55 if k == 0
GraphUtils.clear_all

GraphUtils.set_kmer_size(k)
if ARGV[0].index(".fastq").nil?
  GraphUtils.read_fasta(ARGV[0]) do |read|
    GraphUtils.build_kmers_with_complement(read)
  end
else
  GraphUtils.read_fastq(ARGV[0]) do |read|
    GraphUtils.build_kmers_with_complement(read)
  end
end
#
GraphUtils.simplify!
r = 0
if ARGV[2] == '1'
  r = GraphUtils.remove_tails!
else
  r = GraphUtils.remove_tips!
end
f = File.new("results/_" + ARGV[0] + "_" + ARGV[2] + "_with_k#{k}.txt", 'w')
f.puts "Deleted: #{r}"
f.close

GraphUtils.simplify! if r > 0
#GraphUtils.print_to_file
#GraphUtils.write_to_console
#GraphUtils.write_to_fasta("results/_" + ARGV[0] + ".fasta")
#GraphUtils.write_to_dot("results/_" + ARGV[0] + ".dot")
GraphUtils.write_to_fasta_and_dot("results/_" + ARGV[0] + "_" + ARGV[2] + "_with_k#{k}.fasta", "results/_" + ARGV[0] + "_" + ARGV[2]  + "_with_k#{k}.dot")

GraphViz.parse( "results/_" + ARGV[0] + "_" + ARGV[2] + "_with_k#{k}.dot", :path => File.dirname(__FILE__) ).output(:svg => "results/_" + ARGV[0] + "_" + ARGV[2] + "_with_k#{k}.svg")