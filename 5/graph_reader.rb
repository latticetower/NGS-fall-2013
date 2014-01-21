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
    GraphUtils.build_kmers(read)
    GraphUtils.build_kmers(GraphUtils.get_complement(read))
  end
else
  GraphUtils.read_fastq(ARGV[0]) do |read|
    GraphUtils.build_kmers(read)
    GraphUtils.build_kmers(GraphUtils.get_complement(read))
  end
end
#
GraphUtils.simplify!
#GraphUtils.print_to_file
#GraphUtils.write_to_console
#GraphUtils.write_to_fasta("results/_" + ARGV[0] + ".fasta")
#GraphUtils.write_to_dot("results/_" + ARGV[0] + ".dot")
GraphUtils.write_to_fasta_and_dot("results/_" + ARGV[0] + "_with_k#{k}.fasta", "results/_" + ARGV[0] + "_with_k#{k}.dot")

GraphViz.parse( "results/_" + ARGV[0] + "_with_k#{k}.dot", :path => File.dirname(__FILE__) ).output(:svg => "results/_" + ARGV[0] + "_with_k#{k}.svg")

=begin

node_hashes = []
edge_labels = []
deBruijn = GraphViz.new( :G, :type => :digraph )

node_label1 = genome_string.slice(0, k - 1)
node_hashes << node_label1
node1 = deBruijn.add_nodes("0")

(genome_string.length - (k - 1)).times do |i|
  edge_label = node_label1 + genome_string[k - 1 + i]# now here we have k mer
  node_label2 = edge_label.slice(1, k - 1)
  edge_labels << edge_label unless edge_labels.include? edge_label
  # connect nodes
  if node_hashes.include?(node_label2)
    node2 = deBruijn.find_node(node_hashes.index(node_label2).to_s)
  else 
    node_hashes << node_label2
    node2 = deBruijn.add_nodes(node_hashes.index(node_label2).to_s)
  end
  # node2 = deBruijn.find_node(node_label2)
  # node2 = deBruijn.add_nodes(node_label2) if node2.nil?
  
  if node1.neighbors.include? node2
    #do nothing here, not yet
  else
    deBruijn.add_edges(node1, node2, :label => edge_labels.index(edge_label).to_s)
  end
  
  
  # set node1
  node_label1 = node_label2
  node1 = node2
end

#hello = g.add_nodes( "Hello" )
#world = g.add_nodes( "World" )

# Create an edge between the two nodes
#g.add_edges( hello, world )



# Generate output image
deBruijn.output( :dot => "hello_world.dot" )
=end