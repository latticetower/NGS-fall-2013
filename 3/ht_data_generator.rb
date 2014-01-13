# ruby 1.8.7 compatible
# this codefile generates data hashes from sam files.
# it is used to be executed remotely.
#
# usage: ARGV[0] - reference genome filename
#        ARGV[1] - sam file with reads
#        ARGV[2] - output path prefix
# ARGV[0] : 'DH10B-K12.fasta'
# ARGV[1] : 'tst.fastq.sam'
# ARGV[2] : 'data/bwa_data.txt'
require './helpers/simple_sam_reader.rb'

max_iterations = 0.0
max_iterations = ARGV[3].to_i if ARGV.size > 3
if ARGV.size < 3
  puts "Usage: provide 3 parameters"
  puts " - reference file path"
  puts " - SAM file path"
  puts "- output file path"
  puts "(additional parameter) max number of iterations"
end

=begin
Tasks
(+) 1. coverage
- plot
- percentage of coverage
- average coverage
(+) 2. indels distribution
3. homopolymer line length vs indels count
4. nucleotides quality
- misinsertion
- mismatch
(+) 5. frequencies table
=end

# first, we need to read reference genome
@referencefile = File.new(ARGV[0])
reference_string = ""
while not @referencefile.eof do
  str = @referencefile.gets().strip
  reference_string << str if str.index('>').nil?
end
@referencefile.close

#now reference_string contains reference genome
# puts reference_string
puts "reference string processed"
# -----------
# subtask 1: coverage plot
# -----------
results_array1 = []

# subtask 2: indels
results_array2 = []
#------------
# subtask 5: frequencies table
#------------
# second, we process SAM file
results_array5 = {
  'A' => {'C' => 0, 'G' => 0, 'T' => 0, '_' => 0}, 
  'C' => {'A' => 0, 'G' => 0, 'T' => 0, '_' => 0},
  'G' => {'C' => 0, 'A' => 0, 'T' => 0, '_' => 0},
  'T' => {'C' => 0, 'G' => 0, 'A' => 0, '_' => 0},
  '_' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0}
}

total_length = 0
iter = 0
sam_reader_cycle(ARGV[1]) do |line|
  total_length = line.total_length
  puts "#{iter} * 100 iterations passed" if iter % 100 == 0 #debug
  break if max_iterations > 0 && max_iterations < iter
  # -----------
  # subtask 1
  # -----------
  if line.coverage > 2
    (line.coverage - 2).times do |i|
      if line.start_pos + 1 + i >= 0
        index = (line.start_pos + i ) / 1000
        results_array1[index] ||= [index, 0]
        results_array1[index][1] += 1
      end
    end
  end
  # --------------------
  # subtask 2: indels
  # --------------------
  indel_max_length = line.find_indels().max{|x, y| x[1] <=> y[1]} || [0,0]
  
  results_array2[indel_max_length[1]] ||= [indel_max_length[1], 0]
  results_array2[indel_max_length[1]][1] += 1
  
  #indel_lengths.each do |indel_length|
  #  results_array2[indel_length[1]] ||= [indel_length[1], 0]
  #  results_array2[indel_length[1]][1] += 1
  #end

  # -----------
  # subtask 5: frequency table
  # -----------
  a = line.find_substitutions(reference_string)
  results_array5.keys.each do |k1|
    results_array5[k1].keys.each do |k2|
      results_array5[k1][k2] += a[k1][k2]
    end
  end
  iter = iter + 1
end
# subtask 1 postprocess
results_array1.each_index do |i|
  results_array1[i] ||= [i, 0] #if results_array1[i]
end

results_array1.map!{|x| [x[0], x[1]/1000] }

# ----------------
# subtask 1 output
# ----------------
puts 'printing output'
@outputfile = File.new(ARGV[2], 'w')
@outputfile.puts("{:type=>'plot', :filename => 'coverage_plot.svg', :plot_data => { :name => 'Coverage plot', :x_name => 'position (1K nucleotides)', :y_name => 'Number of reads', :values => [" + 
    results_array1.map{|i| (!i.nil? ? "[#{i[0]}, #{i[1] || 0}, 1]" : "")}.join(', ') + "] }}")
@outputfile.puts("{:coverage_percentage =>  #{ (100.0*1000.0*(results_array1.count{|x| !x.nil? && x[1] > 0 } || 0)/total_length)} }")
@outputfile.puts("{:coverage_avg => #{(results_array1.map{|x| x[1] }.reduce(&:+)||0)/(results_array1.count || 1)} }")
@outputfile.close

# ----------------
# subtask 2 output
# ----------------
puts 'printing output for subtask 2'
@outputfile = File.new(ARGV[2], 'a')

# not so elegant, i know
@outputfile.print "{",
              ":type => 'plot', :filename => 'indels_distribution.svg', " + 
              ":plot_data => { :name => 'Indels distribution plot', ",
                ":x_name => 'Indel length', :y_name => 'Number of appearances', ",
                ":values => [" + 
                  results_array2.compact.map{|i| 
                    (!i.nil? ? "[#{i[0] || 0}, #{i[1] || 0}, 1]" : "")
                  }.join(', ') + 
                "] }" +
              "}"
@outputfile.puts ''

@outputfile.close

#-----------------------------
# subtask 5 output
#-----------------------------

puts 'printing output for subtask 5'
@outputfile = File.new(ARGV[2], 'a')

# not so elegant, i know
@outputfile.puts "{" + ['A', 'C', 'G', 'T', '_'].map { |x1| 
    ['A', 'C', 'G', 'T', '_'].select{|x2| 
      results_array5[x1].has_key?(x2) 
      }.map{|x2| "'#{x1}#{x2}' => #{results_array5[x1][x2]}"
  }.join(', ')}.join(', ') + "}"

@outputfile.close

#puts ARGV[1]

