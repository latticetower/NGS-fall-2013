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
output_file_prefix = ""
output_file_prefix = ARGV[4] if ARGV.size > 4
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


# ------------
# nucleotides quality
# - misinsertion
# - mismatch
# ------------


# subtask 3: indels per homopolymer sequence length

# first, we need to read reference genome

# second, we process SAM file
results_hash3 = {}
# subtask4: hash for mismatches
results_hash4 = {:insertions => {}, :mismatches => {} }
# puts ARGV[1]
total_length = 0
iter = 0
#homopolymer_sequences = SAMString::find_homopolymer_sequences(reference_string).select{|x| x > 1}
puts "homopolymer sequences found"
# puts homopolymer_sequences.size

sam_reader_cycle(ARGV[1]) do |line|
  total_length = line.total_length
  puts "#{iter} iterations passed" if iter % 100 == 0
  break if max_iterations > 0 && iter > max_iterations 
  indels = line.find_indels()
  insertions_and_mismatches_array = line.find_insertions_and_mismatches
  last_char = 0
  homopolymer_length = 0
  # puts "in sam #{line.sequence.size}"
  # homopolymers = line.sequence.scan(/(a+|c+|g+|t+)/)
  sequence_start = line.start_pos
  sequence_length = line.sequence.size
  last_char = reference_string[sequence_start]
  while (sequence_start > 0 && last_char == reference_string[sequence_start - 1]) do 
    sequence_start -= 1
    sequence_length += 1
  end
  last_char = reference_string[sequence_start + sequence_length]
  while(sequence_start + sequence_length < total_length && last_char == reference_string[sequence_start + sequence_length + 1]) do
    sequence_length += 1
  end
  
  sequence_length.times do |index|
  #puts index
    
    if last_char != reference_string[sequence_start + index]
      if homopolymer_length > 0
        c = indels.select{|k, v| k < index && k + v >= index - homopolymer_length}.count
        if c > 0
          results_hash3[homopolymer_length] ||= 0
          results_hash3[homopolymer_length] += c
          
        end
      end
      homopolymer_length = 0
    end
    homopolymer_length += 1
    
    last_char = reference_string[sequence_start + index]
    # ------------------
    # let's do subtask 4
    # ------------------
    if sequence_start + index >= line.start_pos && sequence_start + index < line.start_pos + line.sequence.length
      real_index = sequence_start + index - line.start_pos
      if insertions_and_mismatches_array[real_index] == 'I'
        results_hash4[:insertions][line.qualities[real_index]] ||= 0
        results_hash4[:insertions][line.qualities[real_index]] += 1
      end
      if insertions_and_mismatches_array[real_index] == 'M'
        results_hash4[:mismatches][line.qualities[real_index]] ||= 0
        results_hash4[:mismatches][line.qualities[real_index]] += 1
      end
      
    end
    ##end_of_subtask4
  end
  
  
  iter = iter + 1
end
puts results_hash3.to_a.compact.map{|x| '[' + x.join(' , ') + ']'}.join ', '
puts 'printing indels to output'
@outputfile = File.new(ARGV[2], 'a')

# not so elegant, i know
@outputfile.print "{",
              ":type=>'plot', :filename => '#{output_file_prefix}indels_distribution_homopolymer.svg', " + 
              ":plot_data => {" + 
                ":name => 'Indels distribution per homopolymer sequence length ', ",
                ":x_name => 'homopolymer sequence length', :y_name => 'indels sequences count', ",
                ":values => [" + 
                  results_hash3.select{|k,v| v > 0}.to_a.compact.sort{|x,y| x[0] <=> y[0]}.map{|i| (!i.nil? ? "[#{i[0]}, #{i[1] || 0}, 1]" : "")}.join(', ') + 
                "] " +
              "}}"
@outputfile.puts ''
@outputfile.close

# subtask 4

@outputfile = File.new(ARGV[2], 'a')

# print part 1
@outputfile.print "{",
              ":type=>'plot', :filename => '#{output_file_prefix}qualities_insertion.svg', " + 
              ":plot_data => {" + 
                ":name => 'Qualities distribution for inserted nucleotides', ",
                ":x_name => 'Qualities (phred values)', :y_name => 'Number of inserted nucleotides', ",
                ":values => [" + 
                  results_hash4[:insertions].to_a.compact.sort{|x, y| x[0].ord <=> y[0].ord}.map{|i| (!i.nil? ? "[#{i[0].ord}, #{i[1] || 0}, 1]" : "")}.join(', ') + 
                "], :text_x => ["+ results_hash4[:insertions].to_a.compact.sort{|x, y| x[0].ord <=> y[0].ord}.map{|i| (i.nil? ? "": "%Q{#{i[0]}}")}.join(', ') +"] " +
              "}}"
@outputfile.puts ''
#print part 2
# print part 1
@outputfile.print "{",
              ":type=>'plot', :filename => '#{output_file_prefix}qualities_mismatch.svg', " + 
              ":plot_data => {" + 
                ":name => 'Qualities distribution for mismatched nucleotides', ",
                ":x_name => 'Qualities (phred values)', :y_name => 'Number of mismatched nucleotides', ",
                ":values => [" + 
                  results_hash4[:mismatches].to_a.compact.sort{|x, y| x[0].ord <=> y[0].ord}.map{|i| (!i.nil? ? "[#{i[0].ord}, #{i[1] || 0}, 1]" : "")}.join(', ') + 
                "], " +
                ":text_x => ["+ results_hash4[:mismatches].to_a.compact.sort{|x, y| x[0].ord <=> y[0].ord}.map{|i| (i.nil? ? "": "%Q{#{i[0]}}")}.join(', ') +"] "+
              "}}"
@outputfile.puts ''
@outputfile.close


#sam_reader_cycle('tst.fastq.sam') do |sam_line|
#  puts sam_line
#end