#!~/bin/ruby

#------------
# 
#------------

class SAMString
  attr_accessor :data_hash
  attr_accessor :coverage
  attr_accessor :start_pos, :end_pos, :total_length
  attr_reader :sequence
  attr_reader :qualities
  
  def to_s
    data_hash.to_s
  end
  
  def initialize(arg)
    @data_hash = arg
    # puts arg.to_s
    @coverage = (arg[:seq]||"").length
    #\*|([0-9]+)[MIDNSHPX=]
    aaa = (arg[:cigar]||"").scan(/\*|([0-9]+)[DN]/)
    aaa = aaa.flatten.map{|x| x || 0}.map{|x| x.to_i} 
    @coverage += aaa.inject(&:+)|| 0
    @qualities = arg[:qual] || ""
    @start_pos = arg[:pos]
    @end_pos   = arg[:pos] + @coverage - 1
    
    @total_length = arg[:total_length]
    @sequence = arg[:seq] || ""
  end
  
  def find_substitutions(reference_string)
    aaa = data_hash[:cigar].scan(/\*|[0-9]+[MIDNSHPX=]/)
    current_index = 0
    ref_index = 0
    result_hash = {
      'A' => {'G' => 0, 'C' => 0, 'T' => 0, '_' => 0},
      'C' => {'G' => 0, 'A' => 0, 'T' => 0, '_' => 0},
      'G' => {'A' => 0, 'C' => 0, 'T' => 0, '_' => 0},
      'T' => {'G' => 0, 'C' => 0, 'A' => 0, '_' => 0},
      '_' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0}
      }
    aaa.each do |subst|
      number_of_positions = subst.chop.to_i
      case subst[-1].chr
      when 'M'
        number_of_positions.times do |i|
          next if reference_string.size <= start_pos - 1 + ref_index + i
          if reference_string[start_pos - 1 + ref_index + i].chr.upcase != sequence[current_index + i].chr.upcase && sequence[current_index + i].chr != '=' && sequence[current_index + i].chr != '*'
            letter1 = reference_string[start_pos - 1 + ref_index + i].chr.upcase
            letter2 = sequence[current_index + i].chr.upcase
            
            result_hash[letter1][letter2] += 1 if result_hash[letter1].has_key? letter2
            #result_hash[letter2][letter1] += 1 if result_hash[letter1].has_key? letter1
            
          end
        end
        current_index += number_of_positions
        ref_index += number_of_positions
      when 'I'
        number_of_positions.times do |i|
          letter1 = '_'
          letter2 = sequence[current_index + i].chr.upcase
          result_hash[letter1][letter2] += 1 if result_hash[letter1].has_key? letter2
        end
        current_index += number_of_positions
      when 'D'
        number_of_positions.times do |i|
          if reference_string[start_pos - 1 + ref_index + i].chr.upcase != sequence[current_index + i].chr.upcase && sequence[current_index + i].chr != '=' && sequence[current_index + i].chr != '*'
            letter1 = reference_string[start_pos - 1 + ref_index + i].chr.upcase
            letter2 = '_'
            result_hash[letter1][letter2] += 1 if result_hash[letter1].has_key? letter2
          end
        end
        ref_index += number_of_positions
      when 'N'
        number_of_positions.times do |i|
          letter1 = '_'
          letter2 = sequence[current_index + i].chr.upcase
          result_hash[letter1][letter2] += 1 if result_hash[letter1].has_key? letter2
        end
        current_index += number_of_positions
      when 'S' # do we take into consideration soft clipping or not?
        current_index += number_of_positions
        ref_index += number_of_positions
      when 'H'
        ref_index += number_of_positions
      when 'P'
        ref_index += number_of_positions
      when '='
        current_index += number_of_positions
        ref_index += number_of_positions
      when 'X'
        number_of_positions.times do |i|
          if reference_string[start_pos - 1 + ref_index + i].chr.upcase != sequence[current_index + i].chr.upcase && sequence[current_index + i].chr != '=' && sequence[current_index + i].chr != '*'
            letter1 = reference_string[start_pos - 1 + ref_index + i].chr.upcase
            letter2 = sequence[current_index + i].chr.upcase
            
            result_hash[letter1][letter2] += 1 if result_hash[letter1].has_key? letter2
            #result_hash[letter2][letter1] += 1 if result_hash[letter1].has_key? letter1
          end
        end
        current_index += number_of_positions
        ref_index += number_of_positions
      end
    end
    result_hash
  end
  
  # the following method finds indels (I,D,N sequences) in string, and their lengths. 
  # returns array of indel lengths for the given string
  def find_indels()
    aaa = data_hash[:cigar].scan(/\*|[0-9]+[MIDNSHPX=]/)
    current_index = 0
    ref_index = 0
    result_array = {}
    last_one_was_indel = false
    last_indel_length = 0
    last_indel_start = 0
    aaa.each do |subst|
      number_of_positions = subst.chop.to_i
      if ['I', 'D', 'N'].include? subst[-1].chr
        last_indel_length += number_of_positions
        last_one_was_indel = true
      else
        if last_one_was_indel 
          result_array[last_indel_start]  = last_indel_length
        end
        last_one_was_indel = false
        last_indel_length = 0
      end #end of condition
      
      result_array[last_indel_start] = last_indel_length if last_one_was_indel and last_indel_length > 0
      last_indel_start += number_of_positions
    end #end of each
    result_array
  end
  
  def find_insertions_and_mismatches()
    aaa = data_hash[:cigar].scan(/\*|[0-9]+[MIDNSHPX=]/)
    current_index = 0
    ref_index = 0
    aaa.map{|x| (1..x.chop.to_i).map{ x[-1].chr.upcase }}.flatten
  end
  
  # the following method finds homopolymer sequencies in 
  # we take into consideration only indels in specified positions
  # returns array of indel lengths for the given string

  def self.find_homopolymer_sequences(reference_string)
    homopolymers_array = reference_string.scan(/(A+|C+|G+|T+)/)
    
    result_array = {}
    start_position = 0
    homopolymers_array.each do |homopolymer|
      end_position = start_position + homopolymer.size
      result_array[end_position] = homopolymer.size if homopolymer.size > 1
      start_position = end_position
    end
    
    result_array
  end

  # first parameter: homopolymer_sequences = find_homopolymer_sequences(reference_string)
  # second parameter: indels = find_indels(reference_string)
  #

  def find_in_homopolymers(homopolymer_sequences, indels)
    seq = homopolymer_sequences.select{|k, v| k - v <= end_pos && k > start_pos}
    return [] if seq.nil?
    result_array = []
    seq.each do |key, length|
      c = indels.count{|indel|  indel[0] + indel[1] >= key - length && indel[0] < key } 
      result_array << [ 
                length, c
                      ] if c > 0
    end
    result_array
  end
  
  def find_in_homopolymers2(homopolymer_sequences, indels)
    seq = homopolymer_sequences.select{|hseq| hseq[0] <= end_pos && hseq[0] + hseq[1] >= start_pos}
    return [] if seq.nil?
    result_array = []
    seq.each do |sequence|
      result_array << [ 
                sequence[1], indels.count{|indel|  indel[0] + indel[1] >= sequence[0] && indel[0] <= sequence[0] + sequence[1] } 
                      ]
    end
    result_array
  end


  
  def tlen
    data_hash[:tlen]
  end
  
  def mismatches
  end
 
end

total_length = 0
def sam_reader_cycle(filename, &block)
  return unless block_given?
  @inputfile = File.new(filename, 'r')
  while not @inputfile.eof do
    input_string = @inputfile.gets()
    
    aa = input_string.split("\t")
    total_length = aa[2].scan(/[0-9]+/)[0].to_i if aa[0] == "@SQ"
    next if input_string[0].chr == '@'
    data = {:qname => aa[0],
            :flag  => aa[1].to_i,
            :rname => aa[2], 
            :pos   => aa[3].to_i, 
            :mapq  => aa[4].to_i,
            :cigar => aa[5],
            :rnext => aa[6],
            :pnext => aa[7].to_i,
            :tlen  => aa[8].to_i,
            :seq   => aa[9],
            :qual  => aa[10], 
            :total_length => total_length}
    a_parsed = SAMString.new(data)
    next if a_parsed.start_pos == 0
    yield a_parsed
  end
  @inputfile.close
end


