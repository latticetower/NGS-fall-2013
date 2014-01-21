module GraphUtils
  @@kmers_no_by_str  # this hash for getting kmer labels fast: key - kmer string, value - node number 
  @@kmers_str_by_no  # this array for getting kmer str from label fast, by node number
  @@kmer_coverage
  @@strings
  @@edges
  @@kmer_array_size
  @@node_outdegrees
  @@node_indegrees
  @@node_outedges
  #@@node_inedges
  @@kmer_size
  
  def self.nodes
    @@kmers_no_by_str.each_pair
  end
  
  def self.clear_all  
    @@strings = []
    @@kmers_no_by_str = {}
    @@kmers_str_by_no = []
    @@kmer_coverage = []
    @@edges = {}
    @@kmer_array_size = 0
    @@node_outdegrees = []
    @@node_indegrees = []
    @@node_outedges = []
    #@@node_inedges = []
  end
  
  def self.add_nodes(node_str)
    if @@kmers_no_by_str.has_key? node_str
      @@kmer_coverage[@@kmers_no_by_str[node_str]] += 1
      return @@kmers_no_by_str[node_str]
    end
    
    @@kmers_no_by_str[node_str] = @@kmer_array_size
    @@kmers_str_by_no[@@kmer_array_size] = node_str
    @@kmer_coverage[@@kmer_array_size] = 1
    @@node_outdegrees[@@kmer_array_size] = 0
    @@node_indegrees[@@kmer_array_size] = 0
    @@node_outedges[@@kmer_array_size] = []
    @@kmer_array_size += 1
    @@kmer_array_size - 1    
  end
  
  def self.edge_exists?(first_node_no, second_node_no, label_str)
    return false unless @@edges.has_key? first_node_no
    return false if @@edges[first_node_no].nil?
    return false if @@edges[first_node_no].has_key? first_node_no
    return false if @@edges[first_node_no][first_node_no].nil?
    @@edges[first_node_no][first_node_no].include? label_str
  end
 
  def self.add_edges(first_node_no, second_node_no, label_str)
    @@edges[first_node_no] ||= {}
    @@node_outdegrees[first_node_no] ||= 0
    @@node_indegrees[second_node_no] ||= 0
    @@node_outedges[first_node_no] ||= []

    @@edges[first_node_no][second_node_no] ||= []
    unless @@edges[first_node_no][second_node_no].include?(label_str)
      @@edges[first_node_no][second_node_no] << label_str 
      @@node_outdegrees[first_node_no] += 1
      @@node_indegrees[second_node_no] += 1
      @@node_outedges[first_node_no] << second_node_no
    end
  end
  
  def self.get_edge_labels(first_node_no, second_node_no)
    return [] unless @@edges.has_key? first_node_no
    return [] unless @@edges[first_node_no].has_key? second_node_no
    @@edges[first_node_no][second_node_no] 
  end
  
  def self.check_node_degrees_and_remove_if_orphan(node)
    return if @@node_indegrees[node] > 0
    return if @@node_outdegrees[node] > 0
    # if here- means that node is orphan
    # puts "fesffse #{@@node_indegrees[node]} #{@@node_outdegrees[node]}" if node == 0
    self.remove_nodes(node)
  end
  
  def self.remove_edges(first_node_no, second_node_no)
    # puts "remove edges called: ", "#{first_node_no}, #{second_node_no}" 
    return if @@edges[first_node_no].nil?
    puts "something bad happened in remove_edges" if @@edges[first_node_no][second_node_no].size > 1
    @@edges[first_node_no][second_node_no] = nil
    @@node_outdegrees[first_node_no] -= 1
    @@node_indegrees[second_node_no] -= 1
    
    @@node_outedges[first_node_no].delete(second_node_no) unless @@node_outedges[first_node_no].nil?
    
  end
  
  def self.remove_nodes(node)
    @@edges[node] = nil
    # @@node_outdegrees[node] = nil
    # @@node_indegrees[node] = nil
    @@node_outedges[node] = nil
    kmer_str = @@kmers_str_by_no[node]
    # @@kmers_no_by_str.delete(kmer_str)
    @@kmers_str_by_no[node] = nil
    # @@kmer_coverage[@@kmer_array_size] = 1
  end
  @@compl_pairs = {'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', 'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c'}
  def self.get_complement(str)
    str.reverse.chars.map{|x| @@compl_pairs[x] }.join('')
  end
  
  def self.read_fasta(filename)
    self.clear_all
    
    referencefile = File.new(filename)
    genome_string = ""
    while not referencefile.eof do
      str = referencefile.gets().strip
      if str.index('>').nil?
        genome_string << str 
      else
        if genome_string.size > 0
          yield(genome_string) if block_given?
          @@strings << genome_string
          genome_string = ""
        end
      end
    end
    if genome_string.size > 0
      yield(genome_string) if block_given?
      @@strings << genome_string 
    end
    referencefile.close
  end
  
  def self.read_fastq(filename)
    self.clear_all
    
    referencefile = File.new(filename)
    while not referencefile.eof do
      str = referencefile.gets().strip
      str = referencefile.gets().strip
      yield(str) if block_given?
      @@strings << str
      str = referencefile.gets().strip
      str = referencefile.gets().strip
    end
    referencefile.close
  end
  
  def self.set_kmer_size(k)
    @@kmer_size = k
    # puts "#{@@kmer_size }"
  end
  
  def self.build_kmers(read)
    k = @@kmer_size 
    # puts "#{@@kmer_size }"
    first_node = read.slice(0, k)
    first_node_no = self.add_nodes(first_node)
    
    (read.size - k).times do |i|
      edge_label = first_node + read[k + i].chr
      second_node = edge_label.slice(1, k)
      second_node_no = self.add_nodes(second_node)
      unless self.edge_exists?(first_node_no, second_node_no, edge_label)
        self.add_edges(first_node_no, second_node_no, edge_label)
      end
      
      first_node = second_node
      first_node_no = second_node_no
    end
  end
  
  def self.write_to_fasta(outputfile)
    file = File.new(outputfile, 'w')
    file.puts ">-------------------nodes-----------------------------"
    puts @@kmers_str_by_no[0]
    @@kmer_array_size.times do |i|
      next if @@kmers_str_by_no[i].nil? or @@kmers_str_by_no[i].size < 1
      file.puts ">node number #{i} (#{@@node_outdegrees[i]} outgoing edges, #{@@node_indegrees[i]} ingoing edges):", "#{@@kmers_str_by_no[i]}" 
    end
    file.puts ">---------------------edges:-----------------------"
    @@kmer_array_size.times do |i|
      next if @@edges[i].nil?
      @@kmer_array_size.times do |j|
        next if @@edges[i][j].nil? or @@edges[i][j].size < 1
        @@edges[i][j].each do |str|
          file.puts ">edge from node #{i} to node #{j}:", 
                    "#{str}"
        end
      end
    end
    file.close
  end
  
  def self.get_edge_weight(edge)
    coverage = 0.0
    
    (edge.size - @@kmer_size + 1).times do |i|
      kmer_no = @@kmers_no_by_str[edge.slice(i, @@kmer_size)]
      kmer_no ||= 0
      coverage += @@kmer_coverage[kmer_no] 
    end    
    (100*coverage/(edge.size - @@kmer_size + 1)).round/100.0
  end
  
  def self.write_to_console
    puts "-------------------", "nodes", "-----------------------------"
    @@kmer_array_size.times do |i|
      next if @@kmers_str_by_no[i].nil? or @@kmers_str_by_no[i].size < 1
      puts "#{i}: #{@@kmers_str_by_no[i]}"
    end
    puts "---------------------", "edges:", "-----------------------"
    @@kmer_array_size.times do |i|
      next if @@edges[i].nil?
      @@kmer_array_size.times do |j|
        next if @@edges[i][j].nil? or @@edges[i][j].size < 1
        puts ">edge from #{i} to #{j}:", 
             "#{@@edges[i][j]}"
      end
    end
  end
  
  def self.simplify!
    changes_flag = true
    # begin
    while changes_flag do 
      changes_flag = false
      @@kmer_array_size.times do |index1|
        next if @@edges[index1].nil?
        @@kmer_array_size.times do |index2|
          next if @@edges[index1][index2].nil? or @@edges[index1][index2].size < 1
          # puts "#{index2} #{@@node_outdegrees[index2]}"
          if @@node_outdegrees[index2] == 1 && @@node_indegrees[index2] == 1
            k = @@node_outedges[index2][0]
            #next if @@edges[index2][k].nil? or @@edges[index2][k].size < 1
            label1s = self.get_edge_labels(index1, index2)
            label2s = self.get_edge_labels(index2, k)
            if label1s.size != 1 or label2s.size!=1
              puts "some error while simplifying"
              return
            end
            label1 = label1s[0]
            label2 = label2s[0]
            # puts label1
            new_edge_label = label1.slice(0, label1.size - @@kmer_size + 1) + label2
            # puts new_edge_label
            self.remove_edges(index1, index2)
            self.remove_edges(index2, k)
            unless self.edge_exists?(index1, k, new_edge_label)
              self.add_edges(index1, k, new_edge_label)
            end
        
            self.check_node_degrees_and_remove_if_orphan(index2)
            changes_flag = true
          end
        end
      end
      #end of kmer array cycle
    end
    # end
  end
  
  def self.write_to_dot(outfile)
    file = File.new(outfile, 'w')
    file.puts "digraph g {"
    #begin
    @@kmer_array_size.times do |i|
      next if @@edges[i].nil?
      @@kmer_array_size.times do |j|
        next if @@edges[i][j].nil? or @@edges[i][j].size < 1
        @@edges[i][j].each do |lbl|
          file.puts "#{i} -> #{j} [label=\"#{lbl}\"];"
        end
      end
    end
    #end
    file.puts "}"
    file.close
  end
  
  def self.write_to_fasta_and_dot(outfile_fasta, outfile_dot)
    file_fasta = File.new(outfile_fasta, 'w')
    # file_fasta.puts ">-------------------nodes-----------------------------"
    @@kmer_array_size.times do |i|
      next if @@kmers_str_by_no[i].nil? or @@kmers_str_by_no[i].size < 1
      file_fasta.puts ">node number #{i} (#{@@node_outdegrees[i]} outgoing edges, #{@@node_indegrees[i]} ingoing edges):", "#{@@kmers_str_by_no[i]}" 
    end
    # file_fasta.puts ">---------------------edges:-----------------------"
    
    file_dot = File.new(outfile_dot, 'w')
    file_dot.puts "digraph g {"
    
    edges_iterator = 1 #renumerates edges
    
    @@kmer_array_size.times do |i|
      next if @@edges[i].nil?
      @@kmer_array_size.times do |j|
        next if @@edges[i][j].nil? or @@edges[i][j].size < 1
        @@edges[i][j].each do |str|
          file_fasta.puts ">edge N#{edges_iterator}, from node #{i} to node #{j}:"
          splitter = str.size/78
          (splitter + 1).times do |k| 
            file_fasta.puts str.slice(78*k, [78, str.size - 78*k].min)
          end
          file_dot.puts "#{i} -> #{j} [label=\"N#{edges_iterator}, w=#{self.get_edge_weight(str)}\"];"
          
          edges_iterator += 1
        end
      end
    end
    file_dot.puts "}"
    file_fasta.close
    file_dot.close
  end
  
  def self.print_to_file  
    f = File.new('debug.txt', 'w')
    f.puts "start of graph data"
    f.puts @@kmer_size
    f.puts @@kmers_no_by_str
        # this hash for getting kmer labels fast: key - kmer string, value - node number 
    # this array for getting kmer str from label fast, by node number
    f.puts @@kmer_coverage.join(', ')
    f.puts @@strings.size
    
    f.puts "[ " + @@edges.map{|x| x.nil? ? "[]": "[" + x.map{|y|  (y.nil? ? "0" : y.to_s )}.join(', ') + "]" }.join(', ') + " ]"
    @@kmer_array_size.times do |i|
      f.puts "#{i}: #{@@kmers_str_by_no[i]}, [ #{(@@node_outedges[i] || []).join(', ')} ] #{@@node_indegrees[i]} - #{@@node_outdegrees[i]}"
    end
    @@kmer_array_size.times do |i|
      if @@edges[i].nil? && !(@@node_outedges[i].nil? || @@node_outedges[i].size == 0)
          puts "something wrong with #{i} #{@@node_outedges[i].join ','}"
          next
      end
      # f.puts "#{@@kmers_str_by_no[i]}: #{@@node_indegrees[i]},  #{@@node_outdegrees[i]}, [#{@@node_outedges[i].join(',')}] - outedges"
      @@kmer_array_size.times do |j| 
        
        if !@@edges[i].nil? && @@edges[i][j].nil? && @@node_outedges[i].include?(j)
          puts "something wrong with edge #{i} #{j} #{@@node_outedges[i]}"
        end
        if !@@edges[i].nil? && (!@@edges[i][j].nil?) && (!@@node_outedges[i].include?(j))
          puts "something wrong2 with edge #{i} #{j}"
        end
      end
    end

  
  
  end
  
end