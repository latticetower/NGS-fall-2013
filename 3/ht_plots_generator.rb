# this codefile uses ruby 1.9.3 and is executed locally, with previously generated data hashes as an input.
# data hashes are generated from sam files
# example input file line:
# {:type=>"plot", :filename => 'ff.svg', :plot_data => {:values => [[0,1,1],[1,2,1]], :name => 'Insertion distance between pair reads', :x_name => 'Distance between pair reads', :y_name => 'Number of reads'}}
## ARGV - must contain data files like "bwa_data125", "bowtie_data"
## 
# - this line produces svg image
require './helpers/plotter_utility.rb'


@@data_hash = Hash.new

# process all data generated previously
ARGV.each do |filename|
  File.open(filename) do |datafile|
    while not datafile.eof? do
      inputline = eval(datafile.gets)
      if inputline.has_key?(:type) and inputline[:type] == "plot" #if data string describes plot, get plot name, and plot data and output somewhere...
        PlotterUtility::generate_plot('img/' + inputline[:filename], inputline[:plot_data])
      else
        @@data_hash.merge!(inputline)
      end
    end
  end
end

# test require
require 'mustache'


#test data

filename2 = File.new('hometask3.txt', 'w')
File.open('hometask3.txt.template') do |file|
  while not file.eof? do
    filename2.puts Mustache.render(file.gets, @@data_hash)
  end
end
filename2.close