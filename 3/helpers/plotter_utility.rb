require "rubyvis"

module PlotterUtility
  def self.generate_plot(filename, data)
    self.create_plot(to_openstruct(data[:values]), filename, :name => data[:name],
    :x_name => data[:x_name], :y_name => data[:y_name])
  #gc_distribution, ARGV[1], 
  #:name => 'Insertion distance between pair reads', 
  #:x_name => 'Distance between pair reads', 
  #:y_name => 'Number of reads')
  end
  
  # let's assume that data_values is array of arrays: [[]] each subarray contains 3 elements - x, y and z coordinate. array is unsorted
  def self.to_openstruct(data_values)
      data_values.map{|point| OpenStruct.new({:x => point[0], :y => point[1], :z => point[2] })}.sort{|d1, d2| d1.x <=> d2.x }
  end
  
  def self.to_openstruct2(data_values, text_x)
      data_values.zip(text_x).map{|point| OpenStruct.new({:x => point[0], :y => point[1], :z => point[2], :text_x =>  point[3] })}.sort{|d1, d2| d1.x <=> d2.x }
  end
  
  def self.create_plot(data_values, filename, *options)
    w = 800
    h = 700
    # /* The root panel. */
    vis1 = pv.Panel.new().width(w + 40).height(h + 40).bottom(40).left(100).right(40).top(5)
      

    x = pv.Scale.linear(data_values, lambda {|d| d.x}).range(0, w)
    y = pv.Scale.linear(data_values, lambda {|d| d.y}).range(0, h)
    c = pv.Scale.log(1, 100).range("orange", "brown")
    
    panel = vis1.add(Rubyvis::Panel).left(10).top(20).width(w).height(h)
    panel.anchor('top').add(Rubyvis::Label).top(5).left(100).text(options[0][:name])
    panel.anchor('left').add(Rubyvis::Label).left(5).top(100).text(options[0][:y_name])
    panel.anchor('bottom').add(Rubyvis::Label).right(5).bottom(5).text(options[0][:x_name])
    #/* X-axis and ticks. */
    xRule = panel.add(pv.Rule).data(x.ticks()).visible(lambda {|d| !d.nil? && d > 0}).strokeStyle(lambda {|d|  d!=0 ? "rgba(128,128,128,.2)" : "#000"}).left(x).bottom(-10).anchor("bottom")
    
    xRule.add(pv.Label).text(x.tick_format)
    
    panel.add(pv.Rule).data(x.ticks()).left(x).bottom(-10).stroke_style(lambda {|d|  d!=0 ? "rgba(255,255,255,.3)" : "#000"})
    #/* Y-axis and ticks. */
    panel.add(pv.Rule).data(y.ticks(10)).visible(lambda {|d| !d.nil? && d > 0}).bottom(y).strokeStyle(lambda {|d|  d!=0 ? "rgba(128,128,128,.2)" : "#000"}).anchor("left").add(pv.Label).text(y.tick_format)
    panel.add(pv.Rule).data(y.ticks(10)).left(0).bottom(y).stroke_style(lambda {|d|  d!=0 ? "rgba(255,255,255,.3)" : "#000"})
    
    

    panel.add(Rubyvis::Line).data(data_values).line_width(2).left(lambda {|d| x.scale(d.x)}).bottom(lambda {|d| y.scale(d.y)}).
      interpolate("linear").stroke_style("lightblue")
    panel.add(Rubyvis::Dot).data(data_values).stroke_style(lambda {|d| c.scale(d.z)}).
    fill_style(lambda {|d| c.scale(d.z).alpha(0.5)}).shape_size(lambda {|d| d.z}).
  line_width(1).left(lambda {|d| x.scale(d.x)}).bottom(lambda {|d| y.scale(d.y)}).
      interpolate("linear") #.stroke_style("green")
    vis1.render()

    @outfile1 = File.new(filename, 'w')
    @outfile1.puts vis1.to_svg
    @outfile1.close
  end
  
  def self.print_data(data)
   data.each do |d|
    puts "x: #{d.x}, y:#{d.y}"
   end
  end
  
  def self.save_to_file(data, filename)
  
    @outfile1 = File.new(filename, 'w')
    @outfile1.puts "[" + data.map{|element|  yield(element) }.join(',') + "]" if block_given? 
    @outfile1.close
  end
end