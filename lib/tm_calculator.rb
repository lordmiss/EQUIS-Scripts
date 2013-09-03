require 'bio'

class NAcid < Bio::Sequence::NA
	
	# Returns a hash with :tm, :b_tm, and :error_code.
	#
	#   seq = NAcid.new("agtcctgca")
	#   seq.calc_tm[:tm] # salt adjusted tm value
	#   seq.calc_tm[:b_tm] # basic tm value
	#   seq.calc_tm(:conc_na => 0.10) # 100 nM Na+
	#   NAcid.new("").calc_tm # { :error_code => 1 }
	#
	def calc_tm(opts = {})
		opts = { :conc_na => 0.05 }.merge(opts)
		conc_na = opts[:conc_na]
		size = self.length
		if size < 1
			return {:error_code => 1}
		end
		count = self.composition
		num_a, num_t, num_g, num_c = count['a'], count['t'], count['g'], count['c']
		
		if size < 14
			b_tm = (num_a+num_t)*2 + (num_g+num_c)*4
			tm = (num_a+num_t)*2 + (num_g+num_c)*4 - 16.6 * Math.log10(0.05) + 16.6 * Math.log10(conc_na)
		else
			b_tm = 64.9 + 41*(num_g+num_c-16.4)/size
			tm = 100.5 + (41*self.gc_content) - (820/size) + 16.6 * Math.log10(conc_na)
		end
		return {:error_code=>0, :tm=>tm.round(2), :b_tm=>b_tm.round(2)}
	end
	
	# get front half sequence
	def get_f_half
		half = (self.length) / 2
		self[0..(half - 1)]
	end
	
	# get tail half sequence
	# exclude the center NA if sequence length is odd
	def get_t_half
		half = (self.length) / 2
		if self.length % 2 == 0
			self[half..-1]
		else
			self[(half+1)..-1]
		end
	end
	
	def is_symmetric?
		f_end = self.get_f_half
		t_end = self.get_t_half
		f_end == t_end.complement
	end
	
	def get_symmetric_part
		f_end = self.get_f_half.to_s
		t_end_c = self.get_t_half.complement.to_s # complement sequence of tail half sequence
		longest_common = (f_end.split_string & t_end_c.split_string)[0]
		i = f_end.index /#{longest_common}/
		return [longest_common, i]
	end
	
	def find_seq_at_tm(tm)
		val = nil
		if tm.to_f > 12.0
			max = self.length
			max_tm = self.calc_tm[:tm]
			if max_tm > tm
				i = 3
				while self[0..i].calc_tm[:tm] < tm
					i += 1
				end
				val = self[0..i]
			end
		end
		return val
	end
	
end # class NA

class String
	def split_string
		# return all the substrings with descending size
		output = []
		(0..self.length).each do |i|
			(i..self.length).each do |j|
				output << self[i..j] unless i == j # exclude empty string, ""
			end
		end
		return output.uniq.sort_by{|x| x.size}.reverse
	end
end